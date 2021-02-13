---
title: "ATAC-seq_workflow_2020"
author: "Slava"
date: "2/12/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## ATAC-seq workflow template
The workflow starts with raw reads in fastq format and ends with the detection of differentially bound peaks between experimental groups. The workflow can be separated into 4 major modules: 

1. Initial quality control and adapter trimming
2. Read alignment 
3. Peak calling
4. Detection of differentialy bound intervals

### 1. QC and trimming
1) Initial quality control is done using FastQC https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

```
fastqc <fastq>

```

This step can be parallelized with GNU parallel.

```
parallel "fastqc {}" ::: *.fastq

```

2) Trim the adapters and low quality sequences using Trim Galore! 

```
trim_galore --fastqc -q 30 --nextera --length 20 -j 4 --paired <R1.fastq> <R2.fastq>
```

3) Run fastqc again on trimmed reads and check the results.

### 2. Read alignment.

Bowtie2 and BWA are good choices for read mapping.

1) Map reads with bowtie2, then convert resulting SAM files to BAM, sort and index with samtools.

```
bowtie2 --very-sensitive -p 20 -x <Bowtie2Index> -1 <R1.fastq> -2 <R2.fastq> > <SAM>
samtools view -bS <SAM> > <BAM>
samtools sort <BAM> <sorted>
samtools index <sorted.bam>
```
2) Remove mitochondrial reads from BAM files.

```
samtools idxstats <input.bam> | cut -f 1 | grep -v MT | xargs samtools view -b <input.bam> > <no_mt.bam>
samtools index <no_mt.bam>
```

3) Collect mapping statics with samtools flagstat parsed with this bash script.

```
#! /bin/bash

# Parse and summarize the output for samtools flagstat command
# for paired-end reads

echo -e "SampleID\ttotal_reads\tmapped_reads\tproperly_paired"

for file in <BAM_DIR>/*sorted.bam
do
    filename=`basename $file`
    samplename=${filename%%.*}

    stats[0]=$samplename
      
    samtools flagstat $file > stats.tmp 
    while IFS= read -r line
    do
        if [[ $line =~ "in total"  ]]
        then
             total_reads=`echo $line | cut -f 1 -d ' '`
             stats[1]="$total_reads"
        fi

        if [[ $line =~ 'mapped (' ]]
        then
             mapped=`echo $line | cut -f 1 -d ' '`
             stats[2]="$mapped"
        fi

        if [[ $line =~ "properly paired" ]]
        then
             properly_paired=`echo $line | cut -f 1 -d ' '`
             stats[3]="$properly_paired"
        fi
    done < "stats.tmp"

    echo -e "${stats[0]}\t${stats[1]}\t${stats[2]}\t${stats[3]}"

done

rm stats.tmp
```
4. Check fragment length distribution and make sure that 2 major peaks are present, one at < 100 bp representing areas of open chromatin and 200 bp peak representing mono-nucleosomes; thare may be other smaller peaks at 400, 600, and 800 bp due to the presence of multi-nucleosomes.

Use this one liner to collect the data and an R script to build a plot.

```
samtools view <BAM> | awk '$9>0' | cut -f 9 | sort | uniq -c | sed -e 's/^[ \t]*//' > <fragment_length.txt>
```

```
library(ggplot2)
library(dplyr)

setwd("<DIR>")
list.files()

# Get fragment length count data
plotFragmentLength <- function(file_name, prefix) 
{
  res <- read.table(file_name, sep=" ", header = F, stringsAsFactors = F)
  names(res) <- c("count", "length")
  
  res <- res %>% filter(length < 1200)
  
  # Plot fragment length count vs fragment length
  ggplot(res, aes(y=count, x=length)) + geom_line() + 
    ggtitle(prefix)
  ggsave( paste(prefix, "_length_distr.tiff", sep=""), device = "tiff", dpi=300)
}

plotFragmentLength(file_name = "<fragment_length.txt>", prefix = "<prefix>")

```

### 3. Peak calling.

Call peaks with MACS2 https://github.com/macs3-project/MACS.

```
macs2 callpeak -t <BAM> -f BAMPE -n <sample_name> -g mm -q 0.05 --keep-dup 1 --outdir <out_dir>
```
-BAMPE - force MACS to use full length fragments taking the advantage of paired-end nature of the data.

-g - effective genome siz, in this example I used mouse. Effective genome sizes can be found here: https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html.

-q - FDR adjusted p-value cutoff for peaks.

--keep-dup - what to do with duplicate reads? *1* - ignore duplicate reads. 

*NOTE*: In the future consider using Genrich https://github.com/jsh58/Genrich.

### 4. Detect differentially bound regions using Bioconductor package DiffBind.

```
library(DiffBind)
library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(KEGG.db)
library(reactome.db)
data("TSS.mouse.GRCm38")

## This example is based on Mouse GRCm38 genome and the experiment the includes 2 grops KO and WT with 3 replicates in each.

setwd("<DIR>")

# Create target file
SampleID <- c("ko1", "ko2", "ko3", "wt1", "wt2", "wt3")
Tissue <- rep("SomeTissue", 6) 
Factor <- rep("SomeFactor", 6)
Condition <- rep("SomeCondition", 6)
Treatment <- c("ko", "ko", "ko", "wt", "wt", "wt")
Replicate <- c(1,2,3,1,2,3)
bamReads <- c("no_mt_bam/ko1.no_mt.bam", "no_mt_bam/ko2.no_mt.bam",
              "no_mt_bam/ko3.no_mt.bam", "no_mt_bam/wt1.no_mt.bam",
              "no_mt_bam/wt2.no_mt.bam", "no_mt_bam/wt3.no_mt.bam")

Peaks <- c("peaks/ko1/ko1_peaks.narrowPeak", "peaks/ko2/ko2_peaks.narrowPeak",
           "peaks/ko3/ko3_peaks.narrowPeak", "peaks/wt1/wt1_peaks.narrowPeak",
           "peaks/wt2/wt2_peaks.narrowPeak", "peaks/wt3/wt3_peaks.narrowPeak")
PeakCaller <- rep("narrow", 6)

samples <- data.frame(SampleID=SampleID, Tissue=Tissue, Factor=Factor,
                      Condition=Condition, Treatment=Treatment, Replicate=Replicate,
                      bamReads=bamReads, Peaks=Peaks, PeakCaller=PeakCaller)
write.csv(samples, file="samples.csv")

# Read in hotspot data and create dba object
atac <- dba(sampleSheet=samples)
save(atac, file="atac_DBA.RData")
tiff("Correlation_heatmap_using_occupancy_peakCaller_data.tiff", 
     width = 800, height = 800)
plot(atac)
dev.off()


# Obtain read counts for hotspots
atac <- dba.count(atac, minOverlap=2)
save(atac, file="ATAC_diffind_counts.DBA.RData")

tiff("Correlation_heatmap_using_counts_data.tiff", width = 800, height = 800)
plot(atac)
dev.off()

# Show the number of reads overlapping peaks
info <- dba.show(atac)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, 
                  PeakReads=round(info$Reads * info$FRiP))
write.csv(libsizes, file="FRiP.csv")

# Plot PCA
tiff("ATAC-seq_PCA_based_on_read_counts.tiff", width = 800, height = 800)
dba.plotPCA(atac, DBA_TREATMENT, label=DBA_ID)
dev.off()

# Get masks for specific experimental groups
Mutmask <- dba.mask(atac, DBA_TREATMENT, "ko")
WTmask <- dba.mask(atac, DBA_TREATMENT, "wt")

dba.object <- atac
group1 <- Mutmask
group2 <- WTmask
name1 <- "ko"
name2 <- "wt"

# Function to detect differential binding events, annotate the results by distance to
# promoters. Calculate GO term and pathway enrichment for those differential binding spots
# located no further then 5000 bp up- or down-stream to promoters
detectDiffBind <- function(dba.object, group1, group2, name1, name2) 
{
# Calculate differential binding between two selected cell types
dba.object <- dba.contrast(dba.object, group1=group1, 
group2=group2, 
name1=name1, name2=name2, minMembers=2)
dba.object <- dba.analyze(dba.object, bParallel = F, method = "DBA_EDGER")
dba.object.DB <- dba.report(dba.object, th=1)
save(dba.object.DB, 
file=paste(name1, "_vs_", name2, ".GRanges.RData", sep=""))

# save MA plot
png(paste(name1, "_vs_", name2, ".MAplot.png", sep=""))
dba.plotMA(dba.object)
dev.off()

# Save XY plot
png(paste(name1, "_vs_", name2, ".XYplot.png", sep=""))
dba.plotMA(dba.object, bXY=TRUE)
dev.off()

# Save PCA plot
png(paste(name1, "_vs_", name2, ".PCAplot.png", sep=""))
dba.plotPCA(dba.object, contrast=1,th=.05,label=DBA_TREATMENT)
dev.off()

# annotate by promoters
promoterData <- promoters(TSS.mouse.GRCm38, upstream=5000, downstream=500)
dba.object.DB_annot <- annotatePeakInBatch(dba.object.DB, 
AnnotationData=promoterData)
# Add gene additional gene info
dba.object.DB_annot <- addGeneIDs(dba.object.DB_annot, orgAnn="org.Mm.eg.db", IDs2Add=c("symbol", "entrez_id"))
write.csv(as.data.frame(dba.object.DB_annot), file=paste(name1, "_vs_", name2, 
".diff_bind.csv", sep=""))

# Perform GO and pathway enrichment analysis 
# This will be done on differential hotspots located no further then 5000 bp from the middle of the promoter
# In addition we will do separate analysis for up- and down- regulated hotspots
all_5000 <- dba.object.DB_annot[as.data.frame(dba.object.DB_annot)$shortestDistance <= 5000,]
overGO <- getEnrichedGO(all_5000, orgAnn="org.Mm.eg.db", 
maxP=0.05, minGOterm=10, 
multiAdjMethod="BH")
overPATH <- getEnrichedPATH(all_5000, orgAnn="org.Mm.eg.db", pathAnn="reactome.db",
feature_id_type = "ensembl_gene_id", maxP=0.05, minPATHterm = 10)
# Write out the results
write.table(overGO[["bp"]], file=paste(name1, "_vs_", name2, ".bidirectional_GO_Pval0.05.txt", sep=""),
sep="\t", col.names=T, row.names=F, quote=F)
write.table(overPATH, file=paste(name1, "_vs_", name2, ".bidirectional_Path_Pval0.05.txt", sep=""),
sep="\t", col.names=T, row.names=F, quote=F)

# Up-regulated in group2
up <- dba.object.DB_annot[as.data.frame(dba.object.DB_annot)$shortestDistance <= 5000 & 
as.data.frame(dba.object.DB_annot)$Fold < 0,]
overGO <- getEnrichedGO(up, orgAnn="org.Mm.eg.db", 
maxP=0.05, minGOterm=10, 
multiAdjMethod="BH")
overPATH <- getEnrichedPATH(up, orgAnn="org.Mm.eg.db", pathAnn="reactome.db",
feature_id_type = "ensembl_gene_id", maxP=0.05, minPATHterm = 10)
# Write out the results
write.table(overGO[["bp"]], file=paste(name1, "_vs_", name2, "_up_in_", group2 ,"_GO_Pval0.05.txt", sep=""),
sep="\t", col.names=T, row.names=F, quote=F)
write.table(overPATH, file=paste(name1, "_vs_", name2, "_up_in_", "_PATH_Pval0.05.txt", sep=""),
sep="\t", col.names=T, row.names=F, quote=F)

# Down-regulated in group2
down <- dba.object.DB_annot[as.data.frame(dba.object.DB_annot)$shortestDistance <= 5000 & 
as.data.frame(dba.object.DB_annot)$Fold > 0,]
overGO <- getEnrichedGO(down, orgAnn="org.Mm.eg.db", 
maxP=0.05, minGOterm=10, 
multiAdjMethod="BH")
overPATH <- getEnrichedPATH(down, orgAnn="org.Mm.eg.db", pathAnn="reactome.db",
feature_id_type = "ensembl_gene_id", maxP=0.05, minPATHterm = 10)
# Write out the results
write.table(overGO[["bp"]], file=paste(name1, "_vs_", name2, "_down_in_", group2 ,"_GO_Pval0.05.txt", sep=""),
sep="\t", col.names=T, row.names=F, quote=F)
write.table(overPATH, file=paste(name1, "_vs_", name2, "_down_in_", group2 ,"_PATH_Pval0.05.txt", sep=""),
sep="\t", col.names=T, row.names=F, quote=F)
}

```











