library(DiffBind)
library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(KEGG.db)
library(clusterProfiler)
library(fgsea)
data("TSS.mouse.GRCm38")

setwd("C:/path/to/dir")

# Create target file
SampleID <- c("ko1", "ko2", "ko3", "wt1", "wt2", "wt3")
Tissue <- rep("unknown", 6) 
Factor <- rep("Plagl1", 6)
Condition <- rep("unknown", 6)
Treatment <- c("ko", "ko", "ko", "wt", "wt", "wt")
Replicate <- c(1,2,3,1,2,3)
bamReads <- c("no_mt_bam/ko1.no_mt.bam", "no_mt_bam/ko2.no_mt.bam",
              "no_mt_bam/ko3.no_mt.bam", "no_mt_bam/wt1.no_mt.bam",
              "no_mt_bam/wt2.no_mt.bam", "no_mt_bam/wt3.no_mt.bam")

Peaks <- c("broad_peaks/ko1/ko1_peaks.broadPeak", "broad_peaks/ko2/ko2_peaks.broadPeak",
           "broad_peaks/ko3/ko3_peaks.broadPeak", "broad_peaks/wt1/wt1_peaks.broadPeak",
           "broad_peaks/wt2/wt2_peaks.broadPeak", "broad_peaks/wt3/wt3_peaks.broadPeak")
PeakCaller <- rep("narrow", 6)

samples <- data.frame(SampleID=SampleID, Tissue=Tissue, Factor=Factor,
                      Condition=Condition, Treatment=Treatment, Replicate=Replicate,
                      bamReads=bamReads, Peaks=Peaks, PeakCaller=PeakCaller)
samples
write.csv(samples, file="samples_broad.csv")

# Read in hotspot data and create dba object
atac <- dba(sampleSheet=samples)
save(atac, file="atac_broad_DBA.RData")
tiff("Correlation_heatmap_using_occupancy_peakCaller_data_broad.tiff", 
     width = 800, height = 800)
plot(atac)
dev.off()

# Obtain read counts for hotspots
atac <- dba.count(atac, minOverlap=2, bRemoveDuplicates=T,
                  summits=F, bUseSummarizeOverlaps=T)
atac
save(atac, file="ATAC_diffind_counts_broad.DBA.RData")

# Overlap rate
dba.overlap(atac, atac$masks$ko, mode=DBA_OLAP_RATE)
# [1] 62966 41558 29703
dba.overlap(atac, atac$masks$wt, mode=DBA_OLAP_RATE)
# 67767 48705 37402

# Show the number of reads overlapping peaks
info <- dba.show(atac)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, 
                  PeakReads=round(info$Reads * info$FRiP))
libsizes
write.csv(libsizes, file="FRiP_broad.csv")

tiff("Correlation_heatmap_using_counts_data_broad.tiff", width = 800, height = 800)
plot(atac)
dev.off()

# Plot PCA
tiff("ATAC-seq_PCA_based_on_read_counts_broad.tiff", width = 800, height = 800)
dba.plotPCA(atac, DBA_TREATMENT, label=DBA_ID)
dev.off()

## Specify contrast
atac <- dba.contrast(atac, contrast=c("Treatment", "ko", "wt"))
dba.show(atac, bContrasts = T)

## ---- Native normalization with library size based on reads in background ----- ##
## Normalize atac seq (all methods, library sizes based on background)
load("diffbind_objects/ATAC_diffind_counts_broad.DBA.RData")
atac <- dba.normalize(atac, method=DBA_EDGER, normalize=DBA_NORM_NATIVE,
                      background = T, library=DBA_LIBSIZE_BACKGROUND)

## Compare groups
atac <- dba.analyze(atac, method=DBA_EDGER)
dba.show(atac, bContrasts = T)
save(atac, file="atac_broad_analyzed.RData")

## MA plot 
tiff("native_normalization_background_MA.tiff", width = 500, height = 500)
par(mfrow=c(1,1))
dba.plotMA(atac, method = DBA_EDGER, sub="edgeR:TMM:background", th=0.05)
dev.off()

## Volcano plot 
tiff("native_normalization_background_volcano.tiff", width = 500, height = 500)
dba.plotVolcano(atac, contrast=1, method=DBA_EDGER, th=0.05, bUsePval=F)
dev.off()

## Heatmap plot 
tiff("native_normalization_background_heatmap.tiff", width = 500, height = 500)
dba.plotHeatmap(atac, contrast=1, maxSites = 833, method=DBA_EDGER, th=0.05,
                bUsePval = F, distMethod = "pearson")
dev.off()

## Get the report, convert to GRanges object
atac.DB <- dba.report(atac, method = DBA_EDGER, contrast = 1, th = 1,
                      bUsePval = F)
atac.DB
save(atac.DB, 
     file=paste("KO", "_vs_", "WT", ".GRanges.RData", sep=""))

# annotate by distance to TSS
atac.DB_annot <- annotatePeakInBatch(atac.DB, 
                                     AnnotationData=TSS.mouse.GRCm38,
                                     output="nearestLocation",
                                     PeakLocForDistance = "middle",
                                     featureType = "TSS")
dim(as.data.frame(atac.DB_annot) %>% filter(FDR < 0.05))
dim(as.data.frame(atac.DB_annot) %>% filter(FDR < 0.05) %>% 
      filter(Fold < 0))
# 814 are down in KO

dim(as.data.frame(atac.DB_annot) %>% filter(FDR < 0.05) %>% 
      filter(Fold > 0))
# 20 are up in KO

# Add gene additional gene info
atac.DB_annot <- addGeneIDs(atac.DB_annot, orgAnn="org.Mm.eg.db", 
                                  IDs2Add=c("symbol", "entrez_id"))
write.csv(as.data.frame(atac.DB_annot), 
          file=paste("KO", "_vs_", "WT", ".diff_bind.csv", sep=""), row.names=F)


## Prepare ID lists (all, up, and down) for GeneSCF analysis
all_5000 <- atac.DB_annot[as.data.frame(atac.DB_annot)$FDR < 0.05,]
as.data.frame(all_5000)$distancetoFeature
all_5000 <- all_5000[abs(as.data.frame(all_5000)$distancetoFeature) <= 5000,]
all_5000
write.csv(all_5000, file="all_sig_5000bp.csv", row.names = F)

## Distribution of aggregated peak scores around TSS
all_sig <- atac.DB_annot[as.data.frame(atac.DB_annot)$FDR < 0.05,]
tiff("distr_sig_peak_numbers_tss.tiff", width = 500, height=500)
binOverFeature(all_sig, annotationData=TSS.mouse.GRCm38,
               radius=5000, nbins=10, FUN=length,
               ylab="count",
               main=c("Distribution of aggregated peak numbers around TSS"))
dev.off()

# distribution over genomic features
if(require(TxDb.Mmusculus.UCSC.mm10.knownGene)){
  aCR<-assignChromosomeRegion(all_sig, nucleotideLevel=FALSE, 
                              precedence=c("Promoters", "immediateDownstream", 
                                           "fiveUTRs", "threeUTRs", 
                                           "Exons", "Introns"), 
                              TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene)
  acR <- as.data.frame(aCR$percentage)
  names(acR) <- c("Genomic_feature", "Percentage")
  ggplot(acR, aes(x=reorder(Genomic_feature, -Percentage), y=Percentage)) +
    geom_bar(stat="identity") + theme(axis.text.x = element_text(angle=90)) +
    xlab("Genomic feature") + theme(axis.text = element_text(size=12, face="bold")) +
    theme(axis.title = element_text(size=14)) + ylab("Frequency (%)") +
    ggsave("DARS_distr_genomic_features.tiff", dpi=200)
}


## Run gene set enrichment analysis using clusterProfiler
## Load the results for diffbind analysis
res <- read.csv("diffbind_comparison_broad/KO_vs_WT.diff_bind.csv", header = T)
# limit to the genes located no further than 5000 bp from TSS
res <- res[abs(res$distancetoFeature) <= 5000,]
dim(res)
head(res)
res <- res[!is.na(res$entrez_id),]
head(res$Fold)
res <- res[!duplicated(res$entrez_id),]
gsea_in <- res$Fold 
names(gsea_in) <- res$entrez_id
head(gsea_in)
gsea_in <- sort(gsea_in, decreasing = T)
head(gsea_in)
head(names(gsea_in))
grep(";", test)
names(gsea_in) <- gsub(";\\d+", "", names(gsea_in))

# Run GSEA against GO terms
gse <- gseGO(geneList = gsea_in, ont="BP", keyType = "ENTREZID",
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, 
             pAdjustMethod = "BH", by = "fgsea",
             OrgDb = org.Mm.eg.db)
head(gse@result)
write.csv(gse@result, file="GSEA_BP_FC.csv", row.names = F)

# Run GSEA against GO terms (CC)
gse <- gseGO(geneList = gsea_in, ont="CC", keyType = "ENTREZID",
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, 
             pAdjustMethod = "BH", by = "fgsea",
             OrgDb = org.Mm.eg.db)
head(gse@result)
write.csv(gse@result, file="GSEA_CC_FC.csv", row.names = F)

# Run GSEA against GO terms (MF)
gse <- gseGO(geneList = gsea_in, ont="MF", keyType = "ENTREZID",
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, 
             pAdjustMethod = "BH", by = "fgsea",
             OrgDb = org.Mm.eg.db)
head(gse@result)
write.csv(gse@result, file="GSEA_MF_FC.csv", row.names = F)

# Run GSEA against KEGG pathways
kk <- gseKEGG(geneList = gsea_in, organism = "mmu", nPerm = 1000,
              keyType = "kegg", minGSSize = 10, maxGSSize = 500,
              pAdjustMethod = "BH")
write.csv(kk@result, file="GSEA_KEGG_FC.csv", row.names = F)


# Compare ATAC-seq DARs and differentially expressed genes
# First detect differentially expressed genes
library("dplyr")
setwd("C:/Quarantine_projects/Carol_Schuurmans_projects/Carol_Schuurmans_May2019/")

list.files("comparisons_pooled/")
rna_res <- read.csv("comparisons_pooled/results_pooled_groups.txt", header=T, sep="\t",
                    fill=T, quote="\"") 
head(rna_res)
rna_sig <- rna_res %>% filter(padj < 0.05)
rna_sig <- rna_sig[,c(1, 3)]
# change the sign of log2 fold changes
rna_sig$log2FoldChange <- rna_sig$log2FoldChange * -1
head(rna_sig)
rna_up <- rna_sig %>% filter(log2FoldChange > 0)
rna_down <- rna_sig %>% filter(log2FoldChange < 0)
dim(rna_sig)
dim(rna_down)
dim(rna_up)

# Get up and down DARs
setwd("C:/Quarantine_projects/Carol_ATACseq_JAN2021/diffbind_comparison_broad/")
list.files()
dars <- read.csv("KO_vs_WT.diff_bind.csv", header=T)
head(dars)
# limit the set to significant DARs (FDR < 0.05) located no further than 5000 bp
# from the nearest TSS
dars <- dars %>% filter(FDR < 0.05) %>% filter(abs(distancetoFeature) <= 5000)
dim(dars)
# 311 significant DARs 

# Get up- and down-regulated DARs
dars_up <- dars %>% filter(Fold > 0)
dars_down <- dars %>% filter(Fold < 0)
dim(dars_up)
dim(dars_down)

intersect(dars_up$feature, rna_up$Row.names)
intersect(dars_up$feature, rna_down$Row.names)
reg <- intersect(dars_down$feature, rna_down$Row.names)
reg
opposite <- intersect(dars_down$feature, rna_up$Row.names)
opposite

# annotate overlaps
head(as.data.frame(reg))
reg <- as.data.frame(reg)
names(reg) <- "ensembl_gene_id"
reg_df <- merge(reg, annot, by="ensembl_gene_id")
write.csv(reg_df, file="same_direction.csv")
dim(reg_df)

head(as.data.frame(opposite))
opposite <- as.data.frame(opposite)
names(opposite) <- "ensembl_gene_id"
head(opposite)
opposite_df <- merge(opposite, annot, by="ensembl_gene_id")
opposite_df
write.csv(opposite_df, file="opposite_direction.csv")
dim(opposite_df)


