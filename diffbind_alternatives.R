library(DiffBind)
library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(KEGG.db)
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

## ------------- Native normalization with full library size ---------------- ##
## Normalize atac seq (all methods, native normalization, full libraries)
atac <- dba.normalize(atac, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE,
                      background = F, library=DBA_LIBSIZE_FULL)

## Compare groups
atac <- dba.analyze(atac, method=DBA_ALL_METHODS)
dba.show(atac, bContrasts = T)
# Factor Group Samples Group2 Samples2 DB.edgeR DB.DESeq2
# Treatment    ko       3     wt        3       12         3

## MA plot before and after normalization 
tiff("native_normalization_full_library_size.tiff", width = 500, height = 500)
par(mfrow=c(2,2))
dba.plotMA(atac, method = DBA_EDGER, sub="edgeR:TMM:total library", th=0.05)
dba.plotMA(atac, method = DBA_DESEQ2, sub="DESeq2:RLE:total library", th=0.05)
dba.plotMA(atac, bNormalized = F, sub="Non-normalized", th=0)
dev.off()

## Plot venn diagrams for results
tiff("native_normalization_full_library_edger_deseq2.tiff", width=500, height=500)
dba.plotVenn(atac, contrast = 1, method=DBA_ALL_METHODS, bDB = T)
dev.off()

## ------------------ Loess normalization with background ------------------ ##
## Loess normalization
atac <- dba.normalize(atac, offsets = T, background = T, 
                      library=DBA_LIBSIZE_BACKGROUND, method=DBA_ALL_METHODS)
atac <- dba.analyze(atac)
dba.show(atac, bContrasts = T)
# Factor Group Samples Group2 Samples2 DB.edgeR DB.DESeq2
# 1 Treatment    ko       3     wt        3       12         8

## MA plot before and after normalization
tiff("loess_normalization_background.tiff", width = 500, height = 500)
par(mfrow=c(2,2))
dba.plotMA(atac, method = DBA_EDGER, sub="edgeR:loess", th=0.05)
dba.plotMA(atac, method = DBA_DESEQ2, sub="DESeq2:loess", th=0.05)
dba.plotMA(atac, bNormalized = F, sub="Non-normalized", th=0)
dev.off()

## Plot venn diagrams for results
tiff("loess_normalization_venn_edger_deseq2.tiff", width=500, height=500)
dba.plotVenn(atac, contrast = 1, method=DBA_ALL_METHODS, bDB = T)
dev.off()

## ---- Native normalization with library size based on reads in peaks----- ##
## Normalize atac seq (all methods, native normalization, reads in peaks)
atac <- dba.normalize(atac, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE,
                      background = F, library=DBA_LIBSIZE_PEAKREADS)

## Compare groups
atac <- dba.analyze(atac, method=DBA_ALL_METHODS)
dba.show(atac, bContrasts = T)
# Factor Group Samples Group2 Samples2 DB.edgeR DB.DESeq2
# 1 Treatment    ko       3     wt        3       12         3

## MA plot before and after normalization
tiff("native_normalization_reads_in_peaks.tiff", width = 500, height = 500)
par(mfrow=c(2,2))
dba.plotMA(atac, method = DBA_EDGER, sub="edgeR:reads in peaks", th=0.05)
dba.plotMA(atac, method = DBA_DESEQ2, sub="DESeq2:reads in peaks", th=0.05)
dba.plotMA(atac, bNormalized = F, sub="Non-normalized", th=0)
dev.off()

## Plot venn diagrams for results
tiff("native_normalization_reads_in_peaks_venn_edger_deseq2.tiff", width=500, height=500)
dba.plotVenn(atac, contrast = 1, method=DBA_ALL_METHODS, bDB = T)
dev.off()

## ---- Native normalization with library size based on reads in background ----- ##
## Normalize atac seq (all methods, library sizes based on background)
atac <- dba.normalize(atac, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE,
                      background = T, library=DBA_LIBSIZE_BACKGROUND)

## Compare groups
atac <- dba.analyze(atac, method=DBA_ALL_METHODS)
dba.show(atac, bContrasts = T)
# Factor Group Samples Group2 Samples2 DB.edgeR DB.DESeq2
# 1 Treatment    ko       3     wt        3      629       650

## MA plot before and after normalization 
tiff("native_normalization_background.tiff", width = 500, height = 500)
par(mfrow=c(2,2))
dba.plotMA(atac, method = DBA_EDGER, sub="edgeR:TMM:background", th=0.05)
dba.plotMA(atac, method = DBA_DESEQ2, sub="DESeq2:RLE:background", th=0.05)
dba.plotMA(atac, bNormalized = F, sub="Non-normalized", th=0)
dev.off()

## Plot venn diagrams for results
tiff("native_normalization_background_venn_edger_deseq2.tiff", width=500, height=500)
dba.plotVenn(atac, contrast = 1, method=DBA_ALL_METHODS, bDB = T)
dev.off()



