library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(genomation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

blkList <- import.bed("hg38.blacklist.bed.gz")
FOXP1_peaks <- readNarrowPeak("20230510_FOXP1_CUTTag/macs2/q-0.01/FOXP1/FOXP1/3-FOXP1-1_macs2_peak_q0.01_peaks.narrowPeak")
FOXP4_peaks <- readNarrowPeak("20230510_FOXP1_CUTTag/macs2/q-0.01/FOXP1/FOXP4/3-FOXP4-1_macs2_peak_q0.01_peaks.narrowPeak")
H3K27me3_peaks <- readNarrowPeak("20230510_FOXP1_CUTTag/macs2/q-0.01/FOXP1/H3K27me3/3-H3K27me3-1_macs2_peak_q0.01_peaks.narrowPeak")
FOXP1_peaks <- FOXP1_peaks[seqnames(FOXP1_peaks) %in% unique(seqnames(blkList))]
FOXP4_peaks <- FOXP4_peaks[seqnames(FOXP4_peaks) %in% unique(seqnames(blkList))]
H3K27me3_peaks <- H3K27me3_peaks[seqnames(H3K27me3_peaks) %in% unique(seqnames(blkList))]


FOXP1_qcRes <- ChIPQCsample("20230510_FOXP1_CUTTag/alignment/bam/sorted/FOXP1/FOXP1/3-FOXP1-1_sorted.bam", chromosomes = unique(seqnames(blkList)),
                            peaks = FOXP1_peaks, annotation = "hg38", blacklist = "hg38.blacklist.bed.gz",
                            verboseT = FALSE)
saveRDS(FOXP1_qcRes, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-FOXP1-1_ChIPQCsample.rds")
FOXP4_qcRes <- ChIPQCsample("20230510_FOXP1_CUTTag/alignment/bam/sorted/FOXP1/FOXP4/3-FOXP4-1_sorted.bam", chromosomes = unique(seqnames(blkList)),
                            peaks = FOXP4_peaks, annotation = "hg38", blacklist = "hg38.blacklist.bed.gz",
                            verboseT = FALSE)
saveRDS(FOXP4_qcRes, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-FOXP4-1_ChIPQCsample.rds")
H3K27me3_qcRes <- ChIPQCsample("20230510_FOXP1_CUTTag/alignment/bam/sorted/FOXP1/H3K27me3/3-H3K27me3-1_sorted.bam", chromosomes = unique(seqnames(blkList)),
                               peaks = H3K27me3_peaks, annotation = "hg38", blacklist = "hg38.blacklist.bed.gz",
                               verboseT = FALSE)
saveRDS(H3K27me3_qcRes, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-H3K27me3-1_ChIPQCsample.rds")


QCmetrics(FOXP1_qcRes) %>% t %>% data.frame %>% dplyr:::select(Reads, starts_with(c("Filt")), 
                                                               starts_with(c("RiP")), starts_with(c("RiBL"))) %>% datatable(rownames = NULL)
flagtagcounts(FOXP1_qcRes) %>% t %>% data.frame %>% mutate(Dup_Percent = (DuplicateByChIPQC/Mapped) * 
                                                             100) %>% dplyr:::select(Mapped, Dup_Percent) %>% datatable(rownames = NULL)

#duprate 12.03%

QCmetrics(FOXP4_qcRes) %>% t %>% data.frame %>% dplyr:::select(Reads, starts_with(c("Filt")), 
                                                               starts_with(c("RiP")), starts_with(c("RiBL"))) %>% datatable(rownames = NULL)
flagtagcounts(FOXP4_qcRes) %>% t %>% data.frame %>% mutate(Dup_Percent = (DuplicateByChIPQC/Mapped) * 
                                                             100) %>% dplyr:::select(Mapped, Dup_Percent) %>% datatable(rownames = NULL)

#duprate 24.46%
QCmetrics(H3K27me3_qcRes) %>% t %>% data.frame %>% dplyr:::select(Reads, starts_with(c("Filt")), 
                                                                  starts_with(c("RiP")), starts_with(c("RiBL"))) %>% datatable(rownames = NULL)
flagtagcounts(H3K27me3_qcRes) %>% t %>% data.frame %>% mutate(Dup_Percent = (DuplicateByChIPQC/Mapped) * 
                                                                100) %>% dplyr:::select(Mapped, Dup_Percent) %>% datatable(rownames = NULL)
#duprate 21.31%

#remove blacklisted peaks#
FOXP1_peaks2 <- granges(FOXP1_qcRes)
data.frame(Blacklisted = sum(FOXP1_peaks2 %over% blkList), Not_Blacklisted = sum(!FOXP1_peaks2 %over% blkList))
FOXP1_peaks_filtered <- FOXP1_peaks2[!FOXP1_peaks2 %over% blkList]

FOXP4_peaks2 <- granges(FOXP4_qcRes)
data.frame(Blacklisted = sum(FOXP4_peaks2 %over% blkList), Not_Blacklisted = sum(!FOXP4_peaks2 %over% blkList))
FOXP4_peaks_filtered <- FOXP4_peaks2[!FOXP4_peaks2 %over% blkList]

H3K27me3_peaks2 <- granges(H3K27me3_qcRes)
data.frame(Blacklisted = sum(H3K27me3_peaks2 %over% blkList), Not_Blacklisted = sum(!H3K27me3_peaks2 %over% blkList))
H3K27me3_peaks_filtered <- H3K27me3_peaks2[!H3K27me3_peaks2 %over% blkList]

#####annotating open regions#####
library(ChIPseeker)
FOXP1_Anno <- annotatePeak(FOXP1_peaks_filtered, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
plotAnnoPie(FOXP1_Anno)

FOXP4_Anno <- annotatePeak(FOXP4_peaks_filtered, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
plotAnnoPie(FOXP4_Anno)

H3K27me3_Anno <- annotatePeak(H3K27me3_peaks_filtered, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
plotAnnoPie(H3K27me3_Anno)

FOXP1_df <- as.data.frame(FOXP1_Anno)
FOXP4_df <- as.data.frame(FOXP4_Anno)
H3K27me3_df <- as.data.frame(H3K27me3_Anno)

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes1 <- FOXP1_df$transcriptId
#FOXP1_df$geneName <- NA
G_list1 <- getBM(filters= "ensembl_transcript_id_version", attributes= c("ensembl_transcript_id_version", "ensembl_gene_id",
                                                                         "hgnc_symbol"),values=genes1,mart= mart)
write.table(FOXP1_df, "20230510_FOXP1_CUTTag/Greenleaf/het3/3-FOXP1-1_allpeaks.csv", quote = FALSE, 
            row.names = FALSE, sep = ",")
write.csv(G_list1, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-FOXP1-1_allpeaks_genes.csv")

genes4 <- FOXP4_df$transcriptId
FOXP4_df$geneName <- NA
G_list4 <- getBM(filters= "ensembl_transcript_id_version", attributes= c("ensembl_transcript_id_version", "ensembl_gene_id",
                                                                         "hgnc_symbol"),values=genes1,mart= mart)

write.table(FOXP4_df, "20230510_FOXP1_CUTTag/Greenleaf/het3/3-FOXP4-1_allpeaks.csv", quote = FALSE, 
            row.names = FALSE, sep = ",")
write.csv(G_list4, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-FOXP4-1_allpeaks_genes.csv")

genes <- H3K27me3_df$transcriptId
H3K27me3_df$geneName <- NA
G_list <- getBM(filters= "ensembl_transcript_id_version", attributes= c("ensembl_transcript_id_version", "ensembl_gene_id",
                                                                         "hgnc_symbol"),values=genes1,mart= mart)

write.table(H3K27me3_df, "20230510_FOXP1_CUTTag/Greenleaf/het3/3-H3K27me3-1_allpeaks.csv", quote = FALSE, 
            row.names = FALSE, sep = ",")
write.csv(G_list, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-H3K27me3-1_allpeaks_genes.csv")

#functional analysis#
library(rGREAT)
seqlevelsStyle(FOXP1_peaks_filtered) <- "UCSC"
FOXP1_great_Job <- submitGreatJob(FOXP1_peaks_filtered, species = "hg38")
availableCategories(FOXP1_great_Job)
FOXP1_great_ResultTable = getEnrichmentTables(FOXP1_great_Job, category = "GO")
names(FOXP1_great_ResultTable)
FOXP1_great_ResultTable[["GO Biological Process"]][1:4, ]
save(FOXP1_great_ResultTable, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-FOXP1-1-great.RData")


seqlevelsStyle(FOXP4_peaks_filtered) <- "UCSC"
FOXP4_great_Job <- submitGreatJob(FOXP4_peaks_filtered, species = "hg38")
availableCategories(FOXP4_great_Job)
FOXP4_great_ResultTable = getEnrichmentTables(FOXP4_great_Job, category = "GO")
names(FOXP4_great_ResultTable)
FOXP4_great_ResultTable[["GO Biological Process"]][1:4, ]
save(FOXP4_great_ResultTable, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-FOXP4-1-great.RData")

seqlevelsStyle(H3K27me3_peaks_filtered) <- "UCSC"
H3K27me3_great_Job <- submitGreatJob(H3K27me3_peaks_filtered, species = "hg38")
availableCategories(H3K27me3_great_Job)
H3K27me3_great_ResultTable = getEnrichmentTables(H3K27me3_great_Job, category = "GO")
names(H3K27me3_great_ResultTable)
H3K27me3_great_ResultTable[["GO Biological Process"]][1:4, ]
save(H3K27me3_great_ResultTable, file = "20230510_FOXP1_CUTTag/Greenleaf/het3/3-H3K27me3-1-great.RData")


setwd("/nowakowskilab/data1/kelsey_new/kelsey")
peaks <- dir("20230510_FOXP1_CUTTag/macs2/q-0.01/FOXP1/FOXP1/", pattern = "*.narrowPeak", 
             full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
names(myPeaks) <- c("3-FOXP1-1", "3-FOXP1-2", "7-FOXP1-1", "7-FOXP1-2")
genotype <- factor(c("3", "7"))
consensusToCount <- soGGi:::runConsensusRegions(GRangesList(myPeaks), "none")

#PCA of overlaps
library(tidyr)
myPlot <- as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(-consensusIDs) %>% 
  as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(genotype = gsub("_\\d", "", Samples)) %>% ggplot(aes(x = PC1, y = PC2, 
                                                              colour = genotype)) + geom_point(size = 5)
myPlot

library(Rsubread)
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% rowSums
table(occurrences) %>% rev %>% cumsum

bamsToCount <- dir("20230510_FOXP1_CUTTag/alignment/bam/sorted/FOXP1/FOXP1/", full.names = TRUE, pattern = "*.\\.bam$")
# indexBam(bamsToCount)
regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount), 
                                            start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
                             Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))
fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = TRUE, 
                           countMultiMappingReads = FALSE, maxFragLength = 300)
myCounts <- fcResults$counts
colnames(myCounts) <- c("3-FOXP1-1", "3-FOXP1-2", "7-FOXP1-1", "7-FOXP1-2")
save(myCounts, file = "20230510_FOXP1_CUTTag/Greenleaf/FOXP1_R514H-FOXP4.RData")

library(DESeq2)
genotype2 <- factor(c("R514H", "R514H", "WT", "WT"))
metaData <- data.frame(genotype2, row.names = colnames(myCounts))
DDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~genotype2, rowRanges = consensusToCount)
DDS <- DESeq(DDS)
Rlog <- rlog(DDS)
plotPCA(Rlog, intgroup = "genotype2", ntop = nrow(Rlog))

#difference in signal between groups#
library(DESeq2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tracktables)

cuttag <- results(DDS, c("genotype2", "R514H", "WT"), format = "GRanges")
cuttag <- cuttag[order(cuttag$pvalue)]
cuttag

#subset peaks#
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#toOverLap <- promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, 500, 500)
cuttag <- cuttag[(!is.na(cuttag$padj) & cuttag$padj < 0.05)]
makebedtable(cuttag, "FOXP1_differentialpeaks-2.html", "20230510_FOXP1_CUTTag/Greenleaf/")


