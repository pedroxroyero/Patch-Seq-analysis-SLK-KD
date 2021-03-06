### Quality control

library(vctrs)
library(rlang)
library(scater)
library(scran)
library(edgeR)
library(ggfortify)
library(DEsingle)
library(limma)
library(biomaRt)
library(ggplot2)
library(RColorBrewer)
library(ggfortify)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(Rtsne)
library(BiocSingular)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(VennDiagram)
library(PerformanceAnalytics)
library(uwot)
library(zinbwave)
library(tidyr)
library(ggpubr)
library(cowplot)
library(DescTools)
library(plotrix)
library(EnhancedVolcano)
library(devtools)
library(disgenet2r)
library(reshape2)

# Load the data
sample_table_shSLK_vs_hrGFP <- read.csv(file = 'sample_table_shSLK_vs_hrGFP.csv', check.names = FALSE, ,row.names = 1)
count_exons_introns_shSLK_vs_hrGFP <- read.csv(file = 'count_exons_introns_shSLK_vs_hrGFP.csv', check.names = FALSE, ,row.names = 1)
genes_list_shSLK_vs_hrGFP <- read.csv(file = 'genes_list_shSLK_vs_hrGFP.csv', check.names = FALSE, ,row.names = 1)

# Create dge object
dge <- DGEList(counts=count_exons_introns_shSLK_vs_hrGFP,
               genes=genes_list_shSLK_vs_hrGFP,
               group = sample_table_shSLK_vs_hrGFP$group,
               samples=sample_table_shSLK_vs_hrGFP)

# Calculate % mitocondrial reads and total number of expressed genes for each sample
# Library size (total number of reads) of each sample is already in dge object
mito <- grep("^mt-", dge$genes$mgi_symbol)
percent_mito <- (colSums(dge$counts[mito, ])) / dge$samples$lib.size
num_genes <- colSums(dge$counts != 0)
dge$samples <- cbind(dge$samples, percent_mito=percent_mito, num_genes=num_genes)

# Create histograms of QC variables
# Samples with values deviating at least 3 times from the median absolute deviation (dashed lines) aew classified as low quality and removed from further analysis
qc_matrix <- dge$samples[dge$samples$group == 'hrGFP' | dge$samples$group == 'shSLK',]

ggplot(qc_matrix, aes(x=log(qc_matrix$lib.size), fill=as.factor(qc_matrix$group.1))) +
  geom_histogram(binwidth =0.3 ,position = 'identity', alpha=0.8, color='black') +
  labs(x="Log Libeary size", y = "N° of cells", fill="Group") + ylim(0,10) + theme_bw(base_size = 40)  +
  geom_vline(aes(xintercept=(median(log(dge$samples$lib.size)) - (3*mad(log(dge$samples$lib.size))))), color="blue", linetype="dashed", size=1) +
  scale_fill_manual(values=c("black", "red"))

ggplot(qc_matrix, aes(x=log(qc_matrix$num_genes), fill=as.factor(qc_matrix$group))) +
  geom_histogram(color='black',binwidth =0.3 ,position = 'identity', alpha=0.8) +
  labs(x="Log N° expressed genes", y = "N° of cells", fill="Group") + ylim(0,10) + theme_bw(base_size = 40)  +
  geom_vline(aes(xintercept=(median(log(dge$samples$num_genes)) - (3*mad(log(dge$samples$num_genes))))), color="blue", linetype="dashed", size=1) +
  scale_fill_manual(values=c("black", "red"))

ggplot(qc_matrix, aes(x=log(qc_matrix$percent_mito*100), fill=as.factor(qc_matrix$group))) +
  geom_histogram(binwidth =0.5 ,position = 'identity', alpha=0.8, color='black') +
  labs(x="Log % Mitocondrial genes", y = "N° of cells", fill="Group") + ylim(0,8) + theme_bw(base_size = 40)  +
  geom_vline(aes(xintercept=(median(log(qc_matrix$percent_mito*100)) + (3*mad(log(qc_matrix$percent_mito*100))))), color="blue", linetype="dashed", size=1) +
  scale_fill_manual(values=c("black", "red"))

# Check what samples do not fulfill the cutoffs 
rownames(dge$samples[log(dge$samples$lib.size) < (median(log(dge$samples$lib.size)) - (3*mad(log(dge$samples$lib.size)))),])
rownames(dge$samples[log(dge$samples$num_genes) < (median(log(dge$samples$num_genes)) - (3*mad(log(dge$samples$num_genes)))),])
rownames(dge$samples[log(dge$samples$percent_mito) < (median(log(dge$samples$percent_mito)) - (3*mad(log(dge$samples$percent_mito)))),])

# Remove bad quality samples from count matrix and sample table
libsize_drop <- isOutlier(dge$samples$lib.size, nmads = 3, type= "lower", log=TRUE) 
ngenes_drop <- isOutlier(dge$samples$num_genes, nmads = 3, type= "lower", log=TRUE) 
percentagemit_drop <- isOutlier(dge$samples$percent_mito, nmads = 3, type= "higher", log=TRUE)

count_exons_introns_shSLK_vs_hrGFP <- dge$counts[,!(libsize_drop | ngenes_drop | percentagemit_drop)] 
sample_table_shSLK_vs_hrGFP <- dge$samples[!(libsize_drop | ngenes_drop | percentagemit_drop),]

# Filtering low abundance genes (keeps genes with a read average across samples of at least 1)
ave_counts <- data.frame(rowMeans(dge$counts))
keep_genes <- ave_counts >=1
sum(keep_genes)

count_exons_introns_shSLK_vs_hrGFP <- count_exons_introns_shSLK_vs_hrGFP[keep_genes,] #remove genes that have less than 1 read as average among all samples
genes_list_shSLK_vs_hrGFP <- genes_list_shSLK_vs_hrGFP[keep_genes,]

# Remove predicted genes
predicted_genes <- genes_list_shSLK_vs_hrGFP[grep("^Gm\\d{2}", genes_list_shSLK_vs_hrGFP$mgi_symbol),]
count_exons_introns_shSLK_vs_hrGFP <- count_exons_introns_shSLK_vs_hrGFP[!(rownames(count_exons_introns_shSLK_vs_hrGFP) %in% rownames(predicted_genes)),]
genes_list_shSLK_vs_hrGFP <- genes_list_shSLK_vs_hrGFP[rownames(genes_list_shSLK_vs_hrGFP) %in% rownames(count_exons_introns_shSLK_vs_hrGFP),]

# Remove 'Striatum' and 'Unknown' samples, leaving only shSLK and hrGFP samokes
count_exons_introns_shSLK_vs_hrGFP <- count_exons_introns_shSLK_vs_hrGFP[,sample_table_shSLK_vs_hrGFP$group != "Unknown" & sample_table_shSLK_vs_hrGFP$group!= "Striatum"]
sample_table_shSLK_vs_hrGFP <- sample_table_shSLK_vs_hrGFP[sample_table_shSLK_vs_hrGFP$group!= "Unknown" & sample_table_shSLK_vs_hrGFP$group!= "Striatum",]

table(sample_table_shSLK_vs_hrGFP$group) # number of samples after QC
nrow(count_exons_introns_shSLK_vs_hrGFP) # number of genes after QC

# Save QC count matrix, sample table and gene list for further analysis
write.csv(count_exons_introns_shSLK_vs_hrGFP, file='count_exons_introns_shSLK_vs_hrGFP_QC.csv')
write.csv(genes_list_shSLK_vs_hrGFP, file='genes_list_shSLK_vs_hrGFP_QC.csv')
write.csv(sample_table_shSLK_vs_hrGFP, file='sample_table_shSLK_vs_hrGFP_QC.csv')
