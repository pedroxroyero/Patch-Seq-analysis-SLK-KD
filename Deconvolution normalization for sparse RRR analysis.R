### Deconvolution normalization for sparse RRR analysis

# Load the QC data
sample_table_shSLK_vs_hrGFP_QC <- read.csv(file = 'sample_table_shSLK_vs_hrGFP_QC.csv', check.names = FALSE, ,row.names = 1)
count_exons_introns_shSLK_vs_hrGFP_QC <- read.csv(file = 'count_exons_introns_shSLK_vs_hrGFP_QC.csv', check.names = FALSE, ,row.names = 1)
genes_list_shSLK_vs_hrGFP_QC <- read.csv(file = 'genes_list_shSLK_vs_hrGFP_QC.csv', check.names = FALSE, ,row.names = 1)

# Create a single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = count_exons_introns_shSLK_vs_hrGFP_QC),  colData = sample_table_shSLK_vs_hrGFP_QC)

# Creating cluster to be used in the deconvulution method
clust.sce <- quickCluster(sce, min.size = 5) 

# Compute size factors 
sce <- computeSumFactors(sce, cluster = clust.sce)
summary(sizeFactors(sce))

# Visualize size factors vs library size per sample
sizeFactors <- data.frame(size_factors = sizeFactors(sce), libsize = sce$lib.size, group = sce$group, batch = sce$Batch, row.names = colnames(sce))

ggplot(sizeFactors, aes(x=size_factors, y=libsize/1e6, color=group, shape=sce$Batch)) +
  geom_point(size=4) + labs(x="Size factor", y = "Library size (millions)") + theme_bw(base_size = 30) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(trans='log') +
  scale_y_continuous(trans='log') +
  theme(axis.text.x = element_text(size = 15, face="bold"), axis.text.y = element_text(size = 15, face="bold")) +
  guides(color=guide_legend(title=NULL))

# Apply factors to matrix 
sce <- logNormCounts(sce,log=FALSE)
count_exons_introns_shSLK_vs_hrGFP_QC_Norm <- normcounts(sce)

# Save normalized count matrix
write.csv(count_exons_introns_shSLK_vs_hrGFP_QC_Norm, file='count_exons_introns_shSLK_vs_hrGFP_QC_Norm.csv')

