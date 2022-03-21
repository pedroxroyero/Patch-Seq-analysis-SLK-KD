### Desingle DE analysis

# Load QC count, gene and sample information
sample_table_shSLK_vs_hrGFP_QC <- read.csv(file = 'sample_table_shSLK_vs_hrGFP_QC.csv', check.names = FALSE, ,row.names = 1)
count_exons_introns_shSLK_vs_hrGFP_QC <- read.csv(file = 'count_exons_introns_shSLK_vs_hrGFP_QC.csv', check.names = FALSE, ,row.names = 1)
genes_list_shSLK_vs_hrGFP_QC <- read.csv(file = 'genes_list_shSLK_vs_hrGFP_QC.csv', check.names = FALSE, ,row.names = 1)

# Select genes expressed in at least 28 cells (size of smallest group)
table(sample_table_shSLK_vs_hrGFP_QC$group)
keep28 <- rowSums(count_exons_introns_shSLK_vs_hrGFP_QC >= 1) >= 28
nrow(count_exons_introns_shSLK_vs_hrGFP_QC[keep28,]) # Number of genes for DEsingle analysis

# Run DEsingle analysis
factor_DEsingle <- as.factor(sample_table_shSLK_vs_hrGFP_QC$group)
factor_DEsingle <- droplevels(factor_DEsingle) # remove unused factor levels

DEsingle_results <- DEsingle(counts =  count_exons_introns_shSLK_vs_hrGFP_QC[keep28,], group = factor_DEsingle)
DEsingle_results_classified <- DEtype(results = DEsingle_results, threshold = 0.05)

# Add gene information to Desingle results
DEsingle_results_classified$gene.name <- genes_list_shSLK_vs_hrGFP[rownames(DEsingle_results_classified), 3]
DEsingle_results_classified$gene.description <- genes_list_shSLK_vs_hrGFP[rownames(DEsingle_results_classified), 2]
DEsingle_results_classified$gene.type <- genes_list_shSLK_vs_hrGFP[rownames(DEsingle_results_classified), 4]

# Add column with log2 of fold change (fold change inverted to plot SLK counts with respect to hrGFP)
DEsingle_results_classified$norm_log2FoldChange <- log2(1/(DEsingle_results_classified$norm_foldChange))

# Select genes with p-value < 0.05 (to be labeled)
DEgenes <- DEsingle_results_classified[DEsingle_results_classified$pvalue.adj.FDR <= 0.05,"gene.name", drop=T]

# Create volcano plot 
EnhancedVolcano(toptable= DEsingle_results_classified, 
                lab = DEsingle_results_classified$gene.name, 
                x = 'norm_log2FoldChange', y = 'pvalue.adj.FDR', 
                xlim = c(-5.8, 5.8), ylim = c(-0.05,3.2), 
                pointSize = 3.5, labSize = 5,     labFace = 'bold', colAlpha = 0.6,
                pCutoff = 0.05, FCcutoff = 1, 
                col = c("grey30", "forestgreen", "royalblue", "red2"), 
                ylab = bquote(~-Log[10]~adjusted~italic(P)), 
                xlab =  bquote(~Log[2]~ "Norm fold change"), 
                title = "",
                subtitle = "",
                selectLab = DEgenes,
                legendLabels =  c("NS", "Norm fold change > 2", "Adjusted P < 0.05", "Norm fold change > 2 & Adjusted P < 0.05")) +
  theme(axis.text.x = element_text(size = 30, face="bold"), axis.text.y = element_text( size = 20, face="bold")) +
  theme_bw(base_size = 20) +
  guides(color=guide_legend(title=NULL), shape=guide_legend(title=NULL)) +
  theme(legend.position="top")

### Gene-disease associations

# Select DE genes from DEsingle analysis
DEsingle_results_classified <- data.frame(DEsingle_results_classified)
DEgenesnames <- DEsingle_results_classified[DEsingle_results_classified$pvalue.adj.FDR < 0.05,"gene.name", drop=T]

# Log in information to get the API key (https://www.disgenet.org/static/disgenet2r/disgenet2r.html)  
disgenet_api_key <- get_disgenet_api_key(
  email = "email@email.com", 
  password = "xxxx" )

Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

# Find gene-disease associations using DISGENET package
gene_disease <- gene2disease(gene= DEgenesnames, database = "ALL", verbose=T)
gene_disease <- disgenet2r::extract(gene_disease) # extract data from DISGENET object

#Find the unique disease classes in gene-disease matrix
disease_class_list <- vector()
for (item in c(1:length(gene_disease$disease_class_name))) {
  x = unlist(strsplit(as.character(gene_disease$disease_class_name[item]), ";"))
  disease_class_list <- as.vector(c(disease_class_list,x))
}

disease_class_list <- (levels(as.factor(trimws(disease_class_list)))) # Select unique entries

# Compute number of different diseases per disease class
perc_disease_class <- data.frame(matrix(nrow = length(disease_class_list), ncol = length(levels(as.factor(gene_disease$gene_symbol)))), row.names = disease_class_list) # Create empty data frame (columns = genes; rows = disease classes)
names(perc_disease_class) <- levels(as.factor(gene_disease$gene_symbol))
perc_disease_class[is.na(perc_disease_class)] <- 0.0

for (disease in c(1:length(disease_class_list))){
  for (item in c(1: nrow(gene_disease))) {
    if (grepl(pattern= disease_class_list[disease], x = gene_disease$disease_class_name[item])){
      perc_disease_class[as.character(disease_class_list[disease]), as.character(gene_disease$gene_symbol[item])] <- (perc_disease_class[disease, as.character(gene_disease$gene_symbol[item])] + 1)
    }      
  }
}

#Create a heatmap 
perc_disease_class <- perc_disease_class[order(-rowSums(perc_disease_class)),] #order by disease class with highest values
perc_disease_class <- perc_disease_class[1:10,] # Select the top 10 disease class
row.names(perc_disease_class)[row.names(perc_disease_class) == "Congenital, Hereditary, and Neonatal Diseases and Abnormalities"] <- c("Congenital, Hereditary, and Neonatal Diseases") #Shortening long name of diseases

perc_disease_class$disease_class <- row.names(perc_disease_class)
perc_disease_class2 <- melt(perc_disease_class) # To be used in ggplot

ggplot(perc_disease_class2[,], aes(x=variable, y= reorder(disease_class, value), fill=value)) + 
  scale_fill_gradient(low = "grey95", high = "black") +
  geom_tile() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size = 5, face="bold"), axis.text.y = element_text(size = 5, face="bold"), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_text(size = 20, face="bold")) +
  guides(fill=guide_legend(title="N° diseases")) +
  labs(x=NULL, y=NULL)
