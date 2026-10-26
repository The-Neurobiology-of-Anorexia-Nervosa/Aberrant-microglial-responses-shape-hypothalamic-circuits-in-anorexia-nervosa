library(ggplot2)
library(biomaRt)
library(dplyr)
library(DESeq2)
library(ggrepel)
library(tibble)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

path <- "set/your/path"
setwd(path)


# 1. Data preparation
# 1.1 Import data in correct format
counts_data <- read.csv("gene_counts_all_samples_annot.csv", skip = 1, header = FALSE)
sample_names <- as.character(counts_data[1, -1])
count_matrix <- counts_data[-1, -1]
gene_names <- as.character(counts_data[-1, 1])
count_matrix <- as.matrix(sapply(count_matrix, as.numeric))
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- sample_names
count_matrix <- count_matrix[, colnames(count_matrix) != "blank"]


# 1.2 Make a sample information metadatasheet
#Extract the main condition (AN or CT) and treatment (GLP1, untreated, HTsyn etc.)
main_condition <- sub("([A-Z]+).*", "\\1", colnames(count_matrix))
treatment_condition <- sapply(colnames(count_matrix), function(x) {
  # If no suffix after the last underscore, it's "untreated"
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) < 3) {
    return("untreated")
  }
  return(parts[length(parts)])
})

# Create sample information dataframe 
sample_info <- data.frame(
  sample = colnames(count_matrix),
  main_group = main_condition,
  treatment = treatment_condition,
  group = paste(main_condition, treatment_condition, sep = "_"),
  row.names = colnames(count_matrix)
)

write.csv(sample_info, "sample_metadata.csv", row.names = TRUE)

# 1.3 Create a dataframe and annotate with ensembl gene names and biotypes
ensembl <- useEnsembl(
  biomart = "genes", 
  dataset = "hsapiens_gene_ensembl", 
  mirror = "useast" 
)

gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = gene_names,
  mart = ensembl
) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE)


count_matrix_df <- as.data.frame(count_matrix)
count_matrix_df$ensembl_gene_id <- rownames(count_matrix_df)
count_matrix_annot <- left_join(count_matrix_df, gene_annotations, by = "ensembl_gene_id")
cat("Number of genes:", nrow(count_matrix_annot), "\n")

# 1.4 Sample Filtering 
# 1.4.1 Filter to only include protein coding genes
count_matrix_prot <- count_matrix_annot %>% filter(gene_biotype == "protein_coding")
cat("Number of protein coding genes:", nrow(count_matrix_prot), "\n")

# 1.4.2 Filter treatment condition (in this example only samples which have been GLP-1 stimulated)
treated_samples <- sample_info %>% filter(treatment == "glp")

# 1.4.3 Extract metadata indexed by ensembl_gene_id
gene_metadata <- count_matrix_prot %>%
  select(ensembl_gene_id, external_gene_name, gene_biotype) %>%
  column_to_rownames("ensembl_gene_id")

treated_counts <- count_matrix_prot %>%
  select(all_of(treated_samples$sample)) %>%
  as.matrix()

rownames(treated_counts) <- count_matrix_prot$ensembl_gene_id

treated_sample_metadata <- sample_info[treated_samples$sample, ]
treated_sample_metadata$main_group <- factor(
  treated_sample_metadata$main_group, 
  levels = c("CT", "AN")
)

# 2 Quality control and DESeq2 object creation
# 2.1 DESeq2 object creation
dds <- DESeqDataSetFromMatrix(
  countData = treated_counts,
  colData = treated_sample_metadata,
  design = ~ main_group,
  rowData = gene_metadata)

# 2.2 Filtering low gene counts (lowly expressed genes)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]
cat("Number of genes after filtering:", nrow(dds), "\n")

# 2.3 Normalization and checking size factors 
dds_norm <- DESeq(dds)
sizeFactors(dds_norm)

# 2.4 Diagnostic plots
# Plot size factors
size_factor_df <- data.frame(
  sample = names(sizeFactors(dds_norm)),
  size_factor = sizeFactors(dds_norm)
)

# 2.4.1 size factor plot
ggplot(size_factor_df, aes(x = sample, y = size_factor)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sample Size Factors", 
       y = "Size Factor", 
       x = "Sample")

# 2.4.2 Variance stabilizing transformation
rld_1 <- rlog(dds_norm, blind = FALSE)
# Plot dispersion estimates
plotDispEsts(dds_norm)

# 2.4.3 PCA plot with sample names
p <- plotPCA(rld_1, intgroup = "main_group")
p + geom_text_repel(aes(label = rownames(p$data)), 
                    max.overlaps = Inf) +
  theme_minimal() +
  ggtitle("PCA of Baseline Samples")

library_sizes <- colSums(counts(dds_norm))
log_library_sizes <- log2(library_sizes)

# 2.4.4 Violin plot of library sizes
#Create a data frame for plotting
library_size_df <- data.frame(
  sample = names(library_sizes),
  group = dds_norm$main_group,
  log2_library_size = log_library_sizes
)

# Violin plot of library sizes
ggplot(library_size_df, aes(x = group, y = log2_library_size, fill = group)) +
  geom_violin() +
  geom_jitter(width = 0.2) +
  ggtitle("Distribution of Library Sizes by Group") +
  ylab("Log2 Library Size") +
  theme_minimal()

par(mfrow=c(1,2))

# 2.4.5 Box plots of counts before and after normalization
# Raw Counts Boxplot
boxplot(log2(counts(dds) + 1), 
        main="Raw Counts", 
        las=2,  # rotate labels
        cex.axis=0.7,
        ylab="log2(count + 1)")  # reduce text size

# Normalized Counts Boxplot
boxplot(log2(counts(dds_norm, normalized = TRUE) + 1), 
        main="Normalized Counts", 
        las=2,  # rotate labels
        cex.axis=0.7)  # reduce text size

# After QC use the diagnostic plots to determine outliers for removal (in this case the conclusion is to remove CT17)
# Repeat steps 2.1-2.3 making a new DESeq2 object without the outlier
# Define vector of samples to remove
low_quality_samples <- c("CT17_8_glp", "CT17_9_glp")

# Remove from count matrix
treated_counts2 <- treated_counts[, !colnames(treated_counts) %in% low_quality_samples]

# Remove from sample metadata
treated_sample_metadata2 <- treated_sample_metadata[
  !treated_sample_metadata$sample %in% low_quality_samples, 
]

# New DESeq2 object
dds2 <- DESeqDataSetFromMatrix(
  countData = treated_counts2,
  colData = treated_sample_metadata2,
  design = ~ main_group,
  rowData = gene_metadata)


# filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds2) >= 5) >= smallestGroupSize
dds2 <- dds2[keep,]

# Normalizatrion
dds2_norm <- DESeq(dds2)
cat("Final sample count:", ncol(dds2_norm), "\n")


# 3.1 Differential Gene Expression Analysis
#Extract results with gene names
res_df_with_genes <- as.data.frame(results(dds2_norm, contrast = c("main_group", "AN", "CT")))
res_df_with_genes$ensembl_gene_id <- rownames(res_df_with_genes)
res_df_with_genes$external_gene_name <- mcols(dds2_norm)$external_gene_name
res_df_with_genes$gene_name <- ifelse(
  is.na(res_df_with_genes$external_gene_name) | res_df_with_genes$external_gene_name == "",
  res_df_with_genes$ensembl_gene_id,
  res_df_with_genes$external_gene_name
)
res_df_with_genes <- res_df_with_genes[, c(
  "ensembl_gene_id", "external_gene_name", "gene_name",
  "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
)]

# 3.2 Check number of NA results after DESeq2
cat("Total genes with NA padj:", sum(is.na(res_df_with_genes$padj)), "\n")

# 3.3 Export DEGs
write.csv(res_df_with_genes, "Results_GLP1.csv", row.names = FALSE)

# 3.4 Make an MA plot with DEGs
MAplot <- res_df_with_genes[!is.na(res_df_with_genes$padj), ] %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up Regulated",
      padj < 0.05 & log2FoldChange < 0 ~ "Down Regulated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  ggplot(aes(x = log2(baseMean), y = log2FoldChange, color = Significance)) +
  geom_point(shape = 1, size = 2) +
  geom_text_repel(
    data = . %>% filter(padj < 0.05),
    aes(label = gene_name),  
    max.overlaps = 20,
    size = 3
  ) +
  theme_bw() +
  xlab("Log2 BaseMean") +
  ylab("Log2 Fold-Change") +
  ggtitle("GLP-1 stimulated AN vs HC microglia") +
  scale_color_manual(values = c(
    "Up Regulated" = "red", 
    "Down Regulated" = "blue", 
    "Not Significant" = "grey"
  )) +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  )
print(MAplot)

ggsave("MAplot_GLP1.tiff", 
       volcano_plot, 
       width = 9, 
       height = 7, 
       dpi = 600) 


# 4. GSEA analysis
# 4.1 Read data DEGs
df <- read.csv("Results_GLP1.csv")

# 4.2 Set human (organism) 
organism = "org.Hs.eg.db"

# 4.3 Calculate the combined metric for ranking
df$combined_metric <- -log10(df$pvalue) * sign(df$log2FoldChange)

#Create a named vector with the combined metric and gene names
gene_list <- df$combined_metric
names(gene_list) <- df$gene_name

# 4.4 Omit any NA values
gene_list <- na.omit(gene_list)

# 4.5 Sort the list in decreasing order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)
keytypes(org.Hs.eg.db)

# 4.6 Handle duplicated genes
duplicated_genes <- names(gene_list)[duplicated(names(gene_list))]
print(duplicated_genes)

set.seed(42)  

remove_duplicates <- function(gene_list) {
  duplicate_genes <- names(gene_list)[duplicated(names(gene_list))]
  
  
  keep_genes <- rep(TRUE, length(gene_list))
  
  
  for (gene in unique(duplicate_genes)) {
    
    gene_indices <- which(names(gene_list) == gene)
    
    
    keep_index <- sample(gene_indices, 1)
    
    
    keep_genes[gene_indices[gene_indices != keep_index]] <- FALSE
  }
  
  
  return(gene_list[keep_genes])
}

# Apply the function to remove duplicates
gene_list <- remove_duplicates(gene_list)

# 4.7 Perform GSE analysis
gse <- gseGO(geneList = gene_list, 
              ont = "BP", 
              keyType = "SYMBOL", 
              minGSSize = 5, 
              maxGSSize = 800,
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = organism, 
              pAdjustMethod = "BH")

#Export GSEA results
write.csv(gse, "GSEA_GLP1.csv")


# 4.8 Make a dataframe and format direction of pathways for plotting
gse_df <- as.data.frame(gse)
gse_df$direction <- ifelse(gse_df$NES > 0, "Upregulated", "Downregulated")

# Create top 20 up and downregulated pathways object for plotting
gse_top20 <- bind_rows(
  gse_df %>% filter(direction == "Upregulated") %>% slice_max(NES, n = 10),
  gse_df %>% filter(direction == "Downregulated") %>% slice_min(NES, n = 10)
)

# 4.9 Export GSEA plot
GSEA_plot <- ggplot(gse_top20, 
                        aes(x = NES, 
                            y = reorder(Description, NES),
                            fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(
    title = "GSEA: GO Biological Processes \n GLP-1 stimulated AN vs HC Microglia",
    x = "Normalized Enrichment Score (NES)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

print(GSEA_plot)

#Save plot
ggsave("GSEA_GOBP_GLP1.tiff", plot = GSEA_plot,
       width = 8, 
       height = 5,
       dpi = 600)

