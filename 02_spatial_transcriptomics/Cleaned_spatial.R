###############################################################################
# GeoMx DSP analysis: IBA1+ AOIs (DMH, DA) - WT vs Anx
#
# Pipeline: load -> build SPE -> sample QC -> subset -> negative probe QC ->
#           gene/ROI QC -> PCA/PRINCE -> TMM normalisation -> RUV4 batch
#           correction -> limma-voom differential expression
#
# Expected input files (in the working directory):
#   Counts_File.csv, Sample_Annotations.csv, Feature_Annotations.csv
###############################################################################


# SETUP ------------------------------------------------------------------------

library(SpatialExperiment)
library(standR)
library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggalluvial)
library(scater)
library(edgeR)
library(limma)
library(DT)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(readxl)
library(openxlsx)
library(gridExtra)
library(swamp)

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

# Small helper so every figure is saved the same way
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(filename = file.path(output_dir, filename), plot = plot,
         width = width, height = height, dpi = dpi)
}


# STEP 1: LOAD DATA --------------------------------------------------------------

countFile       <- read.csv("Counts_File.csv",         header = TRUE, check.names = FALSE)
sampleAnnoFile  <- read.csv("Sample_Annotations.csv",  header = TRUE, check.names = FALSE)
featureAnnoFile <- read.csv("Feature_Annotations.csv", header = TRUE, check.names = FALSE)

head(countFile)[, 1:5]
head(sampleAnnoFile)[, 1:5]
head(featureAnnoFile)[, 1:5]


# STEP 2: BUILD SPATIALEXPERIMENT OBJECT -------------------------------------------

spe <- readGeoMx(countFile, sampleAnnoFile, featureAnnoFile)
spe
assayNames(spe)
colData(spe)[1:5, 1:5]
rowData(spe)[1:5, 1:5]

# Round-trip through a DGEList (as in the original workflow)
dge <- edgeR::SE2DGEList(spe)
spe <- readGeoMxFromDGE(dge)
spe


# STEP 3: SAMPLE QC USING NANOSTRING FLAGS -----------------------------------------
# Remove AOIs carrying any of the flags below. AOIs that are flagged ONLY for
# "Low Negative Probe Count" are kept (negative probes are assessed in Step 5).

table(colData(spe)$QCFlags)

flags_to_remove <- c("Low Percent Aligned Reads", "Low Percent Stitched Reads",
                     "Low Raw Reads", "Low Sequencing Saturation")

qc_flags        <- colData(spe)$QCFlags
samples_to_keep <- is.na(qc_flags) | !grepl(paste(flags_to_remove, collapse = "|"), qc_flags)

cat("Removing", sum(!samples_to_keep), "out of", length(samples_to_keep),
    "samples due to QC flags\n")

spe_filtered <- spe[, samples_to_keep]

# Library sizes before / after filtering
lib_size_before  <- colSums(assay(spe, "counts"))
lib_size_after   <- colSums(assay(spe_filtered, "counts"))
lib_size_removed <- lib_size_before[!samples_to_keep]

lib_size_df <- data.frame(
  LibrarySize = c(lib_size_before, lib_size_after, lib_size_removed),
  Group = factor(
    rep(c("All Samples", "Kept Samples", "Removed Samples"),
        times = c(length(lib_size_before), length(lib_size_after), length(lib_size_removed))),
    levels = c("All Samples", "Kept Samples", "Removed Samples"))
)

QCplot <- ggplot(lib_size_df, aes(x = Group, y = LibrarySize, fill = Group)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Library Size Distribution", y = "Library Size (log10 scale)", x = "") +
  theme(legend.position = "none")

print(QCplot)
save_plot(QCplot, "library_size_qc_plot.pdf")


# STEP 4: SUBSET TO IBA1+ AOIs IN DMH, DA ------------------------------------------

neg_probes <- metadata(spe_filtered)$NegProbes
if (!is.null(neg_probes)) {
  cat("Number of negative probes:", length(neg_probes), "\n")
} else {
  cat("No negative probe information found in metadata\n")
}

colnames(colData(spe_filtered))   # available annotations
head(colData(spe_filtered))

subset_idx <- colData(spe_filtered)$iba1 == "True" &
  colData(spe_filtered)$area == "DMH, DA"
sum(subset_idx)

spe_subset_IBA1_DMH <- spe_filtered[, subset_idx]
dim(spe_subset_IBA1_DMH)

# All IBA1+ AOIs regardless of region (used only for the negative probe QC)
spe_subset_IBA1_all_regions <- spe_filtered[, colData(spe_filtered)$iba1 == "True"]


# STEP 5: NEGATIVE PROBE QC (exploratory) ------------------------------------------
# NegProbe-WTX counts are taken directly from the original count file.

wtx_row <- which(countFile$TargetName == "NegProbe-WTX")
print(paste("NegProbe-WTX found at row", wtx_row, "in count file"))

# Negative probe counts (excluding the TargetName column) and their sample names
negative_probe_counts    <- as.numeric(countFile[wtx_row, -1])
sample_names_from_counts <- colnames(countFile)[-1]


## 5a. DMH, DA subset ----------------------------------------------------------

current_sample_names <- colnames(spe_subset_IBA1_DMH)

print("Sample matching:")
print(paste("Samples in count file:", length(sample_names_from_counts)))
print(paste("Samples in current SPE object:", length(current_sample_names)))

# Match the negative probe data to the current samples
matched_negative_data <- data.frame(
  sample_id = current_sample_names,
  negative_signal = NA
)

for (i in 1:length(current_sample_names)) {
  match_idx <- which(sample_names_from_counts == current_sample_names[i])
  if (length(match_idx) > 0) {
    matched_negative_data$negative_signal[i] <- negative_probe_counts[match_idx]
  }
}

successful_matches <- sum(!is.na(matched_negative_data$negative_signal))
print(paste("Successfully matched", successful_matches, "out of", length(current_sample_names), "samples"))

# Remove samples with missing negative probe data
complete_data <- matched_negative_data[!is.na(matched_negative_data$negative_signal), ]
if (nrow(complete_data) == 0) {
  stop("No samples could be matched between count file and SPE object. Check sample naming.")
}

# Metadata and total signal for the matched samples
sample_indices <- match(complete_data$sample_id, colnames(spe_subset_IBA1_DMH))
sample_data <- colData(spe_subset_IBA1_DMH)[sample_indices, ]
total_signal_per_sample <- colSums(assay(spe_subset_IBA1_DMH, "counts")[, sample_indices, drop = FALSE])

plot_data_dmh <- data.frame(
  sample_id = complete_data$sample_id,
  negative_signal = complete_data$negative_signal,
  total_signal = total_signal_per_sample,
  genotype = sample_data$genotype,
  animal_id = sample_data$animal_id,
  aoi_nuclei = sample_data$AOINucleiCount
)

print(paste("Final dataset for plotting:", nrow(plot_data_dmh), "samples"))

# Plot 1: Negative probe signal by genotype
p1_dmh <- ggplot(plot_data_dmh, aes(x = genotype, y = negative_signal, fill = genotype)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "NegProbe-WTX Signal by Genotype",
       x = "Genotype",
       y = "Negative Probe Count",
       subtitle = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
p1_dmh

# Plot 2: Total signal vs negative probe
p2_dmh <- ggplot(plot_data_dmh, aes(x = total_signal, y = negative_signal, color = genotype)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(title = "Total Signal vs NegProbe-WTX Signal",
       x = "Total Signal (All Genes)",
       y = "Negative Probe Signal",
       subtitle = "",
       color = "Genotype") +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10()
p2_dmh

# Plot 3: Negative probe vs AOI nuclei count
p3_dmh <- ggplot(plot_data_dmh, aes(x = aoi_nuclei, y = negative_signal, color = genotype)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(title = "AOI Nuclei Count vs NegProbe-WTX Signal",
       x = "AOI Nuclei Count",
       y = "Negative Probe Signal",
       color = "Genotype") +
  theme_minimal() +
  scale_y_log10()
p3_dmh

# Plot 4: Negative probe by animal (to check for batch effects)
p4_dmh <- ggplot(plot_data_dmh, aes(x = animal_id, y = negative_signal, fill = genotype)) +
  geom_col(alpha = 0.7) +
  labs(title = "NegProbe-WTX Signal by Animal",
       x = "Animal ID",
       y = "Negative Probe Count",
       subtitle = "Check for animal-specific batch effects",
       fill = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4_dmh

# Summary statistics
print("Summary of negative probe signals by genotype:")
summary_stats <- plot_data_dmh %>%
  group_by(genotype) %>%
  summarise(
    neg_mean = mean(negative_signal),
    neg_sd = sd(negative_signal),
    neg_median = median(negative_signal),
    neg_min = min(negative_signal),
    neg_max = max(negative_signal),
    n_samples = n(),
    .groups = "drop"
  )
print(summary_stats)

# Potential filtering thresholds based on negative probe
neg_mean <- mean(plot_data_dmh$negative_signal)
neg_sd <- sd(plot_data_dmh$negative_signal)
threshold_2sd <- neg_mean + (2 * neg_sd)
threshold_3sd <- neg_mean + (3 * neg_sd)
threshold_95th <- quantile(plot_data_dmh$negative_signal, 0.95)

print("Potential filtering thresholds based on negative probe:")
print(paste("Mean negative signal:", round(neg_mean, 2)))
print(paste("Mean + 2*SD:", round(threshold_2sd, 2)))
print(paste("95th percentile:", round(threshold_95th, 2)))

# Flag potentially problematic samples (informational only - nothing is removed here)
high_negative_threshold <- threshold_95th
problematic_samples <- plot_data_dmh[plot_data_dmh$negative_signal > high_negative_threshold, ]

print(paste("Using 95th percentile threshold:", round(high_negative_threshold, 2)))
print("Samples with unusually high negative probe signals:")
if (nrow(problematic_samples) > 0) {
  print(problematic_samples[, c("sample_id", "negative_signal", "genotype", "aoi_nuclei")])
} else {
  print("No samples above threshold")
}

# STEP 6: GENE-LEVEL AND ROI-LEVEL QC (standR) --------------------------------------

spe_subset_gene_qc <- addPerROIQC(spe_subset_IBA1_DMH, rm_genes = TRUE,
                                  min_count = 10, sample_fraction = 0.9)
dim(spe_subset_gene_qc)
metadata(spe_subset_gene_qc) |> names()

# Gene QC
p_geneqc_genotype <- plotGeneQC(spe_subset_gene_qc, ordannots = "genotype", col = genotype, point_size = 2)
p_geneqc_sex      <- plotGeneQC(spe_subset_gene_qc, ordannots = "sex",      col = sex,      point_size = 2)
p_geneqc_genotype
p_geneqc_sex
save_plot(p_geneqc_genotype, "gene_qc_by_genotype.png")
save_plot(p_geneqc_sex,      "gene_qc_by_sex.png")

# ROI QC: nuclei count (threshold = 5) coloured by different annotations
p_roiqc_slide    <- plotROIQC(spe_subset_gene_qc, x_threshold = 5, color = SlideName)
p_roiqc_genotype <- plotROIQC(spe_subset_gene_qc, x_threshold = 5, color = genotype)
p_roiqc_animal   <- plotROIQC(spe_subset_gene_qc, x_threshold = 5, color = animal_id)
p_roiqc_slide
p_roiqc_genotype
p_roiqc_animal
save_plot(p_roiqc_slide,    "AOINucleiCountFrequency_SlideName.png")
save_plot(p_roiqc_genotype, "AOINucleiCountFrequency_genotype.png")
save_plot(p_roiqc_animal,   "AOINucleiCountFrequency_animalID.png")

# ROI QC: surface area vs nuclei count / library size
p_area_nuclei <- plotROIQC(spe_subset_gene_qc, x_axis = "AOISurfaceArea", x_lab = "AreaSize",
                           y_axis = "AOINucleiCount", y_lab = "AOInucleicount", col = genotype)
p_area_lib    <- plotROIQC(spe_subset_gene_qc, x_axis = "AOISurfaceArea", x_lab = "AreaSize",
                           y_axis = "lib_size", y_lab = "Library size", col = animal_id)
p_area_nuclei
p_area_lib
save_plot(p_area_nuclei, "AOISurfaceArea_vs_NucleiCount.png")
save_plot(p_area_lib,    "AOISurfaceArea_vs_LibrarySize.png")


# STEP 7: PCA AND PRINCE (TECHNICAL vs BIOLOGICAL VARIATION) -------------------------

stopifnot("animal_id" %in% colnames(colData(spe_subset_gene_qc)))

# Per-sample QC metadata (computed once, reused below)
cd <- colData(spe_subset_gene_qc)
qc_meta <- data.frame(
  animal_id            = as.factor(cd$animal_id),
  SlideName            = as.factor(cd$SlideName),
  ROILabel             = cd$ROILabel,
  area                 = cd$area,
  sex                  = as.factor(cd$sex),
  genotype             = as.factor(cd$genotype),
  AlignedReads         = cd$AlignedReads,
  DeduplicatedReads    = cd$DeduplicatedReads,
  AlignedReads_log2    = log2(cd$AlignedReads),
  AlignedToDedup_Ratio = cd$AlignedReads / cd$DeduplicatedReads,
  row.names            = colnames(spe_subset_gene_qc)
)
qc_meta

## 7a. PCA on (un-normalised) logcounts, annotated with technical factors --------
pca_result <- prcomp(t(assay(spe_subset_gene_qc, "logcounts")), scale. = TRUE)
pc_var     <- summary(pca_result)$importance[2, ] * 100

pca_plot_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], qc_meta)

pca_qc_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = genotype, shape = sex, size = AlignedToDedup_Ratio)) +
  scale_size_continuous(name = "Aligned/Dedup Ratio") +
  geom_text_repel(aes(label = animal_id), size = 3, max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "PCA with Technical and Biological Factors",
       x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pc_var[2], 1), "%)"))

print(pca_qc_plot)
save_plot(pca_qc_plot, "pca_qc_plot_IBA1_DMH.pdf", width = 10, height = 8)

# Sample most extreme on PC1 (potential outlier)
outlier_index <- which.max(abs(pca_plot_data$PC1))
cat("\nSample with most extreme PC1 value:", as.character(pca_plot_data$animal_id[outlier_index]), "\n")
cat("PC1 value:", pca_plot_data$PC1[outlier_index], "\n")
print(pca_plot_data[outlier_index, c("sex", "genotype", "AlignedReads_log2", "AlignedToDedup_Ratio")])

## 7b. swamp PRINCE plot ---------------------------------------------------------
gene_matrix <- assay(spe_subset_gene_qc, "logcounts")
sample_anno <- qc_meta[, c("animal_id", "sex", "genotype", "SlideName",
                           "AlignedReads_log2", "DeduplicatedReads", "AlignedToDedup_Ratio")]

set.seed(100)  # prince(permute = TRUE) is stochastic
prince_obj <- swamp::prince(gene_matrix, sample_anno, top = 9, permute = TRUE)

pdf(file.path(output_dir, "swamp_prince_plot_with_values.pdf"), width = 10, height = 8)
swamp::prince.plot(
  prince    = prince_obj,
  label     = colnames(prince_obj$o),
  smallest  = -20,
  note      = TRUE,       # print values in each cell
  notecol   = "black",
  notecex   = 0.7,        # small enough for numbers to fit in cells
  breaks    = seq(-20, 0, length.out = 100),
  col       = heat.colors(99),
  margins   = c(8, 10),
  key       = TRUE,
  cexRow    = 0.8,
  cexCol    = 0.8,
  xlab      = "Principal Components (GeoMX IBA1 DMH data)",
  sepcolor  = "black",
  sepwidth  = c(0.05, 0.05),
  Rsquared  = FALSE
)
dev.off()


# STEP 8: TMM NORMALISATION ---------------------------------------------------------

spe_tmm <- geomxNorm(spe_subset_gene_qc, method = "TMM")

# Relative log expression before vs after normalisation
p_before <- plotRLExpr(spe_subset_gene_qc, assay = 2, color = SlideName) +
  ggtitle("Before Normalization (Raw)") +
  theme(plot.title = element_text(hjust = 0.5))

p_after <- plotRLExpr(spe_tmm, assay = 2, color = SlideName) +
  ggtitle("After TMM Normalization") +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- grid.arrange(p_before, p_after, ncol = 2,
                              top = "Relative Log Expression (RLE) Comparison")
ggsave(file.path(output_dir, "normalization_comparison.pdf"),
       plot = combined_plot, width = 12, height = 6)


# STEP 9: BATCH CORRECTION WITH RUV4 -------------------------------------------------

# PCA before correction, to assess the need for batch correction
set.seed(100)
spe_tmm <- scater::runPCA(spe_tmm)
pca_results_tmm <- reducedDim(spe_tmm, "PCA")
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = genotype)

# Negative control genes (stable across SlideName) for RUV4
spe_bc <- findNCGs(spe_tmm, batch_name = "SlideName", top_n = 300)
metadata(spe_bc) |> names()

# k (number of RUV factors) was chosen by inspecting PCA plots
ruv_k <- 4
spe_ruv <- geomxBatchCorrection(spe_bc, factors = "genotype",
                                NCGs = metadata(spe_bc)$NCGs, k = ruv_k)
plotPairPCA(spe_ruv, assay = 2, color = genotype)


# STEP 10: DIFFERENTIAL EXPRESSION (limma-voom) --------------------------------------

# Confirm the RUV weight columns were added to colData
colData(spe_ruv)[, grep("^ruv_W", colnames(colData(spe_ruv)))] |> head()

dge_ruv <- SE2DGEList(spe_ruv)

# Design: genotype (no intercept) + one term per RUV factor
design_formula <- as.formula(
  paste("~0 + genotype +", paste0("ruv_W", seq_len(ruv_k), collapse = " + "))
)
design_ruv <- model.matrix(design_formula, data = colData(spe_ruv))
colnames(design_ruv)

# Contrast is anx - wt: positive logFC = higher in Anx.
contr.matrix <- makeContrasts(
  AnxvWT = genotypeanx - genotypewt,
  levels = colnames(design_ruv))

# Filter lowly expressed genes
keep <- filterByExpr(dge_ruv, design_ruv)
table(keep)
rownames(dge_ruv)[!keep]
dge_all <- dge_ruv[keep, ]

# voom + limma
v             <- voom(dge_all, design_ruv, plot = TRUE)
fit           <- lmFit(v)
fit_contrast  <- contrasts.fit(fit, contrasts = contr.matrix)
efit          <- eBayes(fit_contrast, robust = TRUE)

# Summary of DE genes (BH-adjusted p < 0.05)
results_efit <- decideTests(efit, p.value = 0.05)
print(summary(results_efit))

dge_results_AnxvWT_ruv      <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
dge_genes_toptable_AnxvWT_ruv <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)

# MA plot (average expression vs logFC); genes with adj.P.Val < 0.05 are labelled
ma_data <- dge_results_AnxvWT_ruv %>%
  tibble::rownames_to_column("gene") %>%
  mutate(Label = case_when(logFC > 0 & adj.P.Val < 0.05 ~ "Up Regulated",
                           logFC < 0 & adj.P.Val < 0.05 ~ "Down Regulated",
                           TRUE                         ~ "Not changed"))

ma_plot <- ggplot(ma_data, aes(AveExpr, logFC, col = Label)) +
  geom_point(shape = 1, size = 1) +
  geom_text_repel(data = dplyr::semi_join(ma_data,
                                          tibble::rownames_to_column(dge_genes_toptable_AnxvWT_ruv, "gene"),
                                          by = "gene"),
                  aes(label = gene), max.overlaps = Inf) +
  scale_color_manual(values = c("Down Regulated" = "blue",
                                "Not changed"    = "gray",
                                "Up Regulated"   = "red")) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("IBA1 DMH") +
  theme(text = element_text(size = 15))

print(ma_plot)
save_plot(ma_plot, "MA_Plot.tiff", width = 6, height = 3.5, dpi = 600)

write.csv(dge_results_AnxvWT_ruv, file.path(output_dir, "DEGS.csv"), row.names = TRUE)


# SESSION INFO (for reproducibility) -------------------------------------------------
sessionInfo()
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
