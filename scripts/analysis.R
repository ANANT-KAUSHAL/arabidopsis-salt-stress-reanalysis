
# ============================================================
# PROJECT: Arabidopsis Salt Stress Transcriptomics
# AUTHOR: Anant Kaushal
# DATASET: GSE16765
# PLATFORM: Affymetrix ATH1 (GPL198)
# ANALYSIS: Exploratory differential expression using limma
# ============================================================

# 1. LOAD LIBRARIES
library(GEOquery)
library(limma)

# 2. IDENTIFY PROJECT DIRECTORY AND CREATE OUTPUT DIRECTORIES

# Identify the location of this script when running in RStudio
script_file <- NULL

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  script_file <- rstudioapi::getActiveDocumentContext()$path
}

# Identify the script location when running through Rscript
if (is.null(script_file) || !nzchar(script_file)) {
  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  
  if (length(file_argument) > 0) {
    script_file <- sub("^--file=", "", file_argument[1])
  }
}

# Stop if the script location cannot be identified
if (is.null(script_file) || !nzchar(script_file)) {
  stop("Unable to identify the project directory. Open and run this script from RStudio or use Rscript.")
}

# The project directory is the parent of the scripts folder
project <- normalizePath(
  file.path(dirname(script_file), ".."),
  winslash = "/",
  mustWork = TRUE
)
# Set the working directory to the project root
setwd(project)

# Define project output directories
figure_dir <- file.path(project, "figures")
result_dir <- file.path(project, "results")

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# 3. LOAD DATA
message("Loading GSE16765...")

gse_salt <- getGEO("GSE16765", GSEMatrix = TRUE)
salt_data <- gse_salt[[1]]

# 4. EXTRACT EXPRESSION AND ANNOTATION
expression_matrix <- exprs(salt_data)
annotation <- fData(salt_data)

stopifnot(
  nrow(expression_matrix) == 22810,
  ncol(expression_matrix) == 12,
  "Gene Symbol" %in% colnames(annotation),
  "Gene Title" %in% colnames(annotation)
)

# 5. VERIFY AND SELECT COL SAMPLES
selected_ids <- c(
  "GSM420236", "GSM420237",
  "GSM420238", "GSM420239"
)

sample_metadata <- pData(salt_data)

sample_index <- match(
  selected_ids,
  sample_metadata$geo_accession
)

if (anyNA(sample_index)) {
  stop("One or more selected GSM identifiers were not found.")
}

selected_metadata <- sample_metadata[sample_index, ]

stopifnot(
  all(selected_metadata$`genotype:ch1` == "Col")
)

expected_titles <- c(
  "Col_leaf_0mM NaCl_rep1",
  "Col_leaf_0mM NaCl_rep2",
  "Col_leaf_100mM NaCl_rep1",
  "Col_leaf_100mM NaCl_rep2"
)

if (!identical(
  as.character(selected_metadata$title),
  expected_titles
)) {
  stop("Sample metadata does not match expected treatment order.")
}

col_samples <- expression_matrix[
  ,
  selected_ids,
  drop = FALSE
]

# 6. EXPERIMENTAL DESIGN
groups <- factor(
  c("Control", "Control", "Salt", "Salt"),
  levels = c("Control", "Salt")
)

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# 7. DIFFERENTIAL EXPRESSION USING LIMMA
message("Running differential-expression analysis...")

fit <- lmFit(col_samples, design)

contrast_matrix <- makeContrasts(
  Diff = Salt - Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# 8. INTEGRATE PROBE ANNOTATION
annotation_index <- match(
  rownames(all_results),
  rownames(annotation)
)

if (anyNA(annotation_index)) {
  stop("Some expression probes could not be matched to annotation.")
}

result_table <- data.frame(
  Probe_ID = rownames(all_results),
  annotation[
    annotation_index,
    c(
      "Gene Symbol",
      "Gene Title",
      "AGI",
      "ENTREZ_GENE_ID"
    ),
    drop = FALSE
  ],
  all_results,
  check.names = FALSE,
  row.names = NULL
)

# 9. EXPORT COMPLETE RESULTS
write.csv(
  result_table,
  "results/Salt_Stress_Full_Results.csv",
  row.names = FALSE
)

# 10. GENERATE VOLCANO PLOT
message("Generating volcano plot...")

plot_data <- all_results

plot_data$negLog10AdjP <- -log10(
  pmax(
    plot_data$adj.P.Val,
    .Machine$double.xmin
  )
)

plot_data$Color <- ifelse(
  plot_data$adj.P.Val < 0.05 &
    abs(plot_data$logFC) > 1,
  "red",
  "grey"
)

png(
  file.path(figure_dir, "volcano_plot_salt.png"),
  width = 1200,
  height = 900,
  res = 150
)

plot(
  plot_data$logFC,
  plot_data$negLog10AdjP,
  col = plot_data$Color,
  pch = 16,
  cex = 0.65,
  main = "Differential Expression: Salt vs Control",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value"
)

abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)

tch4_id <- "247925_at"

if (tch4_id %in% rownames(plot_data)) {
  
  tch4 <- plot_data[tch4_id, ]
  
  points(
    tch4$logFC,
    tch4$negLog10AdjP,
    col = "darkgreen",
    pch = 19,
    cex = 1.3
  )
  
  text(
    tch4$logFC,
    tch4$negLog10AdjP,
    labels = "TCH4",
    pos = 2,
    font = 2
  )
}

dev.off()

# 11. IDENTIFY WRKY-ASSOCIATED PROBE SETS
message("Identifying WRKY-associated annotations...")

gene_symbols <- ifelse(
  is.na(annotation[["Gene Symbol"]]),
  "",
  as.character(annotation[["Gene Symbol"]])
)

gene_titles <- ifelse(
  is.na(annotation[["Gene Title"]]),
  "",
  as.character(annotation[["Gene Title"]])
)

wrky_annotation_match <- grepl(
  "WRKY",
  paste(gene_symbols, gene_titles),
  ignore.case = TRUE
)

wrky_probe_ids <- rownames(annotation)[
  wrky_annotation_match
]

wrky_results <- result_table[
  result_table$Probe_ID %in% wrky_probe_ids,
  ,
  drop = FALSE
]

write.csv(
  wrky_results,
  "results/WRKY_Associated_Probe_Results.csv",
  row.names = FALSE
)

message(
  "WRKY-associated annotated probe sets identified: ",
  length(wrky_probe_ids)
)

# 12. GENERATE WRKY EXPRESSION HEATMAP
if (length(wrky_probe_ids) > 0) {
  
  wrky_expression <- col_samples[
    intersect(
      wrky_probe_ids,
      rownames(col_samples)
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(wrky_expression) > 0) {
    
    wrky_scaled <- t(
      scale(t(wrky_expression))
    )
    
    valid_rows <- apply(
      wrky_scaled,
      1,
      function(x) all(is.finite(x))
    )
    
    wrky_scaled <- wrky_scaled[
      valid_rows,
      ,
      drop = FALSE
    ]
    
    if (nrow(wrky_scaled) > 0) {
      
      png(
        file.path(figure_dir, "wrky_heatmap.png"),
        width = 1400,
        height = 1000,
        res = 150
      )
      
      heatmap(
        wrky_scaled,
        scale = "none",
        Colv = NA,
        margins = c(8, 8),
        main = "WRKY-associated Expression Patterns",
        xlab = "Experimental Samples",
        ylab = "Probe Sets"
      )
      
      dev.off()
    }
  }
}

# 13. ANALYSIS SUMMARY
summary_table <- data.frame(
  Metric = c(
    "Total probe sets",
    "Total samples",
    "Selected samples",
    "Control replicates",
    "Salt replicates",
    "WRKY-associated annotated probe sets",
    "Adjusted P-value < 0.05",
    "Adjusted P-value < 0.05 and |logFC| > 1"
  ),
  Value = c(
    nrow(expression_matrix),
    ncol(expression_matrix),
    ncol(col_samples),
    2,
    2,
    length(wrky_probe_ids),
    sum(all_results$adj.P.Val < 0.05, na.rm = TRUE),
    sum(
      all_results$adj.P.Val < 0.05 &
        abs(all_results$logFC) > 1,
      na.rm = TRUE
    )
  )
)

write.csv(
  summary_table,
  "results/Analysis_Summary.csv",
  row.names = FALSE
)

# 14. COMPLETION MESSAGE
message("------------------------------------------")
message("ANALYSIS COMPLETED")
message("Check the figures and results folders.")
message("------------------------------------------")
wrky_results[
  wrky_results$adj.P.Val < 0.05 &
    abs(wrky_results$logFC) > 1,
  c(
    "Probe_ID",
    "Gene Symbol",
    "Gene Title",
    "logFC",
    "P.Value",
    "adj.P.Val"
  )
]
  wrky_expression_check <- col_samples[
    wrky_results$Probe_ID,
    ,
    drop = FALSE
  ]
  
  print(round(wrky_expression_check, 3))
  boxplot(
    col_samples,
    main = "Expression Distribution Across Samples",
    ylab = "Expression Intensity",
    las = 2,
    col = c("lightblue", "lightblue",
            "salmon", "salmon")
  )
  unique(pData(salt_data)$data_processing)
  summary(as.vector(col_samples))
  sample_cor <- cor(col_samples, method = "pearson")
  
  round(sample_cor, 3)
  # Perform PCA on the four samples
  pca_result <- prcomp(t(col_samples), scale. = FALSE)
  
  # Calculate percentage of variance explained
  pca_variance <- (pca_result$sdev^2 /
                     sum(pca_result$sdev^2)) * 100
  
  # Plot the first two principal components
  plot(
    pca_result$x[, 1],
    pca_result$x[, 2],
    col = c("blue", "blue", "red", "red"),
    pch = 19,
    cex = 1.5,
    xlab = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
    ylab = paste0("PC2 (", round(pca_variance[2], 2), "%)"),
    main = "PCA of Arabidopsis Salt-Stress Samples"
  )
  
  text(
    pca_result$x[, 1],
    pca_result$x[, 2],
    labels = colnames(col_samples),
    pos = 3,
    cex = 0.8
  )
  
  legend(
    "topright",
    legend = c("Control", "Salt"),
    col = c("blue", "red"),
    pch = 19
  )
  round(pca_result$x, 3)
  pData(salt_data)[
    c("GSM420236", "GSM420237", "GSM420238", "GSM420239"),
    c("title", "characteristics_ch1",
      "characteristics_ch1.1", "characteristics_ch1.2",
      "treatment_protocol_ch1")
  ]
  # Identify the 20 probe sets contributing most strongly to PC2
pc2_loadings <- pca_result$rotation[, 2]

top_pc2_probes <- names(
  sort(abs(pc2_loadings), decreasing = TRUE)
)[1:20]

# Display their probe IDs and loading values
data.frame(
  Probe_ID = top_pc2_probes,
  PC2_Loading = pc2_loadings[top_pc2_probes]
)
# Retrieve annotation for the 20 highest-loading probe sets
top_pc2_annotation <- annotation[
  match(top_pc2_probes, annotation$ID),
  c("ID", "Gene Symbol", "Gene Title", "AGI")
]

# Retrieve their expression values in all four samples
top_pc2_expression <- col_samples[
  top_pc2_probes, ,
  drop = FALSE
]

# Combine annotation and expression values
top_pc2_details <- cbind(
  top_pc2_annotation,
  round(top_pc2_expression, 3)
)

print(top_pc2_details)
# Percentage of PC2 loading magnitude represented by the top 20 probes
pc2_top20_contribution <- sum(
  pc2_loadings[top_pc2_probes]^2
) * 100

round(pc2_top20_contribution, 2)
# Calculate the mean expression difference between salt replicates
salt_replicate_difference <- data.frame(
  Probe_ID = top_pc2_probes,
  Gene = top_pc2_annotation$`Gene Symbol`,
  Salt_1 = col_samples[top_pc2_probes, "GSM420238"],
  Salt_2 = col_samples[top_pc2_probes, "GSM420239"]
)

salt_replicate_difference$Difference <-
  salt_replicate_difference$Salt_2 -
  salt_replicate_difference$Salt_1

print(salt_replicate_difference)
boxplot(
  top_pc2_expression[, c("GSM420238", "GSM420239")],
  main = "PC2-associated Probe Expression: Salt Replicates",
  ylab = "RMA Expression Intensity",
  names = c("Salt Replicate 1", "Salt Replicate 2"),
  col = c("salmon", "lightblue")
)
# Calculate expression differences across all probe sets
all_probe_differences <-
  col_samples[, "GSM420239"] -
  col_samples[, "GSM420238"]

# Summarise the differences
summary(all_probe_differences)

# Calculate the proportion of probes with higher expression
# in Salt Replicate 2
mean(all_probe_differences > 0) * 100
# Identify the 20 largest absolute expression differences
largest_differences <- sort(
  abs(all_probe_differences),
  decreasing = TRUE
)[1:20]

# Retrieve their probe IDs
top_difference_probes <- names(largest_differences)

# Retrieve their annotations
top_difference_annotation <- annotation[
  match(top_difference_probes, annotation$ID),
  c("ID", "Gene Symbol", "Gene Title", "AGI")
]

# Display the results
data.frame(
  top_difference_annotation,
  Difference = all_probe_differences[top_difference_probes]
)
# Count probe sets with absolute expression differences > 2
sum(abs(all_probe_differences) > 2)

# Calculate their percentage
mean(abs(all_probe_differences) > 2) * 100
# Identify all probe sets with large differences
large_diff_ids <- names(
  all_probe_differences[
    abs(all_probe_differences) > 2
  ]
)

# Retrieve their annotations
large_diff_annotation <- annotation[
  match(large_diff_ids, annotation$ID),
  c("ID", "Gene Symbol", "Gene Title", "AGI")
]

# Count how many have gene symbols
sum(!is.na(large_diff_annotation$`Gene Symbol`) &
      large_diff_annotation$`Gene Symbol` != "")

# Display the 30 largest differences
large_diff_order <- order(
  abs(all_probe_differences[large_diff_ids]),
  decreasing = TRUE
)

head(
  data.frame(
    large_diff_annotation[large_diff_order, ],
    Difference = all_probe_differences[
      large_diff_ids[large_diff_order]
    ]
  ),
  30
)
sum(
  !is.na(large_diff_annotation$`Gene Symbol`) &
    large_diff_annotation$`Gene Symbol` != ""
)
colnames(large_diff_annotation)
sum(
  !is.na(large_diff_annotation$`Gene Symbol`) &
    large_diff_annotation$`Gene Symbol` != ""
)
unique(large_diff_annotation$`Gene Symbol`)
# Combine gene symbols and gene descriptions
large_diff_text <- paste(
  large_diff_annotation$`Gene Symbol`,
  large_diff_annotation$`Gene Title`
)

# Identify chloroplast/photosynthesis-associated annotations
chloroplast_pattern <- paste(
  "photosystem|chloroplast|thylakoid|",
  "photosynthetic|NADH dehydrogenase|",
  "ATP synthase|cytochrome b6",
  sep = ""
)

chloroplast_candidates <- grepl(
  chloroplast_pattern,
  large_diff_text,
  ignore.case = TRUE
)

# Count matching probe sets
sum(chloroplast_candidates)

# Display the matching annotations
large_diff_annotation[chloroplast_candidates, ]
pData(salt_data)[,
                 c("title", "genotype:ch1", "tissue:ch1",
                   "age:ch1")
]
# Prepare metadata for all 12 samples
sample_metadata <- pData(salt_data)

# Extract genotype and treatment information
sample_metadata$Genotype <- factor(
  sample_metadata$`genotype:ch1`,
  levels = c("Ws", "Col", "Col(gl)")
)

sample_metadata$Treatment <- factor(
  ifelse(
    grepl("100mM", sample_metadata$title),
    "Salt",
    "Control"
  ),
  levels = c("Control", "Salt")
)

# Display the experimental design
design_table <- table(
  sample_metadata$Genotype,
  sample_metadata$Treatment
)

print(design_table)
# Create genotype-treatment group
sample_metadata$Group <- interaction(
  sample_metadata$Genotype,
  sample_metadata$Treatment,
  sep = "_",
  drop = TRUE
)

# Arrange expression matrix in the same sample order
all_samples_expression <- exprs(salt_data)[
  ,
  rownames(sample_metadata),
  drop = FALSE
]

# Create the design matrix for the six genotype-treatment groups
design <- model.matrix(
  ~ 0 + Group,
  data = sample_metadata
)

# Assign consistent names to the six experimental groups
colnames(design) <- c(
  "Ws_Control",
  "Col_Control",
  "Colgl_Control",
  "Ws_Salt",
  "Col_Salt",
  "Colgl_Salt"
)

# Assign sample IDs to the design matrix rows
rownames(design) <- colnames(all_samples_expression)

# Verify sample alignment and model structure
stopifnot(
  identical(
    colnames(all_samples_expression),
    rownames(design)
  ),
  qr(design)$rank == ncol(design),
  nrow(design) == 12,
  ncol(design) == 6
)

# Display matrix dimensions and rank
cat("Expression matrix dimensions:\n")
print(dim(all_samples_expression))

cat("\nDesign matrix dimensions:\n")
print(dim(design))

cat("\nDesign matrix rank:\n")
print(qr(design)$rank)

# Fit the linear model across all 22,810 probe sets
fit <- limma::lmFit(
  all_samples_expression,
  design
)

# Define genotype-specific salt responses and interactions
contrast_matrix <- limma::makeContrasts(
  Ws_Salt_Response = Ws_Salt - Ws_Control,
  Col_Salt_Response = Col_Salt - Col_Control,
  Colgl_Salt_Response = Colgl_Salt - Colgl_Control,
  
  Col_vs_Ws_Interaction =
    (Col_Salt - Col_Control) -
    (Ws_Salt - Ws_Control),
  
  Colgl_vs_Ws_Interaction =
    (Colgl_Salt - Colgl_Control) -
    (Ws_Salt - Ws_Control),
  
  levels = design
)

# Apply contrasts and empirical Bayes moderation
fit_contrasts <- limma::contrasts.fit(
  fit,
  contrast_matrix
)

fit_contrasts <- limma::eBayes(fit_contrasts)

# Confirm the five comparisons
colnames(contrast_matrix)
# Count significant probe sets in each comparison
comparison_names <- colnames(contrast_matrix)

de_summary <- data.frame(
  Comparison = comparison_names,
  Significant_FDR_0.05 = sapply(comparison_names, function(x) {
    results <- limma::topTable(
      fit_contrasts,
      coef = x,
      number = Inf,
      adjust.method = "BH"
    )
    sum(results$adj.P.Val < 0.05, na.rm = TRUE)
  })
)

print(de_summary)
# Count statistically significant and substantial expression changes
de_summary_effect <- data.frame(
  Comparison = comparison_names,
  FDR_Significant = integer(length(comparison_names)),
  FDR_and_FoldChange = integer(length(comparison_names))
)

for (i in seq_along(comparison_names)) {

  results <- limma::topTable(
    fit_contrasts,
    coef = comparison_names[i],
    number = Inf,
    adjust.method = "BH"
  )

  de_summary_effect$FDR_Significant[i] <-
    sum(results$adj.P.Val < 0.05, na.rm = TRUE)

  de_summary_effect$FDR_and_FoldChange[i] <-
    sum(
      results$adj.P.Val < 0.05 &
      abs(results$logFC) > 1,
      na.rm = TRUE
    )
}

print(de_summary_effect)
# Create output directory if it does not exist
# Define the project results directory
result_dir <- file.path(project, "results")

# Create the results directory if it does not exist
dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# Obtain platform annotation
annotation <- as.data.frame(fData(salt_data))

# Generate and save complete results for all five comparisons
for (comparison in comparison_names) {
  
  results <- limma::topTable(
    fit_contrasts,
    coef = comparison,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  results$Probe_ID <- rownames(results)
  
  annotation_index <- match(
    results$Probe_ID,
    annotation$ID
  )
  
  results$Gene_Symbol <- annotation$`Gene Symbol`[annotation_index]
  results$Gene_Title <- annotation$`Gene Title`[annotation_index]
  results$AGI <- annotation$AGI[annotation_index]
  
  results <- results[, c(
    "Probe_ID",
    "Gene_Symbol",
    "Gene_Title",
    "AGI",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )]
  
  write.csv(
    results,
    file = file.path(
      result_dir,
      paste0(comparison, "_Full_Results.csv")
    ),
    row.names = FALSE
  )
}

cat("Saved five complete differential-expression result tables.\n")
# Create figures directory if it does not exist
dir.create("figures", showWarnings = FALSE)

# Perform PCA using all 12 samples
pca_result <- prcomp(
  t(all_samples_expression),
  center = TRUE,
  scale. = FALSE
)

# Calculate variance explained
pca_variance <- (pca_result$sdev^2 /
                   sum(pca_result$sdev^2)) * 100

# Prepare plotting data
pca_data <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Genotype = sample_metadata$Genotype,
  Treatment = sample_metadata$Treatment
)

# Generate and save PCA figure
png(
  filename = file.path(figure_dir, "PCA_All_12_Samples.png"),
  width = 3000,
  height = 2200,
  res = 300
)

plot(
  pca_data$PC1,
  pca_data$PC2,
  col = as.numeric(pca_data$Genotype),
  pch = ifelse(pca_data$Treatment == "Salt", 17, 19),
  xlab = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
  ylab = paste0("PC2 (", round(pca_variance[2], 2), "%)"),
  main = "PCA of Arabidopsis Salt-Stress Transcriptomes",
  cex = 1.3
)

text(
  pca_data$PC1,
  pca_data$PC2,
  labels = pca_data$Sample,
  pos = 3,
  cex = 0.75
)

legend(
  "topright",
  legend = levels(pca_data$Genotype),
  col = seq_along(levels(pca_data$Genotype)),
  pch = 19,
  title = "Genotype"
)

legend(
  "bottomright",
  legend = c("Control", "Salt"),
  pch = c(19, 17),
  title = "Treatment"
)

# Close PNG device automatically
dev.off()

# Display variance explained
cat("PC1 variance explained:", round(pca_variance[1], 2), "%\n")
cat("PC2 variance explained:", round(pca_variance[2], 2), "%\n")

cat("\nPCA figure saved to:", file.path(figure_dir, "PCA_All_12_Samples.png"), "\n")
list.files("figures", pattern = "\\.png$")
# Generate volcano plots for all three genotype-specific salt responses

dir.create("figures", showWarnings = FALSE)

for (comparison in c(
  "Ws_Salt_Response",
  "Col_Salt_Response",
  "Colgl_Salt_Response"
)) {
  
  results <- limma::topTable(
    fit_contrasts,
    coef = comparison,
    number = Inf,
    adjust.method = "BH"
  )
  
  results$Significance <- ifelse(
    results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
    "Significant",
    "Not significant"
  )
  
  png(
    filename = file.path(
      figure_dir,
      paste0("Volcano_", comparison, ".png")
    ),
    width = 2400,
    height = 2000,
    res = 300
  )
  
  plot(
    results$logFC,
    -log10(pmax(results$adj.P.Val, .Machine$double.xmin)),
    col = ifelse(
      results$Significance == "Significant",
      "firebrick",
      "grey60"
    ),
    pch = 20,
    cex = 0.65,
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    main = paste("Volcano Plot:", comparison)
  )
  
  abline(v = c(-1, 1), lty = 2, col = "blue")
  abline(h = -log10(0.05), lty = 2, col = "blue")
  
  legend(
    "topright",
    legend = c("Significant", "Not significant"),
    col = c("firebrick", "grey60"),
    pch = 20,
    bty = "n"
  )
  
  dev.off()
}

cat("All three genotype-specific volcano plots saved successfully.\n")
list.files("figures", pattern = "\\.png$")
# Identify WRKY-associated probe sets from platform annotation

wrky_index <- grepl(
  "WRKY",
  paste(
    annotation$`Gene Symbol`,
    annotation$`Gene Title`
  ),
  ignore.case = TRUE
)

wrky_annotation <- annotation[wrky_index, ]

cat("Total WRKY-associated probe sets:",
    nrow(wrky_annotation), "\n")

# Inspect available WRKY annotations
print(
  wrky_annotation[, c(
    "ID",
    "Gene Symbol",
    "Gene Title",
    "AGI"
  )]
)
# Extract WRKY-associated differential-expression results
wrky_results_all <- list()

for (comparison in c(
  "Ws_Salt_Response",
  "Col_Salt_Response",
  "Colgl_Salt_Response"
)) {
  
  results <- limma::topTable(
    fit_contrasts,
    coef = comparison,
    number = Inf,
    adjust.method = "BH"
  )
  
  results$Probe_ID <- rownames(results)
  
  wrky_results <- results[
    results$Probe_ID %in% wrky_annotation$ID,
  ]
  
  annotation_index <- match(
    wrky_results$Probe_ID,
    annotation$ID
  )
  
  wrky_results$Gene_Symbol <-
    annotation$`Gene Symbol`[annotation_index]
  
  wrky_results$Gene_Title <-
    annotation$`Gene Title`[annotation_index]
  
  wrky_results$AGI <-
    annotation$AGI[annotation_index]
  
  wrky_results$Comparison <- comparison
  
  wrky_results_all[[comparison]] <- wrky_results
}

# Combine the three comparisons
wrky_results_combined <- do.call(
  rbind,
  wrky_results_all
)

# Save the complete WRKY results
write.csv(
  wrky_results_combined,
  "results/WRKY_All_Three_Genotypes.csv",
  row.names = FALSE
)

# Summarise significant WRKY-associated probe sets
wrky_summary <- aggregate(
  adj.P.Val ~ Comparison,
  data = wrky_results_combined,
  FUN = function(x) sum(x < 0.05, na.rm = TRUE)
)

colnames(wrky_summary)[2] <- "Significant_WRKY_Probe_Sets"

print(wrky_summary)
# Generate WRKY expression heatmap using all 12 samples

wrky_expression <- all_samples_expression[
  rownames(all_samples_expression) %in% wrky_annotation$ID,
  ,
  drop = FALSE
]

# Match annotation order to expression matrix
wrky_expression <- wrky_expression[
  match(wrky_annotation$ID, rownames(wrky_expression)),
  ,
  drop = FALSE
]

# Use gene symbols as row labels
wrky_labels <- make.unique(
  as.character(wrky_annotation$`Gene Symbol`)
)

rownames(wrky_expression) <- wrky_labels

# Standardise expression per probe set for visualisation
wrky_expression_scaled <- t(scale(t(wrky_expression)))

# Save high-resolution heatmap
png(
  file.path(figure_dir, "WRKY_Heatmap_All_12_Samples.png"),
  width = 3200,
  height = 3600,
  res = 300
)

heatmap(
  wrky_expression_scaled,
  scale = "none",
  margins = c(8, 10),
  main = "WRKY-Associated Expression Across 12 Arabidopsis Samples",
  cexRow = 0.65,
  cexCol = 0.8
)

dev.off()

cat("Updated 12-sample WRKY heatmap saved successfully.\n")
list.files(
  "figures",
  pattern = "\\.png$",
  full.names = TRUE
)
list.files(
  "results",
  pattern = "\\.csv$",
  full.names = TRUE
)
# Quality control for the five differential-expression tables

result_files <- list.files(
  "results",
  pattern = "_Full_Results\\.csv$",
  full.names = TRUE
)

qc_summary <- do.call(
  rbind,
  lapply(result_files, function(file) {
    
    data <- read.csv(file)
    
    data.frame(
      File = basename(file),
      Probe_Sets = nrow(data),
      Missing_logFC = sum(is.na(data$logFC)),
      Missing_Adjusted_P = sum(is.na(data$adj.P.Val)),
      Invalid_Adjusted_P = sum(
        !is.na(data$adj.P.Val) &
          (data$adj.P.Val < 0 | data$adj.P.Val > 1)
      )
    )
  })
)

print(qc_summary)
# Record software and package versions

cat("R version:\n")
print(R.version.string)

cat("\nlimma version:\n")
print(as.character(packageVersion("limma")))

cat("\nGEOquery version:\n")
print(as.character(packageVersion("GEOquery")))

cat("\nPlatform annotation:\n")
print(annotation$`Species Scientific Name`[1])
# Save reproducibility information

software_info <- data.frame(
  Component = c(
    "R",
    "limma",
    "GEOquery",
    "Organism",
    "Platform"
  ),
  Version_or_Details = c(
    R.version.string,
    as.character(packageVersion("limma")),
    as.character(packageVersion("GEOquery")),
    "Arabidopsis thaliana",
    "Affymetrix ATH1 Genome Array (GPL198)"
  )
)

write.csv(
  software_info,
  "results/Software_Versions.csv",
  row.names = FALSE
)

print(software_info)
browseURL(
  normalizePath(file.path(figure_dir, "PCA_All_12_Samples.png"))
)
cat(
  "ggplot2 available:", requireNamespace("ggplot2", quietly = TRUE), "\n",
  "ggrepel available:", requireNamespace("ggrepel", quietly = TRUE), "\n"
)
# Create improved PCA visualisation using ggplot2

pca_plot <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    colour = Genotype,
    shape = Treatment
  )
) +
  ggplot2::geom_point(size = 4, alpha = 0.9) +
  
  ggplot2::geom_text(
    ggplot2::aes(label = Sample),
    vjust = -1.1,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  
  ggplot2::scale_colour_manual(
    values = c(
      "Ws" = "#2878B5",
      "Col" = "#D9534F",
      "Col(gl)" = "#299D68"
    )
  ) +
  
  ggplot2::scale_shape_manual(
    values = c(
      "Control" = 16,
      "Salt" = 17
    )
  ) +
  
  ggplot2::labs(
    title = "Principal Component Analysis of Arabidopsis Salt-Stress Transcriptomes",
    x = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
    y = paste0("PC2 (", round(pca_variance[2], 2), "%)"),
    colour = "Genotype",
    shape = "Treatment"
  ) +
  
  ggplot2::theme_classic(base_size = 13) +
  
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    legend.position = "right",
    axis.title = ggplot2::element_text(face = "bold")
  )
# Final structural validation of the five experimental contrasts
stopifnot(
  nrow(contrast_matrix) == ncol(design),
  ncol(contrast_matrix) == 5,
  identical(
    colnames(contrast_matrix),
    c(
      "Ws_Salt_Response",
      "Col_Salt_Response",
      "Colgl_Salt_Response",
      "Col_vs_Ws_Interaction",
      "Colgl_vs_Ws_Interaction"
    )
  ),
  all(is.finite(contrast_matrix))
)

cat("\nAll five experimental contrasts passed structural validation.\n")
# Save the improved PCA visualisation
ggplot2::ggsave(
  filename = file.path(figure_dir, "PCA_All_12_Samples.png"),
  plot = pca_plot,
  width = 11,
  height = 8,
  dpi = 300
)

cat("Improved PCA figure saved successfully.\n")
