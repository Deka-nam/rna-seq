#!/usr/bin/env Rscript
# Complete DESeq2 analysis for catheterization study
# Works with featureCounts output

# Load required libraries
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)  # Mouse gene annotations

# 1. LOAD COUNT DATA

cat("Step 1: Loading count data...\n")

# Read featureCounts output
count_data <- read.table("counts/all_samples_counts.txt", 
                         header = TRUE, 
                         sep = "\t", 
                         comment.char = "#",
                         row.names = 1)

# Keep only count columns (remove Chr, Start, End, Strand, Length)
count_matrix <- count_data[, 6:ncol(count_data)]

# Clean up sample names (remove the .bam extension and _sorted suffix)
colnames(count_matrix) <- gsub("_sorted\\.bam$", "", colnames(count_matrix))
colnames(count_matrix) <- gsub("^bam_files/", "", colnames(count_matrix))

# Check dimensions
cat("Count matrix dimensions:", dim(count_matrix), "\n")
cat("Samples:", paste(colnames(count_matrix), collapse = ", "), "\n")

# 2. CREATE SAMPLE METADATA


# Create sample information based on sample naming
# catheterized_n1, catheterized_n2, catheterized_n3
# naive_n1, naive_n2, naive_n3

sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(c(rep("Catheterized", 3), rep("Naive", 3))),
  replicate = factor(rep(1:3, 2)),
  row.names = colnames(count_matrix)
)

cat("\nSample information:\n")
print(sample_info)


# 3. QUALITY CONTROL CHECKS


cat("\nStep 2: Performing quality checks...\n")

# Calculate basic statistics
total_counts <- colSums(count_matrix)
cat("Total counts per sample:\n")
print(total_counts)

# Check for low-count genes
genes_with_counts <- rowSums(count_matrix > 0)
cat("\nNumber of genes detected (count > 0):", sum(genes_with_counts > 0), "\n")
cat("Number of genes with zero counts:", sum(genes_with_counts == 0), "\n")

# 4. RUN DESEQ2 DIFFERENTIAL EXPRESSION


cat("\nStep 3: Running DESeq2 analysis...\n")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# Filter low-count genes (keep genes with at least 10 reads in total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# Run DESeq2
dds <- DESeq(dds)

# Results
res <- results(dds, contrast = c("condition", "Catheterized", "Naive"))
res <- res[order(res$padj), ]  # Sort by adjusted p-value

cat("\nDifferential expression summary:\n")
summary(res)


# 5. ADD GENE ANNOTATIONS

cat("\nStep 4: Adding gene annotations...\n")

# Convert to data frame
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)

# Get gene symbols from RefSeq IDs (like NM_001...)
# First, check if GeneID looks like RefSeq
if (any(grepl("^N[MR]_", res_df$GeneID[1:10]))) {
  cat("Detected RefSeq IDs, converting to gene symbols...\n")
  
  # Map RefSeq to gene symbols using org.Mm.eg.db
  library(org.Mm.eg.db)
  
  # Get mapping (RefSeq to Symbol)
  map <- select(org.Mm.eg.db, 
                keys = res_df$GeneID, 
                keytype = "REFSEQ", 
                columns = c("SYMBOL", "ENTREZID", "GENENAME"))
  
  # Merge with results
  res_df <- merge(res_df, map, by.x = "GeneID", by.y = "REFSEQ", all.x = TRUE)
} else {
  # If already gene symbols
  res_df$SYMBOL <- res_df$GeneID
}

# 6. SAVE RESULTS

cat("\nStep 5: Saving results...\n")

# Save full results
write.csv(res_df, "results/DE_results_full_Catheterized_vs_Naive.csv", row.names = FALSE)
cat("Full results saved to: results/DE_results_full_Catheterized_vs_Naive.csv\n")

# Save significant results (padj < 0.05)
sig_results <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
write.csv(sig_results, "results/DE_results_significant_padj_0.05.csv", row.names = FALSE)
cat("Significant results (padj < 0.05):", nrow(sig_results), "genes\n")

# 7. SPECIFIC ANALYSIS FOR GAG/CS GENES

cat("\nStep 6: Analyzing GAG/CS genes...\n")

# Common GAG/CS pathway genes: dated Jan 8, 2026 (expand with Vitus's list later)
gag_cs_genes <- c(
  # Chondroitin sulfate synthesis
  "Chsy1", "Chsy3", "Chpf", "Chpf2", 
  "Csgalnact1", "Csgalnact2", "Csglca1",
  # Sulfotransferases
  "Chst3", "Chst7", "Chst11", "Chst12", "Chst13", "Chst14", "Chst15",
  "Galnac4s6st", "Ust",
  # Epimerases
  "Dse", "Ds-el",
  # Other GAG-related
  "Xylt1", "Xylt2", "B4galt7", "B3galt6",
  "Ext1", "Ext2", "Extl3",  # Heparan sulfate
  "Ndst1", "Ndst2", "Ndst3", "Ndst4"
)

# Find which GAG genes are in our results
gag_in_data <- gag_cs_genes[gag_cs_genes %in% res_df$SYMBOL]
cat("GAG/CS genes found in data:", length(gag_in_data), "/", length(gag_cs_genes), "\n")

if (length(gag_in_data) > 0) {
  gag_results <- res_df[res_df$SYMBOL %in% gag_in_data, ]
  
  # Sort by significance
  gag_results <- gag_results[order(gag_results$padj, na.last = TRUE), ]
  
  write.csv(gag_results, "results/DE_results_GAG_CS_genes.csv", row.names = FALSE)
  cat("GAG/CS results saved to: results/DE_results_GAG_CS_genes.csv\n")
  
  # Print top GAG genes
  cat("\nTop 10 GAG/CS genes by significance:\n")
  print(gag_results[1:min(10, nrow(gag_results)), 
                    c("SYMBOL", "baseMean", "log2FoldChange", "padj")])
}


# 8. CREATE VISUALIZATIONS

cat("\nStep 7: Creating visualizations...\n")

# Create output directory for plots
dir.create("results/plots", showWarnings = FALSE)

# 8.1 Volcano plot
volcano_data <- res_df
volcano_data$significant <- !is.na(volcano_data$padj) & volcano_data$padj < 0.05
volcano_data$color <- ifelse(volcano_data$significant, "red", "gray")
volcano_data$color <- ifelse(volcano_data$SYMBOL %in% gag_cs_genes, "blue", volcano_data$color)

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = color), alpha = 0.6, size = 1) +
  scale_color_identity() +
  theme_minimal() +
  labs(title = "Volcano Plot: Catheterized vs Naive",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)

ggsave("results/plots/volcano_plot.png", p_volcano, width = 8, height = 6)

# 8.2 Heatmap of top genes
top_genes <- head(res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ], 50)
if (nrow(top_genes) > 5) {
  # Get normalized counts
  vsd <- vst(dds, blind = FALSE)
  top_counts <- assay(vsd)[top_genes$GeneID, ]
  
  # Create annotation for samples
  annotation_col <- data.frame(
    Condition = sample_info$condition,
    row.names = colnames(top_counts)
  )
  
  png("results/plots/heatmap_top_genes.png", width = 800, height = 1000)
  pheatmap(top_counts,
           annotation_col = annotation_col,
           show_rownames = TRUE,
           scale = "row",
           main = "Top 50 DEGs: Catheterized vs Naive")
  dev.off()
}

# 8.3 MA plot
png("results/plots/MA_plot.png", width = 800, height = 600)
plotMA(res, main = "MA Plot: Catheterized vs Naive", ylim = c(-5, 5))
dev.off()


# 9. CREATE SUMMARY REPORT

cat("\nStep 8: Generating summary report...\n")

# Create a summary text file
sink("results/analysis_summary.txt")
cat("RNA-Seq Analysis Summary\n")
cat("=======================\n\n")
cat("Project: Bladder Catheterization in B6 Mice\n")
cat("Date:", date(), "\n\n")

cat("Sample Information:\n")
print(sample_info)
cat("\n")

cat("Read Count Summary:\n")
cat("Total genes analyzed:", nrow(dds), "\n")
cat("Significant DEGs (padj < 0.05):", sum(!is.na(res$padj) & res$padj < 0.05), "\n")
cat("Upregulated in Catheterized (padj < 0.05 & LFC > 1):", 
    sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1), "\n")
cat("Downregulated in Catheterized (padj < 0.05 & LFC < -1):", 
    sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1), "\n\n")

cat("Top 10 DEGs:\n")
top10 <- head(res_df[!is.na(res_df$padj), ], 10)
print(top10[, c("SYMBOL", "log2FoldChange", "padj")])
cat("\n")

if (exists("gag_results") && nrow(gag_results) > 0) {
  cat("GAG/CS Genes Analysis:\n")
  cat("Total GAG/CS genes tested:", length(gag_in_data), "\n")
  cat("Significantly changed GAG/CS genes (padj < 0.05):", 
      sum(gag_results$padj < 0.05, na.rm = TRUE), "\n")
  
  sig_gag <- gag_results[gag_results$padj < 0.05, ]
  if (nrow(sig_gag) > 0) {
    cat("\nSignificant GAG/CS genes:\n")
    print(sig_gag[, c("SYMBOL", "log2FoldChange", "padj")])
  }
}
sink()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Check the 'results' directory for:\n")
cat("1. DE_results_full_Catheterized_vs_Naive.csv - All genes\n")
cat("2. DE_results_GAG_CS_genes.csv - Your genes of interest\n")
cat("3. Various plots in results/plots/\n")
cat("4. analysis_summary.txt - Summary report\n")