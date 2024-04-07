#!/usr/bin/Rscript
### RNASeq Analysis

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Load the DESeq2 package
library(DESeq2)
library(dplyr)
# Read the count data into R
count_data <- read.csv("./data/count_datat.csv", header = TRUE, row.names = 1)
head(count_data)

# Read the sample metadata
metadata <- read.csv("./data/meta_data.csv", header = TRUE, row.names = 1)

##match the row and column names of both datasets
all(colnames(count_data) %in% as.character(rownames(metadata)))

## they are in same order
all(colnames(count_data) == (rownames(metadata)))


# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              design = ~ condition)
dds

# Pre-filtering (optional)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds

##mention which is control and traeted 
dds$condition <- relevel(dds$condition, ref = "control")
dds$condition

# Differential expression analysis
dds <- DESeq(dds)

# Get differential expression results
results <- results(dds)

# Explore results (optional)
head(results)
results
summary(results)
# Write the significant results to a CSV file
write.csv(results, file = "results.csv", row.names = TRUE)


# Filter results for significant genes (optional)
results_sig <- results[which(results$padj < 0.05), ]
results_sig
summary(results_sig)
# Write the significant results to a CSV file
write.csv(results_sig, file = "significant_results.csv", row.names = TRUE)


## Visualiz the results

library(ggplot2)

# NA Plot
ma_plot <- ggplot(results, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "MA Plot",
       x = "Mean of normalized counts (baseMean)",
       y = "Log2 fold change") +
  theme_minimal()

# Show plots
print(ma_plot)


# Volcano Plot
volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "signigficant", "Non-significant")), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Significance") +
  theme_minimal()

# Show plots
print(volcano_plot)


# Volcano Plot with Custom Colors
volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, 
                                ifelse(log2FoldChange > 0, "upregulated", "downregulated"), 
                                "Non-significant")), alpha = 0.7) +
  scale_color_manual(values = c("red", "grey", "blue")) +  # Custom colors for each category
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Significance") +
  theme_minimal()

# Show plots
print(volcano_plot)


# Sort results by absolute log2 fold change
sorted_results <- results[order(abs(results$log2FoldChange), decreasing = TRUE),]

# Identify top three upregulated and downregulated genes
top_upregulated <- head(sorted_results[sorted_results$log2FoldChange > 0,], 5)
top_downregulated <- head(sorted_results[sorted_results$log2FoldChange < 0,], 5)

# Volcano Plot with Custom Colors and Labels
volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, 
                                ifelse(log2FoldChange > 0, "upregulated", "downregulated"), 
                                "Non-significant")), alpha = 0.7) +
  scale_color_manual(values = c("red", "blue", "grey")) +  # Custom colors for each category
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_text_repel(data = rbind(top_upregulated, top_downregulated), 
                  aes(label = rownames(rbind(top_upregulated, top_downregulated))), 
                  size = 3, nudge_y = 0.2, color = "black") +  # Label top up/downregulated genes
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Significance") +
  theme_minimal()

# Show plots
print(volcano_plot)


library(ggrepel)  # Load ggrepel for geom_text_repel function

# Your data and plot code

# Volcano Plot with Custom Colors and Labels
volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, 
                                ifelse(log2FoldChange > 0, "upregulated", "downregulated"), 
                                "Non-significant")), alpha = 0.7) +
  scale_color_manual(values = c("red", "grey", "blue")) +  # Custom colors for each category
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_text_repel(data = rbind(top_upregulated, top_downregulated), 
                  aes(label = rownames(rbind(top_upregulated, top_downregulated))), 
                  size = 3, nudge_y = 0.2, color = "black") +  # Label top up/downregulated genes
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Significance") +
  theme_minimal()

# Show plots
print(volcano_plot)


