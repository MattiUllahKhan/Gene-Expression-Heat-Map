# Install and load necessary libraries
##BiocManager::install("GEOquery")
##install.packages("pheatmap")

# Load libraries
library(GEOquery)
library(pheatmap)

##install.packages("data.table")
library(data.table)



# Set the path to the downloaded file
file_path <- "GSE10245_series_matrix.txt.gz"

# Load the data
geo_data <- getGEO(filename = file_path)
expression_data <- exprs(geo_data)

# Convert the expression data to a data frame
expression_data_df <- as.data.frame(expression_data)

# Save the expression data as a CSV file
write.csv(expression_data_df, "expression_data.csv", row.names = TRUE)

# Confirm the file is saved
print("Expression data saved to expression_data.csv")

# Display first few rows of the dataset
head(expression_data)

# Log2 transform the data if necessary
expression_data <- log2(expression_data + 1)

# Convert the expression data to a data.table
expression_data_dt <- as.data.table(expression_data)

# Perform the same filtering steps as before
gene_variances <- expression_data_dt[, apply(.SD, 1, var)]
top_genes <- order(gene_variances, decreasing = TRUE)[1:1000]
expression_data_subset <- expression_data_dt[top_genes, ]

# Create a heatmap with the subset of data
pheatmap(as.matrix(expression_data_subset), scale = "row", show_rownames = TRUE, show_colnames = TRUE,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         clustering_method = "complete", color = colorRampPalette(c("navy", "white", "firebrick3"))(50))














# Create a heatmap with standardized data
pheatmap(expression_data_standardized,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = TRUE,
         main = "Gene Expression Heatmap 2nd",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# Create a heatmap with more color stops
pheatmap(expression_data_normalized,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = TRUE,
         main = "Gene Expression Heatmap 2nd",
         color = colorRampPalette(c("navy", "blue", "white", "orange", "firebrick3"))(50))



# Create a heatmap with more color stops
pheatmap(expression_data_log, scale = "column",
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,cutree_rows = 4, cutree_cols = 2,
         main = "Gene Expression Heatmap 2nd",
         color = colorRampPalette(c("navy", "blue", "pink", "orange", "firebrick3"))(50))


# Normalize the data (optional)
expression_data_normalized <- scale(expression_data)


# Apply log transformation
expression_data_log <- log1p(expression_data_normalized)

# Create a heatmap with log-transformed data
pheatmap(expression_data_log,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = TRUE,
         main = "Gene Expression Heatmap 2nd",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))



# Set the row names to gene symbols for better visualization
rownames(expression_data_normalized) <- data$gene_symbol

# Create a heatmap
pheatmap(expression_data_normalized, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = TRUE,
         main = "Gene Expression Heatmap 2nd", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

