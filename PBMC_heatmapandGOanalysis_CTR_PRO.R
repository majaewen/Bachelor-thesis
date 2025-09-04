"R Script for the analysis of significant proteins at T2 in PBMC samples separately in PRO and CTR.
Output: Heatmap and GO analysis"


install.packages(c(
  "tidyverse",
  "writexl",
  "readxl",
  "pheatmap",
  "openxlsx"
))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "DOSE",
  "enrichplot"
))
library(tidyverse)
library(writexl)
library(readxl)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(dplyr)
library(stringr)
library(openxlsx)

###CTR###
###CTR###
###CTR###
data <- read_xls("/Users/majaewen/Documents/Uni/6_Semester/Daten/PBMC_proteomics/T_test_Controls.xls", sheet = 1, range = "A1:CH8311")
data_ordered <- data

# sort by subject and time
timepoints <- as.character(data_ordered[2, ])
subjects   <- as.character(data_ordered[4, ])
time_num <- as.integer(sub("T", "", timepoints)) 
order_idx <- order(time_num, subjects)
data_ordered <- data[, order_idx]

new_colnames <- paste0(subjects[order_idx], timepoints[order_idx])
colnames(data_ordered) <- new_colnames

data_ordered <- data_ordered[ , c(1:61, 84)]
colnames(data_ordered)[61] <- "Welch's T-test significant"
names(data_ordered) <- make.unique(names(data_ordered))


# filter rows with significant regulation at T2"
data_sig_rows <- data_ordered %>% filter(`Welch's T-test significant` == "+")
head(data_sig_rows)

data_sig_df <- as.data.frame(data_sig_rows)
rownames(data_sig_df) <- data_sig_rows[[62]]
data_sig_df <- data_sig_df[, -c(61,62)]
colnames(data_sig_df) <- new_colnames[1:60] 


# perform corrected Z-score normalization for each row
data_normalized <- t(apply(data_sig_df, 1, function(x) {
  x_numeric <- as.numeric(x)
  
  if (all(is.na(x_numeric))) return(rep(NA, length(x_numeric)))
  if (length(unique(na.omit(x_numeric))) == 1) return(rep(0, length(x_numeric)))
  
  return(as.numeric(scale(x_numeric, center = TRUE, scale = TRUE)))
}))

data_normalized <- as.data.frame(data_normalized)
colnames(data_normalized) <- colnames(data_sig_df)


# save normalized data
write_xlsx(data_normalized, path = "ProteinData_Normalized_CTR.xlsx")
write_delim(data_normalized, "ProteinData_Normalized_CTR.txt", delim = "\t")


# generate heatmap 
pheatmap_result <- pheatmap(
  as.matrix(data_normalized),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "average",
  color = colorRampPalette(c("royalblue3", "white", "firebrick3"))(200),
  breaks = seq(-3, 3, length.out = 201),
  show_rownames = FALSE
)

# extract cluster assignments 
clusters <- cutree(pheatmap_result$tree_row, k = 2)
data_sig_df$Cluster <- factor(clusters)
data_sig_df$ID <- rownames(data_normalized)

annotation_row <- data.frame(Cluster = data_sig_df$Cluster)
rownames(annotation_row) <- rownames(data_normalized)

ann_colors <- list(
  Timepoint = c(
    T0 = "#deebf7",
    T1 = "#c6dbef",
    T2 = "#9ecae1",
    T3 = "#6baed6",
    T4 = "#3182bd"
  ),
  Cluster = c(
    "1" = "#66C2A5",  
    "2" = "#FFD92F"  
  )
)
timepoints <- gsub(".*T", "T", colnames(data_normalized))
annotation_col <- data.frame(Timepoint = timepoints)
rownames(annotation_col) <- colnames(data_normalized)

pheatmap_result <- pheatmap(
  as.matrix(data_normalized),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "average",
  annotation_col = annotation_col,  # Apply timepoint annotation
  annotation_colors = ann_colors,
  color = colorRampPalette(c("royalblue3", "white", "firebrick3"))(200),
  breaks = seq(-3, 3, length.out = 201),
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  annotation_row = annotation_row,
  main = "CTR"
)

# extract gene symbols
gene_column_name <- "ID"
cluster1_genes <- unique(na.omit(data_sig_df$ID[data_sig_df$Cluster == 1]))
cluster2_genes <- unique(na.omit(data_sig_df$ID[data_sig_df$Cluster == 2]))
background_genes <- unique(na.omit(data[["PG.Genes"]]))

cat("Total unique background genes:", length(background_genes), "\n")

cluster1_entrez <- bitr(cluster1_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>% drop_na(ENTREZID) #all gene names from Cluster 1, unique, without NA
cluster2_entrez <- bitr(cluster2_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>% drop_na(ENTREZID)
background_entrez <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# cluster1: 37
# cluster2: 1


# GO enrichment analysis (only Cluster 1, adapt for Cluster 2)
go_enrich_cluster1 <- enrichGO(
  gene = cluster1_entrez$ENTREZID,
  universe = background_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# save GO enrichment results
write.xlsx(
  as.data.frame(go_enrich_cluster1),
  file = "/Users/majaewen/Documents/Uni/6_Semester/Daten/PBMC_proteomics/Cluster1_BP_results_CTR.xlsx",
  rowNames = FALSE
)


# simplify results 
simplified_results <- simplify(
  go_enrich_cluster1,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)


# generate barplot
p <- barplot(simplified_results, showCategory = 5)
p + ggtitle("Cluster 1")



###PRO###
###PRO###
###PRO###

data <- read_xls("/Users/majaewen/Documents/Uni/6_Semester/Daten/PBMC_proteomics/T_test_Profis.xls", sheet = 1, range = "A1:CG8311")
data_ordered <- data

#sort by subject and time
timepoints <- as.character(data_ordered[2, ])
subjects   <- as.character(data_ordered[4, ])
time_num <- as.integer(sub("T", "", timepoints)) 
order_idx <- order(time_num, subjects)
data_ordered <- data[, order_idx]

new_colnames <- paste0(subjects[order_idx], timepoints[order_idx])
colnames(data_ordered) <- new_colnames

data_ordered <- data_ordered[ , -c(61:64, 66:83, 85)]
colnames(data_ordered)[61] <- "Welch's T-test significant"

names(data_ordered) <- make.unique(names(data_ordered))

#Filter rows with significant regulation at T2"
data_sig_rows <- data_ordered %>% filter(`Welch's T-test significant` == "T2_T0")
head(data_sig_rows)

data_sig_df <- as.data.frame(data_sig_rows)
rownames(data_sig_df) <- data_sig_rows[[62]]
data_sig_df <- data_sig_df[, -c(61,62)]
colnames(data_sig_df) <- new_colnames[1:60] 

#perform corrected Z-score normalization for each row
data_normalized <- t(apply(data_sig_df, 1, function(x) {
  x_numeric <- as.numeric(x)
  
  if (all(is.na(x_numeric))) return(rep(NA, length(x_numeric)))
  if (length(unique(na.omit(x_numeric))) == 1) return(rep(0, length(x_numeric)))
  
  return(as.numeric(scale(x_numeric, center = TRUE, scale = TRUE)))
}))

data_normalized <- as.data.frame(data_normalized)
colnames(data_normalized) <- colnames(data_sig_df)

# Save normalized data
write_xlsx(data_normalized, path = "ProteinData_Normalized.xlsx")
write_delim(data_normalized, "ProteinData_Normalized.txt", delim = "\t")

# Generate heatmap
pheatmap_result <- pheatmap(
  as.matrix(data_normalized),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "average",
  color = colorRampPalette(c("royalblue3", "white", "firebrick3"))(200),
  breaks = seq(-3, 3, length.out = 201),
  show_rownames = FALSE
)


# Extract cluster assignments 
clusters <- cutree(pheatmap_result$tree_row, k = 2)

clusters <- ifelse(clusters == 1, 2, 1)

data_sig_df$Cluster <- factor(clusters)
data_sig_df$ID <- rownames(data_normalized)

annotation_row <- data.frame(Cluster = data_sig_df$Cluster)
rownames(annotation_row) <- rownames(data_normalized)

ann_colors <- list(
  Timepoint = c(
    T0 = "#deebf7",
    T1 = "#c6dbef",
    T2 = "#9ecae1",
    T3 = "#6baed6",
    T4 = "#3182bd"
  ),
  Cluster = c(
    "1" = "#66C2A5",  
    "2" = "#FFD92F"   
  )
)

timepoints <- gsub(".*T", "T", colnames(data_normalized))
annotation_col <- data.frame(Timepoint = timepoints)
rownames(annotation_col) <- colnames(data_normalized)

pheatmap_result <- pheatmap(
  as.matrix(data_normalized),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "average",
  annotation_col = annotation_col,  # Apply timepoint annotation
  annotation_colors = ann_colors,
  color = colorRampPalette(c("royalblue3", "white", "firebrick3"))(201),
  breaks = seq(-3, 3, length.out = 201),
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  main = "PRO"
)

# extract gene symbols
gene_column_name <- "ID"
cluster1_genes <- unique(na.omit(data_sig_df$ID[data_sig_df$Cluster == 1]))
cluster2_genes <- unique(na.omit(data_sig_df$ID[data_sig_df$Cluster == 2]))
background_genes <- unique(na.omit(data[["PG.Genes"]]))

cat("Total unique background genes:", length(background_genes), "\n")

cluster1_entrez <- bitr(cluster1_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>% drop_na(ENTREZID) #all gene names from Cluster 1, unique, without NA
cluster2_entrez <- bitr(cluster2_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>% drop_na(ENTREZID)
background_entrez <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


# GO enrichment analysis (only Cluster 1, adapt for Cluster 2)
go_enrich_cluster1 <- enrichGO(
  gene = cluster1_entrez$ENTREZID,
  universe = background_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Save GO enrichment results
write.xlsx(
  as.data.frame(go_enrich_cluster2),
  file = "/Users/majaewen/Documents/Uni/6_Semester/Daten/PBMC_proteomics/Cluster2_MF_results.xlsx",
  rowNames = FALSE
)

# Simplify GO results
simplified_results <- simplify(
  go_enrich_cluster1,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# generate barplot, adapt font size
p <- barplot(simplified_results, showCategory = 5)
p + ggtitle("Cluster 2 - PRO") + theme(
  plot.title = element_text(size = 20, face = "bold"),
  axis.text.y = element_text(size = 15),              
  axis.text.x = element_text(size = 15),   
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16)
)

#see all GOs
GOsPRO <- as.data.frame(go_enrich_cluster1)