install.packages("tidyverse") 
library(tidyverse)
R.version.string
install.packages("tidyverse")
library(tidyverse)
install.packages("BiocManager")
install("GEOquery")
library(GEOquery)

gse <- getGEO("GSE2034", GSEMatrix = TRUE)
dim(gse[[1]])
head(exprs(gse[[1]]))
pData(gse[[1]])[1:5, 1:5]
colnames(pData(gse[[1]]))
table(pData(gse[[1]])$`bone relapses (1=yes, 0=no):ch1`)
expr_data <- exprs(gse[[1]])
relapse <- pData(gse[[1]])$`bone relapses (1=yes, 0=no):ch1`
table(relapse)
mean(expr_data[1, relapse == 0])
mean(expr_data[1, relapse == 1])
t.test(expr_data[1, ] ~ relapse)
p_values <- apply(expr_data, 1, function(x) {
  t.test(x ~ relapse)$p.value
})
head(sort(p_values))
p_adj <- p.adjust(p_values, method = "BH")
head(sort(p_adj))
sum(p_adj < 0.05)
top_features <- head(sort(p_adj), 20)
top_features
feature_data <- fData(gse[[1]])
feature_data <- fData(gse[[1]])
colnames(feature_data)
top_genes <- feature_data[top_features, "Gene Symbol"]
top_genes
head(rownames(feature_data))
top_genes <- feature_data[match(names(top_features), feature_data$ID), "Gene Symbol"]
top_genes
top_table <- data.frame(
  Probe_ID = names(top_features),
  Gene = top_genes,
  FDR = as.numeric(top_features)
)
top_table
library(ggplot2)
results <- data.frame(
  Probe_ID = names(p_values),
  P_value = p_values,
  FDR = p_adj
)
results$Mean_No_Relapse <- apply(expr_data[, relapse == 0], 1, mean)
results$Mean_Relapse <- apply(expr_data[, relapse == 1], 1, mean)
results$Log2_FC <- log2(results$Mean_Relapse / results$Mean_No_Relapse)
results$neg_log10_FDR <- -log10(results$FDR)
ggplot(results, aes(x = Log2_FC, y = neg_log10_FDR)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10(FDR)"
  )
top_results <- results[order(results$FDR), ][1:20, ]
top_expr <- expr_data[top_results$Probe_ID, ]

top_results
top_results$Gene <- feature_data[
  match(top_results$Probe_ID, feature_data$ID),
  "Gene Symbol"
]

top_results
top_expr <- expr_data[top_results$Probe_ID, ]

dim(top_expr)
rownames(top_expr) <- top_results$Gene

head(top_expr)
install.packages("pheatmap")
library(pheatmap)
top_expr_scaled <- t(scale(t(top_expr)))
annotation_col <- data.frame(
  Relapse = factor(relapse, levels = c(0, 1),
                   labels = c("No Relapse", "Relapse"))
)
annotation_col <- data.frame(
  Relapse = factor(
    relapse,
    levels = c(0, 1),
    labels = c("No Relapse", "Relapse")
  )
)
library(pheatmap)

pheatmap(
  top_expr_scaled,
  annotation_col = annotation_col,
  main = "Top 20 Differentially Expressed Probes"
)
rownames(annotation_col) <- colnames(top_expr_scaled)
rownames(annotation_col) <- colnames(top_expr_scaled)
pheatmap
rownames(top_expr)
rownames(top_expr) <- paste(top_results$Gene, top_results$Probe_ID, sep = " | ")
rownames(top_expr)
pheatmap(
  top_expr_scaled,
  annotation_col = annotation_col,
  main = "Top 20 Differentially Expressed Probes"
)
rownames(top_expr) <- top_results$Probe_ID
pca <- prcomp(t(expr_data), scale. = TRUE)
library(GEOquery)
ls()
probe_var <- apply(expr_data, 1, var)

top_variable <- names(
  sort(probe_var, decreasing = TRUE)
)[1:500]

pca <- prcomp(
  t(expr_data[top_variable, ]),
  scale. = TRUE
)
probe_var <- apply(expr_data, 1, var)
top_variable <- names(
  sort(probe_var, decreasing = TRUE)
)[1:500]
pca <- prcomp(
  t(expr_data[top_variable, ]),
  scale. = TRUE
)
pca <- prcomp(
  t(expr_data[top_variable, ]),
  scale. = TRUE
)
pca <- prcomp(
  t(expr_data[top_variable, ]),
  scale. = TRUE
)
pca <- prcomp(
  t(expr_data[top_variable, ]),
  scale. = TRUE
)
pca_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Relapse = factor(
    relapse,
    levels = c(0, 1),
    labels = c("No Relapse", "Relapse")
  )
)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Relapse)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "PCA of Breast Cancer Samples",
    x = "PC1",
    y = "PC2",
    color = "Relapse Status"
  )
summary(pca)
colnames(top_genes)
colnames(feature_data)
top_genes <- merge(
  top_results,
  feature_data,
  by.x = "Probe_ID",
  by.y = "ID",
  all.x = TRUE
)
head(top_genes[, c("Probe_ID", "Gene Symbol", "Gene Title", "FDR")])
top_genes[, c(
  "Gene Symbol",
  "Gene Title",
  "Mean_No_Relapse",
  "Mean_Relapse",
  "Log2_FC",
  "FDR"
)]
write.csv(
  top_genes,
  "top_genes_results.csv",
  row.names = FALSE
)
write.csv(
  top_genes,
  "top_genes_results.csv",
  row.names = FALSE
)
file.exists("top_genes_results.csv")