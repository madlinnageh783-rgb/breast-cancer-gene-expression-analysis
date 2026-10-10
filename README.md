Breast Cancer Gene Expression Analysis
Project Overview
This project analyzes gene expression data from breast cancer patients to identify probes associated with relapse (recurrence of cancer).
The analysis uses the GSE2034 dataset from the Gene Expression Omnibus (GEO) database.
Steps of Analysis
1. Data collection from the GEO database (GSE2034).
2. Comparison of gene expression between patients with and without relapse.
3. Identification of probes with statistically significant differences using FDR correction.
4. Visualization of results using:
o Volcano Plot
o Heatmap of the top 20 probes
o Principal Component Analysis (PCA)

Results
* A total of 916 probes showed statistically significant differences in gene expression between patients with and without relapse (FDR < 0.05).
* The top 20 probes with the lowest FDR values were selected for focused visualization and further examination.
* A Volcano Plot was used to visualize differences in gene expression and statistical significance.
* A Heatmap was generated to explore expression patterns of the top 20 probes across samples.
* Principal Component Analysis (PCA) was performed using the 500 most variable probes to explore sample-level variation and potential separation between relapse groups.
Tools Used
* R
* RStudio
* GEO database (GSE2034)
* ggplot2
* pheatmap
* Bioconductor
