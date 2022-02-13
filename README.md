# SeuratAnalysis_Shiny
R code for a shiny app to run the standard Seurat workflow with consequent data visualization.

Can be used for scRNA datasets with different treatments. 
1. Quality Control:
Loaded data is visualized and QC threshholds can be set based on: Number of Counts, number of features, mitochondrial percentage and ribosomal percentage.

2. SCTransform:
Data is Normalized and scaled using SCTransform, variable Features are found.
https://satijalab.org/seurat/articles/sctransform_vignette.html

3. PCA
PCA is run for dimensional reduction.

4. Clustering
Clusters are then found using the "FindClusters()" function (https://www.rdocumentation.org/packages/Seurat/versions/4.1.0/topics/FindClusters).

5. Visualization
A variety of plots can be used to analyze the data, including several featurePlots, Violin Plots, a Heatmap, and a Scatterplot (to compare treatments).
To save plots: Right Click and save.
