# BiocManager install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# install Bioconductor packages
packages <- c("tximport", "GenomicFeatures", "DESeq2", "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "EnhancedVolcano", "AnnotationDbi", "ReactomePA", "GEOquery")
BiocManager::install(packages)


# install CRAN package
install.packages(c("readr", "ggplot2", "msigdbr"))
install.packages(packages)





