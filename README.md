
# RNA_seq_analysis_for_ALS

- [My Project](file:///C:/Users/ramis/OneDrive/Desktop/BIOINFORMATICS%20202/RNA_seq_analysis/RNA_seq_analysis.html)

- [The Project](file:///C:/Users/ramis/OneDrive/Desktop/bioinformatics%201/Bioinformatics1_2025-20250716T162722Z-1-001/Bioinformatics1_2025/EXAM%20project/project_RNAseq_2025.html#Software_to_install_in_your_conda_environment)

## Disease sample

ALS( Amyotrophic Lateral Sclerosis) is a fatal neurodegenerative disorder in which motor neurons progressively die leading to muscle weakness, paralysis, and eventually breathing failure
C9ORF72 repeat which is expansions of a hexanucleotide (GGGGCC) repeat in the C9ORF72 gene. It is the most common genetic cause, linking ALS to toxic RNA and protein mechanisms.
Cognitive functions are usually preserved, but in some patients, frontotemporal dementia (FTD) co-occurs.

## The Experiment
[Targeting RNA foci in iPSC-derived motor neurons from ALS patients with C9ORF72 repeat expansion](https://www.ebi.ac.uk/gxa/experiments/E-GEOD-52202/Results?specific=false&geneQuery=%255B%255D&filterFactors=%257B%257D&cutoff=%257B%2522foldChange%2522%253A1%252C%2522pValue%2522%253A0.05%257D&regulation=%2522UP%2522)

We Use iPSC-derived motor neurons from C9ORF72-ALS patients to model the disease. Determine if ALS pathology is due to:
- Loss of function (reduced C9ORF72 expression)
- Gain of function (toxic RNA/protein from repeat expansions)
  
RNA-binding proteins (RBPs) are crucial for RNA processing (splicing, transport, stability, translation, localization)

In C9ORF72-ALS, repeat RNAs form RNA foci that trap RBPs, leading to Mis-splicing of many genes,Disrupted RNA transport and Neuronal dysfunction.

Gene Expression also changes Patient motor neurons show altered expression of genes controlling membrane excitability (e.g., DPP6, KCNQ3). These neurons have reduced ability to fire continuous action potentials, impairing communication with muscles and contributing to ALS symptoms. * ASO Treatment Antisense oligonucleotides (ASOs) targeting C9ORF72 RNA suppressed RNA foci, partially reversed abnormal gene expression
C9ORF72-ALS is mainly caused by toxic gain-of-function RNA (not just loss of gene expression). ASOs are a promising therapy to block toxic RNA effects and restore normal gene expression in ALS patient neurons.

## DESeq2
Differential Expression (DE) analysis using DESeq2 It’s the process of finding which genes are expressed differently between experimental conditions (e.g., treatment vs control, diseased vs healthy, knockout vs wild-type).

```
load library
suppressMessages(suppressWarnings(library(DESeq2)))
```

## The read count table:

```counts <- read.table('../RNA_seq_analysis/data/E-GEOD-52202-raw-counts.txt', header = TRUE, row.names = 1, sep = "\t")  
#remove the "Gene.Name" column
countData <- as.matrix(subset(counts, select = c(-Gene.Name)))
head(countData,10) #it is a numerical matrix
```
```
##                 SRR1027602 SRR1027592 SRR1027598 SRR1027591 SRR1027600
## ENSG00000000003       1141       1067       1119       1111        800
## ENSG00000000005          7          3        105          2          6
## ENSG00000000419        226        269        273        264        191
## ENSG00000000457        160        244        127        254        160
## ENSG00000000460         89        158         89        162         68
## ENSG00000000938          0          0          0          0          0
## ENSG00000000971          7          8         11         13          7
## ENSG00000001036        275        492        406        491        317
## ENSG00000001084        383        591        580        580        349
## ENSG00000001167        605        666        603        693        525
##                 SRR1027596 SRR1027601 SRR1027599 SRR1027605 SRR1027604
## ENSG00000000003       1209       1115        822       1104        660
## ENSG00000000005          8          2          5          1          1
## ENSG00000000419        241        231        220        388        218
## ENSG00000000457        184        184        180        242        128
## ENSG00000000460         80         80         60        105         52
## ENSG00000000938          0          0          0          0          0
## ENSG00000000971         33          3          4         16          1
## ENSG00000001036        335        254        274        452        348
## ENSG00000001084        442        340        375        529        500
## ENSG00000001167        607        609        531        850        852
##                 SRR1027595 SRR1027594 SRR1027593 SRR1027597 SRR1027603
## ENSG00000000003       1076        754        763       1116        659
## ENSG00000000005         20          3          2        106          2
## ENSG00000000419        232        226        227        241        257
## ENSG00000000457        166         96         83        141        132
## ENSG00000000460         72         49         54         86         56
## ENSG00000000938          0          2          1          0          0
## ENSG00000000971         24          7          3         11          0
## ENSG00000001036        327        547        451        403        346
## ENSG00000001084        444        407        366        582        504
## ENSG00000001167        588        905        855        598        913
##                 SRR1027606
## ENSG00000000003       1112
## ENSG00000000005          2
## ENSG00000000419        400
## ENSG00000000457        255
## ENSG00000000460         88
## ENSG00000000938          0
## ENSG00000000971         12
## ENSG00000001036        513
## ENSG00000001084        524
## ENSG00000001167        809
dim(countData)
## [1] 58735    16
summary(countData)
##    SRR1027602        SRR1027592        SRR1027598         SRR1027591     
##  Min.   :    0.0   Min.   :    0.0   Min.   :     0.0   Min.   :    0.0  
##  1st Qu.:    0.0   1st Qu.:    0.0   1st Qu.:     0.0   1st Qu.:    0.0  
##  Median :    0.0   Median :    0.0   Median :     0.0   Median :    0.0  
##  Mean   :  168.7   Mean   :  201.9   Mean   :   187.6   Mean   :  202.9  
##  3rd Qu.:   18.0   3rd Qu.:   22.0   3rd Qu.:    19.5   3rd Qu.:   23.0  
##  Max.   :76265.0   Max.   :92349.0   Max.   :169567.0   Max.   :92420.0  
##    SRR1027600      SRR1027596      SRR1027601        SRR1027599     
##  Min.   :    0   Min.   :    0   Min.   :    0.0   Min.   :    0.0  
##  1st Qu.:    0   1st Qu.:    0   1st Qu.:    0.0   1st Qu.:    0.0  
##  Median :    0   Median :    0   Median :    0.0   Median :    0.0  
##  Mean   :  150   Mean   :  181   Mean   :  169.8   Mean   :  151.2  
##  3rd Qu.:   16   3rd Qu.:   20   3rd Qu.:   17.0   3rd Qu.:   16.0  
##  Max.   :60530   Max.   :88405   Max.   :77346.0   Max.   :60695.0  
##    SRR1027605        SRR1027604         SRR1027595        SRR1027594      
##  Min.   :    0.0   Min.   :     0.0   Min.   :    0.0   Min.   :     0.0  
##  1st Qu.:    0.0   1st Qu.:     0.0   1st Qu.:    0.0   1st Qu.:     0.0  
##  Median :    0.0   Median :     0.0   Median :    0.0   Median :     0.0  
##  Mean   :  219.1   Mean   :   229.8   Mean   :  182.2   Mean   :   192.5  
##  3rd Qu.:   23.0   3rd Qu.:    23.0   3rd Qu.:   20.0   3rd Qu.:    19.0  
##  Max.   :64512.0   Max.   :135382.0   Max.   :89332.0   Max.   :140025.0  
##    SRR1027593         SRR1027597         SRR1027603         SRR1027606     
##  Min.   :     0.0   Min.   :     0.0   Min.   :     0.0   Min.   :    0.0  
##  1st Qu.:     0.0   1st Qu.:     0.0   1st Qu.:     0.0   1st Qu.:    0.0  
##  Median :     0.0   Median :     0.0   Median :     0.0   Median :    0.0  
##  Mean   :   193.4   Mean   :   188.7   Mean   :   231.8   Mean   :  217.6  
##  3rd Qu.:    18.0   3rd Qu.:    19.0   3rd Qu.:    23.0   3rd Qu.:   23.0  
##  Max.   :140136.0   Max.   :172097.0   Max.   :135697.0   Max.   :64054.0
```

## A colData table:
it’s essentially for the experimental setup that tells DESeq2 how to compare samples. DESeq2 uses colData together with the design formula to know what comparisons to make.

```
# read in sample info
colData <- read.table("../RNA_seq_analysis/data/E-GEOD-52202-experiment-design.tsv", header = TRUE, row.names=1 , sep = '\t', 
                      stringsAsFactors = TRUE)
colData<-colData[,c(1,3,5,7,9,11)] #clean the data frame
colData
```

```
##            Sample.Characteristic.cell.type. Sample.Characteristic.disease.
## SRR1027591                     motor neuron                         normal
## SRR1027592                     motor neuron                         normal
## SRR1027593                     motor neuron                         normal
## SRR1027594                     motor neuron                         normal
## SRR1027595                     motor neuron                         normal
## SRR1027596                     motor neuron                         normal
## SRR1027597                     motor neuron                         normal
## SRR1027598                     motor neuron                         normal
## SRR1027599                     motor neuron                            ALS
## SRR1027600                     motor neuron                            ALS
## SRR1027601                     motor neuron                            ALS
## SRR1027602                     motor neuron                            ALS
## SRR1027603                     motor neuron                            ALS
## SRR1027604                     motor neuron                            ALS
## SRR1027605                     motor neuron                            ALS
## SRR1027606                     motor neuron                            ALS
##            Sample.Characteristic.organism.
## SRR1027591                    Homo sapiens
## SRR1027592                    Homo sapiens
## SRR1027593                    Homo sapiens
## SRR1027594                    Homo sapiens
## SRR1027595                    Homo sapiens
## SRR1027596                    Homo sapiens
## SRR1027597                    Homo sapiens
## SRR1027598                    Homo sapiens
## SRR1027599                    Homo sapiens
## SRR1027600                    Homo sapiens
## SRR1027601                    Homo sapiens
## SRR1027602                    Homo sapiens
## SRR1027603                    Homo sapiens
## SRR1027604                    Homo sapiens
## SRR1027605                    Homo sapiens
## SRR1027606                    Homo sapiens
##            Sample.Characteristic.progenitor.cell.type. Factor.Value.disease.
## SRR1027591               induced pluripotent stem cell                normal
## SRR1027592               induced pluripotent stem cell                normal
## SRR1027593               induced pluripotent stem cell                normal
## SRR1027594               induced pluripotent stem cell                normal
## SRR1027595               induced pluripotent stem cell                normal
## SRR1027596               induced pluripotent stem cell                normal
## SRR1027597               induced pluripotent stem cell                normal
## SRR1027598               induced pluripotent stem cell                normal
## SRR1027599               induced pluripotent stem cell                   ALS
## SRR1027600               induced pluripotent stem cell                   ALS
## SRR1027601               induced pluripotent stem cell                   ALS
## SRR1027602               induced pluripotent stem cell                   ALS
## SRR1027603               induced pluripotent stem cell                   ALS
## SRR1027604               induced pluripotent stem cell                   ALS
## SRR1027605               induced pluripotent stem cell                   ALS
## SRR1027606               induced pluripotent stem cell                   ALS
##            Analysed
## SRR1027591      Yes
## SRR1027592      Yes
## SRR1027593      Yes
## SRR1027594      Yes
## SRR1027595      Yes
## SRR1027596      Yes
## SRR1027597      Yes
## SRR1027598      Yes
## SRR1027599      Yes
## SRR1027600      Yes
## SRR1027601      Yes
## SRR1027602      Yes
## SRR1027603      Yes
## SRR1027604      Yes
## SRR1027605      Yes
## SRR1027606      Yes
colnames(colData)
## [1] "Sample.Characteristic.cell.type."           
## [2] "Sample.Characteristic.disease."             
## [3] "Sample.Characteristic.organism."            
## [4] "Sample.Characteristic.progenitor.cell.type."
## [5] "Factor.Value.disease."                      
## [6] "Analysed"
rownames(colData)
##  [1] "SRR1027591" "SRR1027592" "SRR1027593" "SRR1027594" "SRR1027595"
##  [6] "SRR1027596" "SRR1027597" "SRR1027598" "SRR1027599" "SRR1027600"
## [11] "SRR1027601" "SRR1027602" "SRR1027603" "SRR1027604" "SRR1027605"
## [16] "SRR1027606"
```
```
dim(colData)
## [1] 16  6
```
```
str(colData)
## 'data.frame':    16 obs. of  6 variables:
##  $ Sample.Characteristic.cell.type.           : Factor w/ 1 level "motor neuron": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Sample.Characteristic.disease.             : Factor w/ 2 levels "ALS","normal": 2 2 2 2 2 2 2 2 1 1 ...
##  $ Sample.Characteristic.organism.            : Factor w/ 1 level "Homo sapiens": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Sample.Characteristic.progenitor.cell.type.: Factor w/ 1 level "induced pluripotent stem cell": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Factor.Value.disease.                      : Factor w/ 2 levels "ALS","normal": 2 2 2 2 2 2 2 2 1 1 ...
##  $ Analysed                                   : Factor w/ 1 level "Yes": 1 1 1 1 1 1 1 1 1 1 ...
```

## A design formula:
The design formula requires categorical variables to be factors. With a factor, DESeq2 knows how to build comparisons (baseline vs other levels).
```
#define the design formula
designFormula <- "~Sample.Characteristic.disease."
```
```
# making sure the row names in colData matches to column names in countData
all(colnames(countData) %in% rownames(colData))
```
```
## [1] TRUE
```
```
# Reorder colData to match the column names of countData
colData<- colData[colnames(countData),]
all(colnames(countData) == rownames(colData))
## [1] TRUE
```

## construct a DESeqDataSet object:
```
#create a DESeq dataset object from the countData and the colData 
dds <- DESeqDataSetFromMatrix(countData = countData,   colData = colData,  design = as.formula(designFormula))

print(dds)
```
```
## class: DESeqDataSet 
## dim: 58735 16 
## metadata(1): version
## assays(1): counts
## rownames(58735): ENSG00000000003 ENSG00000000005 ... ENSG00000285993
##   ENSG00000285994
## rowData names(0):
## colnames(16): SRR1027602 SRR1027592 ... SRR1027603 SRR1027606
## colData names(6): Sample.Characteristic.cell.type.
##   Sample.Characteristic.disease. ... Factor.Value.disease. Analysed
#dim: 58735 16 
#assays(1) This means my object has 1 assay raw counts.
```

## Filtering non-expressed genes:
we calculate the sum across a particular row it will give the number of samples that gene is expressed in.
```
is_expressed <- assay(dds) >= 1
head(is_expressed) #logical vector
```
```
##                 SRR1027602 SRR1027592 SRR1027598 SRR1027591 SRR1027600
## ENSG00000000003       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000005       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000419       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000457       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000460       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000938      FALSE      FALSE      FALSE      FALSE      FALSE
##                 SRR1027596 SRR1027601 SRR1027599 SRR1027605 SRR1027604
## ENSG00000000003       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000005       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000419       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000457       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000460       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000938      FALSE      FALSE      FALSE      FALSE      FALSE
##                 SRR1027595 SRR1027594 SRR1027593 SRR1027597 SRR1027603
## ENSG00000000003       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000005       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000419       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000457       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000460       TRUE       TRUE       TRUE       TRUE       TRUE
## ENSG00000000938      FALSE       TRUE       TRUE      FALSE      FALSE
##                 SRR1027606
## ENSG00000000003       TRUE
## ENSG00000000005       TRUE
## ENSG00000000419       TRUE
## ENSG00000000457       TRUE
## ENSG00000000460       TRUE
## ENSG00000000938      FALSE
sampleCount<-apply(is_expressed,1,function(x) sum(x) )
head(sampleCount,10)
## ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
##              16              16              16              16              16 
## ENSG00000000938 ENSG00000000971 ENSG00000001036 ENSG00000001084 ENSG00000001167 
##               2              15              16              16              16
```
filtering removing rows with low gene counts and keeping rows that have at least 1 read total
```
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
dds
```
```
## class: DESeqDataSet 
## dim: 36184 16 
## metadata(1): version
## assays(1): counts
## rownames(36184): ENSG00000000003 ENSG00000000005 ... ENSG00000285991
##   ENSG00000285994
## rowData names(0):
## colnames(16): SRR1027602 SRR1027592 ... SRR1027603 SRR1027606
## colData names(6): Sample.Characteristic.cell.type.
##   Sample.Characteristic.disease. ... Factor.Value.disease. Analysed
nrow(dds)
## [1] 36184
table(keep)
## keep
## FALSE  TRUE 
## 22551 36184
```
## Visualizing count distributions:
- plot a histogram
```
hist(sampleCount, main="Number of samples a gene is expressed in",xlab="Sample Count")
```
<img width="870" height="600" alt="image" src="https://github.com/user-attachments/assets/04ffcccb-7d16-4dff-8eb2-1a939ec8b541" />

```
dev.copy(pdf, "images/count distribution/histogram.pdf", width = 7, height = 5)
## pdf 
##   3
dev.off()
## png 
##   2
```

- box plot
```
boxplot(log10(assay(dds)+1))
```
<img width="810" height="495" alt="image" src="https://github.com/user-attachments/assets/1ca9ea54-bff4-4334-9add-a255324bba11" />

```
dev.copy(pdf, "images/count distribution/boxplot.pdf", width = 7, height = 5)
## pdf 
##   3
dev.off()
## png 
##   2
```

## DESeq2 pipeline

### set the factor level:
```
dds$Sample.Characteristic.disease. <- relevel(dds$Sample.Characteristic.disease., ref = "normal")
DESeq()
dds <- DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## -- replacing outliers and refitting for 133 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
## estimating dispersions
## fitting model and testing
```
## DE Results
```
DEresults <- results(dds, contrast = c("Sample.Characteristic.disease.", 'ALS', 'normal'))
DEresults
## log2 fold change (MLE): Sample.Characteristic.disease. ALS vs normal 
## Wald test p-value: Sample.Characteristic.disease. ALS vs normal 
## DataFrame with 36184 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE       stat      pvalue
##                 <numeric>      <numeric> <numeric>  <numeric>   <numeric>
## ENSG00000000003  973.2153      -0.179020  0.170229  -1.051640 0.292964562
## ENSG00000000005   18.3155      -3.278939  0.865188  -3.789858 0.000150734
## ENSG00000000419  251.5193       0.036135  0.113309   0.318906 0.749798058
## ENSG00000000457  167.2970       0.139029  0.207523   0.669945 0.502892559
## ENSG00000000460   82.4290      -0.333785  0.231292  -1.443132 0.148983400
## ...                   ...            ...       ...        ...         ...
## ENSG00000285982  1.693894      0.0456475  1.061789  0.0429912    0.965709
## ENSG00000285985  0.113852     -0.0459999  3.060383 -0.0150308    0.988008
## ENSG00000285987  0.130454     -0.7673390  3.060383 -0.2507330    0.802021
## ENSG00000285991  2.701249      0.7924319  0.644693  1.2291613    0.219011
## ENSG00000285994  0.730520      0.3319741  1.165535  0.2848254    0.775778
##                       padj
##                  <numeric>
## ENSG00000000003 0.69196438
## ENSG00000000005 0.00677323
## ENSG00000000419 0.93237588
## ENSG00000000457 0.83443415
## ENSG00000000460 0.51440606
## ...                    ...
## ENSG00000285982         NA
## ENSG00000285985         NA
## ENSG00000285987         NA
## ENSG00000285991   0.617014
## ENSG00000285994         NA
summary(DEresults)
```
```
## out of 36184 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 594, 1.6%
## LFC < 0 (down)     : 1237, 3.4%
## outliers [1]       : 0, 0%
## low counts [2]     : 14031, 39%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```
with the treshold of pdaj<0.05 and log2FC of 0

```
sig <- subset(DEresults, padj < 0.05 & abs(log2FoldChange) > 0)
sig
```
```
## log2 fold change (MLE): Sample.Characteristic.disease. ALS vs normal 
## Wald test p-value: Sample.Characteristic.disease. ALS vs normal 
## DataFrame with 1242 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE      stat      pvalue
##                 <numeric>      <numeric> <numeric> <numeric>   <numeric>
## ENSG00000000005  18.31550      -3.278939  0.865188  -3.78986 1.50734e-04
## ENSG00000002746 287.81600       0.561139  0.103276   5.43339 5.52935e-08
## ENSG00000002933 111.76963      -2.203240  0.689736  -3.19432 1.40159e-03
## ENSG00000003402  98.84312      -1.119168  0.221858  -5.04451 4.54678e-07
## ENSG00000005102   7.93967      -3.119792  0.925492  -3.37096 7.49077e-04
## ...                   ...            ...       ...       ...         ...
## ENSG00000280255  26.54421       1.556753  0.504789   3.08397 0.002042578
## ENSG00000280351  33.30459       0.588089  0.171367   3.43174 0.000599717
## ENSG00000281406  51.00396       0.670018  0.193594   3.46094 0.000538288
## ENSG00000281880  41.24604       0.617419  0.185146   3.33478 0.000853681
## ENSG00000284681   3.26689      -2.365138  0.646353  -3.65921 0.000252997
##                        padj
##                   <numeric>
## ENSG00000000005 6.77323e-03
## ENSG00000002746 1.91393e-05
## ENSG00000002933 3.19110e-02
## ENSG00000003402 9.01557e-05
## ENSG00000005102 2.07688e-02
## ...                     ...
## ENSG00000280255  0.04069175
## ENSG00000280351  0.01812002
## ENSG00000281406  0.01708409
## ENSG00000281880  0.02267578
## ENSG00000284681  0.00995498
```
summary(sig)
```
## 
## out of 1242 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 388, 31%
## LFC < 0 (down)     : 854, 69%
## outliers [1]       : 0, 0%
## low counts [2]     : 0, 0%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```
## Diagnostic plots
### MAplot
- X-axis mean of normalized counts across all samples (in log10 scale).
- Y-axis estimated log2 fold change between conditions

blue points highlight **significant DE genes.**

Helps us quickly see whether we have many up-/down-regulated genes Most points are expected to be on the horizontal 0 line (most genes are not expected to be differentially expressed).
```
plotMA(DEresults,ylim = c(-5, 5))
```
<img width="867" height="580" alt="image" src="https://github.com/user-attachments/assets/700f30aa-08cf-4703-9d0b-cd98bad39c11" />

```

dev.copy(pdf, "images/Diagnostic plots/MA plot.pdf", width = 7, height = 5)
## pdf 
##   3
dev.off()
## png 
##   2
```

### P-value distribution

Peak at low p-values (< 0.05) this indicates that some genes are truly differentially expressed, they show strong evidence of change between conditions.

Uniform distribution at high p-values (> 0.1) this suggests that the rest of the genes are not significantly different and behave randomly, as expected under the null hypothesis.

after this analysis the padj values are meaningful
```
library(ggplot2)
ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) + 
  geom_histogram(colour='black',fill='white', bins = 100)
```
<img width="906" height="621" alt="image" src="https://github.com/user-attachments/assets/f31bb33b-957b-47e6-9322-edc6e7ddf775" />

```
dev.copy(pdf, "images/Diagnostic plots/p-value distribution.pdf", width = 7, height = 5)
## pdf 
##   3
dev.off()
## png 
##   2
```

### Principal Component Analysis(PCA) plot
we userlog()to transformed our counts from dds. It stabilizes the variance across the range of counts, making the data more suitable for PCA and visualization. PCA plot helps us to see if our samples is cluster by group ntop = 500 means the plot uses the 500 most variable genes (those that differ the most between samples).
```
rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'Sample.Characteristic.disease.') + 
   theme_bw()
```
<img width="855" height="485" alt="image" src="https://github.com/user-attachments/assets/765697ed-eb15-4a4e-a0df-0101166c0b26" />

```
dev.copy(pdf, "images/Diagnostic plots/PCA plot.pdf", width = 7, height = 5)
## pdf 
##   3
dev.off()
## png 
##   2
```

## Heatmap of the sample-to-sample distances
we useassay()to extract our rld object from dds which is log of our normalized count in the count matrix ,then we transpose the matrix with t()it convert samples into rows so we can apply dist() function that compute the distance between the samples based on the gene expression value which reflect which samples are similar based on expression

```
sampleDists <- dist(t(assay(rld)))
sampleDists # Calculate distances between samples
##            SRR1027602 SRR1027592 SRR1027598 SRR1027591 SRR1027600 SRR1027596
## SRR1027592   61.73333                                                       
## SRR1027598   91.89799   84.67398                                            
## SRR1027591   61.76745   27.60546   84.76367                                 
## SRR1027600   46.72881   62.88743   94.77368   62.99777                      
## SRR1027596   56.92073   71.46845   83.16821   71.43202   55.48370           
## SRR1027601   28.33582   61.77847   92.06042   61.79718   46.89798   56.89357
## SRR1027599   46.57142   62.80072   94.81726   62.71245   28.75762   55.42064
## SRR1027605   47.06897   67.61397   93.61219   67.75857   46.63264   57.41148
## SRR1027604   81.57682   91.47860  114.81154   91.40393   80.90335   85.75781
## SRR1027595   56.76214   71.16970   82.85766   71.10791   55.28436   27.84424
## SRR1027594   87.73167   99.19695  117.07802   99.29571   93.87009   96.51548
## SRR1027593   87.53123   99.27191  117.13823   99.21076   93.77802   96.44761
## SRR1027597   91.55230   84.57789   28.45763   84.56347   94.66050   83.06521
## SRR1027603   81.55297   91.29552  114.72211   91.30034   80.76218   85.60008
## SRR1027606   47.13116   67.54769   93.56547   67.73480   46.67217   57.60802
##            SRR1027601 SRR1027599 SRR1027605 SRR1027604 SRR1027595 SRR1027594
## SRR1027592                                                                  
## SRR1027598                                                                  
## SRR1027591                                                                  
## SRR1027600                                                                  
## SRR1027596                                                                  
## SRR1027601                                                                  
## SRR1027599   46.56235                                                       
## SRR1027605   46.91697   46.60942                                            
## SRR1027604   81.57478   80.86706   82.59335                                 
## SRR1027595   56.73077   55.28678   57.40060   85.63377                      
## SRR1027594   87.82767   93.74380   87.21553   68.73593   96.57338           
## SRR1027593   87.66951   93.57336   87.08501   68.43427   96.43501   30.03463
## SRR1027597   91.83993   94.66865   93.34724  114.48042   82.68683  116.97809
## SRR1027603   81.50189   80.72417   82.52343   28.08934   85.47890   68.92698
## SRR1027606   46.96302   46.67026   26.02910   82.47418   57.53511   86.98132
##            SRR1027593 SRR1027597 SRR1027603
## SRR1027592                                 
## SRR1027598                                 
## SRR1027591                                 
## SRR1027600                                 
## SRR1027596                                 
## SRR1027601                                 
## SRR1027599                                 
## SRR1027605                                 
## SRR1027604                                 
## SRR1027595                                 
## SRR1027594                                 
## SRR1027593                                 
## SRR1027597  116.92219                      
## SRR1027603   68.74422  114.39830           
## SRR1027606   86.83102   93.39868   82.45092
```
```
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
```
```
#colData(dds)
annotation_col <- as.data.frame(
  colData(dds)[, 'Sample.Characteristic.disease.', drop = FALSE] #Add sample annotations
)
print(annotation_col)
```
```
##            Sample.Characteristic.disease.
## SRR1027602                            ALS
## SRR1027592                         normal
## SRR1027598                         normal
## SRR1027591                         normal
## SRR1027600                            ALS
## SRR1027596                         normal
## SRR1027601                            ALS
## SRR1027599                            ALS
## SRR1027605                            ALS
## SRR1027604                            ALS
## SRR1027595                         normal
## SRR1027594                         normal
## SRR1027593                         normal
## SRR1027597                         normal
## SRR1027603                            ALS
## SRR1027606                            ALS
```
```
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
      annotation_col = annotation_col,  
  main = "Sample-to-sample distances" )
```
<img width="881" height="620" alt="image" src="https://github.com/user-attachments/assets/16dde1ef-1943-4e70-9e08-ae0473602627" />

```
dev.copy(pdf, "images/Diagnostic plots/Heatmap of the sample-to-sample distances.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

## Volcano plot
```
res <- as.data.frame(DEresults)

# Create a column for significance
res$threshold <- "Not significant"
res$threshold[res$padj < 0.1 & res$log2FoldChange > 0] <- "up"
res$threshold[res$padj < 0.1 & res$log2FoldChange < 0] <- "down"

# Volcano plot
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "Not significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +     # FC cutoffs
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + # padj cutoff
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_minimal()
```
<img width="916" height="617" alt="image" src="https://github.com/user-attachments/assets/2ace4f70-b094-4c2f-82cd-299103cb54c3" />

```
dev.copy(pdf, "images/Diagnostic plots/Volcano plot.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```
## Heatmap
Heatmap of the top 50 most differentially expressed genes

```
vsd <- vst(dds, blind = TRUE) ## Variance stabilizing transform 

top50_vsd<-head(order(rowVars(assay(vsd)), decreasing = TRUE), 50) #: Select the top 50 variable genes

#Extract and scale the matrices
mat_vsd <- assay(vsd)[top50_vsd, ]
mat_vsd <- t(scale(t(mat_vsd)))


annotation_col <- as.data.frame(colData(dds)[, 'Sample.Characteristic.disease.', drop = FALSE])

pheatmap(mat_vsd,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,)
```

<img width="842" height="632" alt="image" src="https://github.com/user-attachments/assets/7e68eaff-62f7-404d-9f97-e4ad0a576926" />

```
dev.copy(pdf, "images/Diagnostic plots/Heatmap.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```
## Select DE
1838 genes are expressed at the treshold of 0.1 and log2FC of 0

```
#remove genes with NA values 
DEG<- DEresults[!is.na(DEresults$padj),]
#select genes with adjusted p-values below 0.1
DEG <- DEG[DEG$padj < 0.1,]
#select genes with absolute log2 fold change above 1 (two-fold change)
DEG <- DEG[abs(DEG$log2FoldChange) > 0,]
DEG
```

```
## log2 fold change (MLE): Sample.Characteristic.disease. ALS vs normal 
## Wald test p-value: Sample.Characteristic.disease. ALS vs normal 
## DataFrame with 1831 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE      stat      pvalue
##                 <numeric>      <numeric> <numeric> <numeric>   <numeric>
## ENSG00000000005  18.31550      -3.278939  0.865188  -3.78986 1.50734e-04
## ENSG00000001036 384.06620      -0.385792  0.131641  -2.93063 3.38270e-03
## ENSG00000002745   5.68631      -1.821268  0.648892  -2.80674 5.00464e-03
## ENSG00000002746 287.81600       0.561139  0.103276   5.43339 5.52935e-08
## ENSG00000002933 111.76963      -2.203240  0.689736  -3.19432 1.40159e-03
## ...                   ...            ...       ...       ...         ...
## ENSG00000284681   3.26689      -2.365138  0.646353  -3.65921 0.000252997
## ENSG00000284770  63.58416      -0.439279  0.148071  -2.96669 0.003010284
## ENSG00000285331  17.56958       0.674047  0.251333   2.68189 0.007320795
## ENSG00000285410 401.44092       0.367307  0.127400   2.88311 0.003937755
## ENSG00000285525   3.80209       1.889286  0.700857   2.69568 0.007024540
##                        padj
##                   <numeric>
## ENSG00000000005 6.77323e-03
## ENSG00000001036 5.69431e-02
## ENSG00000002745 7.34223e-02
## ENSG00000002746 1.91393e-05
## ENSG00000002933 3.19110e-02
## ...                     ...
## ENSG00000284681  0.00995498
## ENSG00000284770  0.05259213
## ENSG00000285331  0.09369011
## ENSG00000285410  0.06251938
## ENSG00000285525  0.09151428
```

summary(DEG)

```
## 
## out of 1831 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 594, 32%
## LFC < 0 (down)     : 1237, 68%
## outliers [1]       : 0, 0%
## low counts [2]     : 0, 0%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

1242 genes are expressed at the treshold of padj<0.05 and logFC2 of 0

```
#remove genes with NA values 
DE<- DEresults[!is.na(DEresults$padj),]
#select genes with adjusted p-values below 0.05
DE <- DE[DE$padj < 0.05,]

DE <- DE[abs(DE$log2FoldChange) > 0,]
DE
```

```
## log2 fold change (MLE): Sample.Characteristic.disease. ALS vs normal 
## Wald test p-value: Sample.Characteristic.disease. ALS vs normal 
## DataFrame with 1242 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE      stat      pvalue
##                 <numeric>      <numeric> <numeric> <numeric>   <numeric>
## ENSG00000000005  18.31550      -3.278939  0.865188  -3.78986 1.50734e-04
## ENSG00000002746 287.81600       0.561139  0.103276   5.43339 5.52935e-08
## ENSG00000002933 111.76963      -2.203240  0.689736  -3.19432 1.40159e-03
## ENSG00000003402  98.84312      -1.119168  0.221858  -5.04451 4.54678e-07
## ENSG00000005102   7.93967      -3.119792  0.925492  -3.37096 7.49077e-04
## ...                   ...            ...       ...       ...         ...
## ENSG00000280255  26.54421       1.556753  0.504789   3.08397 0.002042578
## ENSG00000280351  33.30459       0.588089  0.171367   3.43174 0.000599717
## ENSG00000281406  51.00396       0.670018  0.193594   3.46094 0.000538288
## ENSG00000281880  41.24604       0.617419  0.185146   3.33478 0.000853681
## ENSG00000284681   3.26689      -2.365138  0.646353  -3.65921 0.000252997
##                        padj
##                   <numeric>
## ENSG00000000005 6.77323e-03
## ENSG00000002746 1.91393e-05
## ENSG00000002933 3.19110e-02
## ENSG00000003402 9.01557e-05
## ENSG00000005102 2.07688e-02
## ...                     ...
## ENSG00000280255  0.04069175
## ENSG00000280351  0.01812002
## ENSG00000281406  0.01708409
## ENSG00000281880  0.02267578
## ENSG00000284681  0.00995498
```

summary(DE)

```
## 
## out of 1242 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 388, 31%
## LFC < 0 (down)     : 854, 69%
## outliers [1]       : 0, 0%
## low counts [2]     : 0, 0%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
write.table(DE, file = "filtered_DEGs.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

## Add annotation for genes
Here we are adding a new column with official human gene symbol which crosspondes to the Ensemble ID(rowname) in our DE. mapIds()= ID translator for genes. It takes our input IDs (keys), knows what kind they are (keytype), and returns another ID type (column) using the chosen organism database.

```
library("AnnotationDbi")
library("org.Hs.eg.db")

DE$symbol <- mapIds(
  org.Hs.eg.db,            # human annotation database
  keys = rownames(DE),     # Ensembl IDs
  column = "SYMBOL",       # return gene symbols
  keytype = "ENSEMBL",     # input IDs are Ensembl
  multiVals = "first"      # take the first match if multiple
)
head(DE$symbol,10)
```

```
## ENSG00000000005 ENSG00000002746 ENSG00000002933 ENSG00000003402 ENSG00000005102 
##          "TNMD"         "HECW1"      "TMEM176A"         "CFLAR"         "MEOX1" 
## ENSG00000005243 ENSG00000005436 ENSG00000005471 ENSG00000005893 ENSG00000007384 
##         "COPZ2"         "GCFC2"         "ABCB4"         "LAMP2"        "RHBDF1"
```

here we are mapping Ensembl IDs to Entrez IDs,a new column ENTREZID gets added to our DE .Entrez IDs are the preferred input for many enrichment tools like GO, KEGG, Reactome.

```
DE$entrez <- mapIds(org.Hs.eg.db,
                    keys=rownames(DE),
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
head(DE$entrez, 10)
```

```
## ENSG00000000005 ENSG00000002746 ENSG00000002933 ENSG00000003402 ENSG00000005102 
##         "64102"         "23072"         "55365"          "8837"          "4222" 
## ENSG00000005243 ENSG00000005436 ENSG00000005471 ENSG00000005893 ENSG00000007384 
##         "51226"          "6936"          "5244"          "3920"         "64285"
```

colnames(DE)

```
## [1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
## [5] "pvalue"         "padj"           "symbol"         "entrez"
```

```
dim(DE)
## [1] 1242    8
```

# GO enrichment analysis
enrichGO()tests whether our gene list of DE is enriched in GO terms compared to the whole genome background to get our list of over represented pathways and this type of analysis is called over-representation analysis (ORA)

ont specifies which branch of GO to use:

- BP = Biological Process

- MF = Molecular Function

- CC = Cellular Component
```
suppressMessages(suppressWarnings(library(clusterProfiler)))

GO_BP <- enrichGO(DE$symbol, OrgDb = "org.Hs.eg.db", 
                  keyType = "SYMBOL", ont = "BP")
#head(GO_BP)
dim(GO_BP)
## [1] 672   9
```
```
GO_MF <- enrichGO(DE$symbol, OrgDb = "org.Hs.eg.db", 
                  keyType = "SYMBOL", ont = "MF")
dim(GO_MF)
## [1] 37  9
```
```
GO_CC <- enrichGO(DE$symbol, OrgDb = "org.Hs.eg.db", 
                  keyType = "SYMBOL", ont = "CC")
dim(GO_CC)
## [1] 74  9
```
### Plotting the results
The dotplot is one of the most popular ways to visualize enrichment results This makes it easy to spot the strongest enrichment.
```
dotplot(GO_BP, title = "Biological Process")
```
<img width="572" height="346" alt="image" src="https://github.com/user-attachments/assets/64a1294c-6436-427b-887e-cc1445500252" />

```
dev.copy(pdf, "images/Dot Plots/dot plot of biological process.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

```
dotplot(GO_MF, title = "Molecular Functions")
```
<img width="480" height="342" alt="image" src="https://github.com/user-attachments/assets/027fae30-c59f-4fac-b2cc-8129b56650a1" />

```
dev.copy(pdf, "images/Dot Plots/dot plot of molecular function.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```
```
dotplot(GO_CC, title = "Cellular Component")
```
<img width="480" height="335" alt="image" src="https://github.com/user-attachments/assets/37374eb8-6ed4-40e3-a8cd-fa611680daa5" />
```
dev.copy(pdf, "images/Dot Plots/dot plot of cellular component.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```
## Analyze UP and DOWN regulated genes

these two lines are filtering our differential expression (DE) results into upregulated and downregulated gene
```
DE_up <- DE[DE$log2FoldChange > 0,]
DE_up#significantly upregulated genes.
```
```
## log2 fold change (MLE): Sample.Characteristic.disease. ALS vs normal 
## Wald test p-value: Sample.Characteristic.disease. ALS vs normal 
## DataFrame with 388 rows and 8 columns
##                  baseMean log2FoldChange     lfcSE      stat      pvalue
##                 <numeric>      <numeric> <numeric> <numeric>   <numeric>
## ENSG00000002746   287.816       0.561139 0.1032760   5.43339 5.52935e-08
## ENSG00000008277   701.882       0.258166 0.0771115   3.34796 8.14101e-04
## ENSG00000013523   467.195       0.490027 0.1507970   3.24958 1.15575e-03
## ENSG00000019995   547.966       0.268969 0.0786590   3.41943 6.27515e-04
## ENSG00000020129   889.834       0.395754 0.1223306   3.23512 1.21594e-03
## ...                   ...            ...       ...       ...         ...
## ENSG00000280135   11.1616       1.142246  0.374192   3.05257 0.002268944
## ENSG00000280255   26.5442       1.556753  0.504789   3.08397 0.002042578
## ENSG00000280351   33.3046       0.588089  0.171367   3.43174 0.000599717
## ENSG00000281406   51.0040       0.670018  0.193594   3.46094 0.000538288
## ENSG00000281880   41.2460       0.617419  0.185146   3.33478 0.000853681
##                        padj      symbol      entrez
##                   <numeric> <character> <character>
## ENSG00000002746 1.91393e-05       HECW1       23072
## ENSG00000008277 2.19135e-02      ADAM22       53616
## ENSG00000013523 2.83516e-02      ANGEL1       23357
## ENSG00000019995 1.86595e-02      ZRANB1       54764
## ENSG00000020129 2.93983e-02        NCDN       23154
## ...                     ...         ...         ...
## ENSG00000280135   0.0436538          NA          NA
## ENSG00000280255   0.0406917          NA          NA
## ENSG00000280351   0.0181200          NA          NA
## ENSG00000281406   0.0170841     BLACAT1   101669762
## ENSG00000281880   0.0226758    PAX6-AS1      440034
```
```
DE_down <- DE[DE$log2FoldChange < 0,]#significantly downregulated genes.
DE_down
```
```
## log2 fold change (MLE): Sample.Characteristic.disease. ALS vs normal 
## Wald test p-value: Sample.Characteristic.disease. ALS vs normal 
## DataFrame with 854 rows and 8 columns
##                  baseMean log2FoldChange     lfcSE      stat      pvalue
##                 <numeric>      <numeric> <numeric> <numeric>   <numeric>
## ENSG00000000005  18.31550       -3.27894  0.865188  -3.78986 1.50734e-04
## ENSG00000002933 111.76963       -2.20324  0.689736  -3.19432 1.40159e-03
## ENSG00000003402  98.84312       -1.11917  0.221858  -5.04451 4.54678e-07
## ENSG00000005102   7.93967       -3.11979  0.925492  -3.37096 7.49077e-04
## ENSG00000005243  49.75207       -1.36503  0.412920  -3.30579 9.47102e-04
## ...                   ...            ...       ...       ...         ...
## ENSG00000278910   9.02932      -3.878523  1.121934  -3.45700 0.000546228
## ENSG00000279233  21.72931      -2.141276  0.688603  -3.10959 0.001873456
## ENSG00000279821  19.18041      -0.843605  0.262763  -3.21052 0.001324964
## ENSG00000280219   7.11920      -2.779411  0.764305  -3.63652 0.000276345
## ENSG00000284681   3.26689      -2.365138  0.646353  -3.65921 0.000252997
##                        padj      symbol      entrez
##                   <numeric> <character> <character>
## ENSG00000000005 6.77323e-03        TNMD       64102
## ENSG00000002933 3.19110e-02    TMEM176A       55365
## ENSG00000003402 9.01557e-05       CFLAR        8837
## ENSG00000005102 2.07688e-02       MEOX1        4222
## ENSG00000005243 2.43401e-02       COPZ2       51226
## ...                     ...         ...         ...
## ENSG00000278910  0.01720604       BANCR   100885775
## ENSG00000279233  0.03851389          NA          NA
## ENSG00000279821  0.03086428          NA          NA
## ENSG00000280219  0.01060984          NA          NA
## ENSG00000284681  0.00995498          NA          NA
```
# Enrichment analysis
It contains enrichment results for both UP and DOWN lists, tagged with a group label.

## analysis with BP
This code compares GO Biological Process enrichment between upregulated and downregulated genes, helping us understand if different biological processes are driving up vs down changes.

number of enriched terms found for each gene cluster: - UP: 0 (No GO terms were significantly enriched for the UP genes.) - DOWN: 20 (20 GO terms were significantly enriched for the DOWN genes)
```
GO_up_down <- compareCluster(list(UP = DE_up$symbol, DOWN = DE_down$symbol),  # two gene lists 
                             fun = "enrichGO",  # enrichment function
                             OrgDb = "org.Mm.eg.db", # organism annotation database
                             keyType = "SYMBOL",  # input gene IDs
                             ont = "BP", # GO ontology 
      pvalueCutoff = 0.2)
                         
GO_up_down
```
```
## #
## # Result of Comparing 2 gene clusters 
## #
## #.. @fun      enrichGO 
## #.. @geneClusters    List of 2
##  $ UP  : Named chr [1:388] "HECW1" "ADAM22" "ANGEL1" "ZRANB1" ...
##   ..- attr(*, "names")= chr [1:388] "ENSG00000002746" "ENSG00000008277" "ENSG00000013523" "ENSG00000019995" ...
##  $ DOWN: Named chr [1:854] "TNMD" "TMEM176A" "CFLAR" "MEOX1" ...
##   ..- attr(*, "names")= chr [1:854] "ENSG00000000005" "ENSG00000002933" "ENSG00000003402" "ENSG00000005102" ...
## #...Result   'data.frame':   20 obs. of  10 variables:
##  $ Cluster    : Factor w/ 2 levels "UP","DOWN": 2 2 2 2 2 2 2 2 2 2 ...
##  $ ID         : chr  "GO:0040015" "GO:0044030" "GO:0071514" "GO:0006883" ...
##  $ Description: chr  "negative regulation of multicellular organism growth" "regulation of DNA methylation" "genomic imprinting" "intracellular sodium ion homeostasis" ...
##  $ GeneRatio  : chr  "1/2" "1/2" "1/2" "1/2" ...
##  $ BgRatio    : chr  "17/28891" "19/28891" "20/28891" "27/28891" ...
##  $ pvalue     : num  0.00118 0.00131 0.00138 0.00187 0.00291 ...
##  $ p.adjust   : num  0.00923 0.00923 0.00923 0.00934 0.01008 ...
##  $ qvalue     : logi  NA NA NA NA NA NA ...
##  $ geneID     : chr  "H19" "H19" "H19" "C7" ...
##  $ Count      : int  1 1 1 1 1 1 1 1 1 1 ...
## #.. number of enriched terms found for each gene cluster:
## #..   UP: 0 
## #..   DOWN: 20 
## #
## #...Citation
## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, 
## W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. 
## clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. 
## The Innovation. 2021, 2(3):100141
```

```
head(as.data.frame(GO_up_down))
```

```
##   Cluster         ID                                          Description
## 1    DOWN GO:0040015 negative regulation of multicellular organism growth
## 2    DOWN GO:0044030                        regulation of DNA methylation
## 3    DOWN GO:0071514                                   genomic imprinting
## 4    DOWN GO:0006883                 intracellular sodium ion homeostasis
## 5    DOWN GO:0055078                               sodium ion homeostasis
## 6    DOWN GO:0006305                                       DNA alkylation
##   GeneRatio  BgRatio      pvalue    p.adjust qvalue geneID Count
## 1       1/2 17/28891 0.001176511 0.009227060     NA    H19     1
## 2       1/2 19/28891 0.001314879 0.009227060     NA    H19     1
## 3       1/2 20/28891 0.001384059 0.009227060     NA    H19     1
## 4       1/2 27/28891 0.001868253 0.009341266     NA     C7     1
## 5       1/2 42/28891 0.002905417 0.010078446     NA     C7     1
## 6       1/2 51/28891 0.003527456 0.010078446     NA    H19     1
```

This will plot enriched biological processes for both UP and DOWN genes, so we can immediately see which functions are enriched in each set.

```
dotplot(GO_up_down)
```

<img width="872" height="608" alt="image" src="https://github.com/user-attachments/assets/81022e0f-f77c-440b-b347-53a4ab49b75f" />

```
dev.copy(pdf, "images/Dot Plots/GO_UP_DOWN plot for BP.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

## Analysis with MF
```
GO_Up_Down <- compareCluster(list(UP = DE_up$symbol, DOWN = DE_down$symbol),  # two gene lists 
                             fun = "enrichGO",  # enrichment function
                             OrgDb = "org.Mm.eg.db", # organism annotation database
                             keyType = "SYMBOL",  # input gene IDs
                             ont = "MF", # GO ontology 
      pvalueCutoff = 0.7)
                         
GO_Up_Down
```

```
## NULL
#as.data.frame(GO_up_down)
```

# Pathways databases
## REACTOME
we are doing a comparative Reactome pathway enrichment analysis. REACTOMEupdown is a compareClusterResult object containing enriched Reactome pathways for both UP and DOWN lists. Number of enriched terms found for each gene cluster: - UP: 6 (UP-regulated genes: 6 Reactome pathways enriched) - DOWN: 67 (DOWN-regulated genes: 67 Reactome pathways enriched.)
```
suppressMessages(suppressWarnings(library(ReactomePA)))
REACTOMEupdown <-
  compareCluster(list(UP = DE_up$entrez, DOWN = DE_down$entrez),
    fun = "enrichPathway", # Reactome enrichment function
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE  # convert Entrez IDs back to gene symbols
    )
REACTOMEupdown
```

```
## #
## # Result of Comparing 2 gene clusters 
## #
## #.. @fun      enrichPathway 
## #.. @geneClusters    List of 2
##  $ UP  : Named chr [1:388] "23072" "53616" "23357" "54764" ...
##   ..- attr(*, "names")= chr [1:388] "ENSG00000002746" "ENSG00000008277" "ENSG00000013523" "ENSG00000019995" ...
##  $ DOWN: Named chr [1:854] "64102" "55365" "8837" "4222" ...
##   ..- attr(*, "names")= chr [1:854] "ENSG00000000005" "ENSG00000002933" "ENSG00000003402" "ENSG00000005102" ...
## #...Result   'data.frame':   73 obs. of  10 variables:
##  $ Cluster    : Factor w/ 2 levels "UP","DOWN": 1 1 1 1 1 1 2 2 2 2 ...
##  $ ID         : chr  "R-HSA-112316" "R-HSA-2022923" "R-HSA-6794362" "R-HSA-6794361" ...
##  $ Description: chr  "Neuronal System" "Dermatan sulfate biosynthesis" "Protein-protein interactions at synapses" "Neurexins and neuroligins" ...
##  $ GeneRatio  : chr  "19/156" "4/156" "8/156" "6/156" ...
##  $ BgRatio    : chr  "410/11009" "11/11009" "86/11009" "55/11009" ...
##  $ pvalue     : num  5.22e-06 1.18e-05 2.83e-05 1.20e-04 3.11e-04 ...
##  $ p.adjust   : num  0.00313 0.00356 0.00566 0.01805 0.03742 ...
##  $ qvalue     : num  0.00301 0.00342 0.00545 0.01736 0.03598 ...
##  $ geneID     : chr  "APBA2/KCNK2/SLC32A1/LIN7B/SLC38A1/GAD2/PPFIA2/GRIP1/SLC6A1/GJD2/SHANK2/KCNJ9/GRM5/ADCY5/SLITRK1/GRIN2A/KCNQ3/GNG2/LRRTM3" "CSPG5/NCAN/BCAN/DSEL" "APBA2/LIN7B/PPFIA2/SHANK2/GRM5/SLITRK1/GRIN2A/LRRTM3" "APBA2/LIN7B/SHANK2/GRM5/GRIN2A/LRRTM3" ...
##  $ Count      : int  19 4 8 6 3 3 83 34 41 29 ...
## #.. number of enriched terms found for each gene cluster:
## #..   UP: 6 
## #..   DOWN: 67 
## #
## #...Citation
## 1.  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for
##   reactome pathway analysis and visualization. Molecular BioSystems
##   2016, 12(2):477-479 
## 
## 2.T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, 
## W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. 
## clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. 
## The Innovation. 2021, 2(3):100141
```

```
head(as.data.frame(REACTOMEupdown),10)
```

```
##    Cluster            ID                                 Description GeneRatio
## 1       UP  R-HSA-112316                             Neuronal System    19/156
## 2       UP R-HSA-2022923               Dermatan sulfate biosynthesis     4/156
## 3       UP R-HSA-6794362    Protein-protein interactions at synapses     8/156
## 4       UP R-HSA-6794361                   Neurexins and neuroligins     6/156
## 5       UP R-HSA-9832991     Formation of the posterior neural plate     3/156
## 6       UP R-HSA-9823739      Formation of the anterior neural plate     3/156
## 7     DOWN R-HSA-1474244           Extracellular matrix organization    83/535
## 8     DOWN R-HSA-1474290                          Collagen formation    34/535
## 9     DOWN R-HSA-1474228     Degradation of the extracellular matrix    41/535
## 10    DOWN R-HSA-1650814 Collagen biosynthesis and modifying enzymes    29/535
##      BgRatio       pvalue     p.adjust       qvalue
## 1  410/11009 5.215851e-06 3.134727e-03 3.014213e-03
## 2   11/11009 1.184921e-05 3.560688e-03 3.423798e-03
## 3   86/11009 2.827264e-05 5.663951e-03 5.446202e-03
## 4   55/11009 1.201649e-04 1.805478e-02 1.736067e-02
## 5   10/11009 3.113375e-04 3.742276e-02 3.598406e-02
## 6   11/11009 4.236455e-04 4.243516e-02 4.080375e-02
## 7  300/11009 6.146981e-41 6.226892e-38 5.694046e-38
## 8   90/11009 4.433995e-22 2.245818e-19 2.053640e-19
## 9  140/11009 1.569154e-21 5.298509e-19 4.845106e-19
## 10  67/11009 5.449515e-21 1.380090e-18 1.261993e-18
##                                                                                          geneID
## 1                                                                             APBA2/KCNK2/SLC32A1/LIN7B/SLC38A1/GAD2/PPFIA2/GRIP1/SLC6A1/GJD2/SHANK2/KCNJ9/GRM5/ADCY5/SLITRK1/GRIN2A/KCNQ3/GNG2/LRRTM3
## 2                                                                                          CSPG5/NCAN/BCAN/DSEL                                                      
## 3                                                                                         APBA2/LIN7B/PPFIA2/SHANK2/GRM5/SLITRK1/GRIN2A/LRRTM3
## 4   APBA2/LIN7B/SHANK2/GRM5/GRIN2A/LRRTM3                      
## 5  OTX2/GBX2/POU3                                               
## 6  BX2/POU3F1
## 7  DCN/VCAN/TNC/LTBP1/ELN/LAMC3/COL11A1/ACTN1/P4HA2/CAPN6/FBLN1/CEACAM1/COL16A1/MMP2/NID2/LTBP4/ICAM1/LAMB1/COL9A3/MMP11/TIMP1/PCOLCE/SERPINE1/COL1A1/P3H3/COL12A1/COL9A1/LAMA4/SPARC/FN1/NID1/MFAP2/P3H1/LTBP2/TNN/SDC4/BMP4/PXDN/COL5A1/LOXL2/COL4A2/LAMC1/THBS1/ITGA11/LOXL4/ADAMTS14/FBN2/COL2A1/LUM/COL6A2/HSPG2/ITGA9/COL8A1/SERPINH1/ADAMTS1/ADAMTS3/MMP14/ITGA5/DDR2/CTSS/COL6A3/ADAMTS9/COL1A2/FBN1/MFAP4/LTBP3/BMP1/COL3A1/LAMB2/BSG/EFEMP2/A2M/CD151/BGN/COL4A1/COL14A1/COL4A5/LAMA2/COL13A1/COL5A2/COL15A1/COL6A6/ITGA1
## 8  COL11A1/P4HA2/COL16A1/COL9A3/PCOLCE/COL1A1/P3H3/COL12A1/COL9A1/P3H1/PXDN/COL5A1/LOXL2/COL4A2/LOXL4/ADAMTS14/COL2A1/COL6A2/COL8A1/SERPINH1/ADAMTS3/CTSS/COL6A3/COL1A2/BMP1/COL3A1/CD151/COL4A1/COL14A1/COL4A5/COL13A1/COL5A2/COL15A1/COL6A6
## 9    DCN/ELN/COL11A1/CAPN6/COL16A1/MMP2/LAMB1/COL9A3/MMP11/TIMP1/COL1A1/COL12A1/COL9A1/FN1/NID1/COL5A1/COL4A2/LAMC1/FBN2/COL2A1/COL6A2/HSPG2/COL8A1/ADAMTS1/MMP14/CTSS/COL6A3/ADAMTS9/COL1A2/FBN1/BMP1/COL3A1/BSG/A2M/COL4A1/COL14A1/COL4A5/COL13A1/COL5A2/COL15A1/COL6A6
## 10   COL11A1/P4HA2/COL16A1/COL9A3/PCOLCE/COL1A1/P3H3/COL12A1/COL9A1/P3H1/COL5A1/COL4A2/ADAMTS14/COL2A1/COL6A2/COL8A1/SERPINH1/ADAMTS3/COL6A3/COL1A2/BMP1/COL3A1/COL4A1/COL14A1/COL4A5/COL13A1/COL5A2/COL15A1/COL6A6                                                                                                                       
##    Count
## 1     19
## 2      4
## 3      8
## 4      6
## 5      3
## 6      3
## 7     83
## 8     34
## 9     41
## 10    29
```

```
dim(REACTOMEupdown)
## [1] 73 10
```

```
dotplot(REACTOMEupdown)
```

<img width="866" height="625" alt="image" src="https://github.com/user-attachments/assets/b27768c8-19eb-4c38-ad23-43c27cad0a76" />

```
dev.copy(pdf, "images/Dot Plots/REACTOMEupdown.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

# GSEA: Gene Set Enrichment Analysis:
Hallmark(H) collection have 50 curated gene sets and together they contain 7,333 human genes

```
library(msigdbr)
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") #Get Hallmark ("H") gene sets for human
dim(h_gene_sets)
## [1] 7333   23
```

```
# 7333   23

#build TERM2GENE table for GSEA

#This keeps two columns:
#gs_name (the pathway name)
#ensembl_gene ID

msigdbr_t2g <-
  as.data.frame(dplyr::distinct(h_gene_sets, gs_name, ensembl_gene)) 

head(msigdbr_t2g)
```

```
##                 gs_name    ensembl_gene
## 1 HALLMARK_ADIPOGENESIS ENSG00000165029
## 2 HALLMARK_ADIPOGENESIS ENSG00000197150
## 3 HALLMARK_ADIPOGENESIS ENSG00000167315
## 4 HALLMARK_ADIPOGENESIS ENSG00000115361
## 5 HALLMARK_ADIPOGENESIS ENSG00000117054
## 6 HALLMARK_ADIPOGENESIS ENSG00000122971
dim(msigdbr_t2g) # 7333   2
## [1] 7333    2
```

## Running GSEA
we have to create ranking list as our input for GSEA result because unlike other enrichment analysis GSEA will not take DE genes as input

```
# create the gene rank 
rank <- DEG$stat
names(rank) <- rownames(DEG) #attaches gene IDs as names to the statistic values.

# run GSEA
gsea <- clusterProfiler::GSEA(geneList = sort(rank, decreasing = TRUE), TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.1)

gsea_res <- as.data.frame(gsea)
dim(gsea_res) #28 11
## [1] 28 11
View(gsea_res)
write.table(gsea_res , file = "GSEA gene sets.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

## GSEA Enrichment plots

```
# 1st enriched Hallmark geneset:
library(enrichplot)
gseaplot(gsea, geneSetID = 1, by = "runningScore", title = "GSEA ES plot")
```

<img width="898" height="617" alt="image" src="https://github.com/user-attachments/assets/6dcd8e90-c30a-492f-a5a8-083c2a8f4967" />

```
dev.copy(pdf, "images/GSEA Plots/ GSEA plot 1.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

```
#2nd enriched Hallmark geneset.
gseaplot(gsea, geneSetID = 2, by = "runningScore", title = "GSEA ES plot")
```

<img width="903" height="622" alt="image" src="https://github.com/user-attachments/assets/bd84d82f-bc06-47a7-b5ba-1427ce46aafc" />

```
dev.copy(pdf, "images/GSEA Plots/ GSEA plot 2.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

```
#3rd enriched Hallmark geneset.
gseaplot(gsea, geneSetID = 3, by = "runningScore", title = "GSEA ES plot")
```

<img width="888" height="618" alt="image" src="https://github.com/user-attachments/assets/1f736742-6e7b-4f14-8650-6efbd17e868f" />

```
dev.copy(pdf, "images/GSEA Plots/ GSEA plot 3.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

## bar plot of the GSEA results plotting NES values
```
library(ggplot2)
ggplot(gsea, aes(NES, reorder(Description, NES), fill=qvalue),
       showCategory=20) + 
  geom_col(orientation='y') +
  scale_fill_continuous(low='red', 
                        high='blue', 
                        guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) + 
  ggtitle("GSEA Hallmark") +
  theme_minimal()
```

<img width="906" height="626" alt="image" src="https://github.com/user-attachments/assets/e89ec02d-ca64-4c13-b601-89d517153dc0" />

```
dev.copy(pdf, "images/GSEA Plots/ GSEA barplot.pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

## Heatmap of the core enriched genes of GSEA gene set
coregenes will be a vector of the genes in the leading edge subset for the first enriched pathway in our GSEA results. This accesses the first gene set in our GSEA results and retrieves the core enrichment genes, the genes that contributed most to the enrichment score.
```
coregenes <- as.vector(stringr::str_split(string = gsea@result$core_enrichment[1], pattern = "/", simplify = TRUE))
head(coregenes,10)
##  [1] "ENSG00000172638" "ENSG00000106483" "ENSG00000011465" "ENSG00000078098"
##  [5] "ENSG00000087245" "ENSG00000196924" "ENSG00000198959" "ENSG00000166147"
##  [9] "ENSG00000084636" "ENSG00000135862"
length(coregenes)
## [1] 55
```

```
library(pheatmap)
pheatmap(assay(rld)[coregenes, ], main = "Core genes")
```

<img width="887" height="623" alt="image" src="https://github.com/user-attachments/assets/129b879f-7f7b-4c6e-bac2-51010891f970" />

```
dev.copy(pdf, "images/GSEA Plots/ core genes from 1st gene set .pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

we are mapping your core genes from GSEA into readable gene symbols with the help of mapIds() function

```
library("AnnotationDbi")
library('org.Hs.eg.db')
coresymbol <- mapIds(org.Hs.eg.db,
                    keys=coregenes,
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
## 'select()' returned 1:1 mapping between keys and columns
as.vector(coresymbol)
##  [1] "EFEMP2"   "SFRP4"    "DCN"      "FAP"      "MMP2"     "FLNA"    
##  [7] "TGM2"     "FBN1"     "COL16A1"  "LAMC1"    "EMP3"     "LGALS1"  
## [13] "TPM1"     "SERPINH1" "INHBA"    "TAGLN"    "FOXC2"    "COL11A1" 
## [19] "ACTA2"    "PCOLCE"   "POSTN"    "ECM1"     "PRRX1"    "IGFBP3"  
## [25] "ITGA5"    "MYLK"     "MMP14"    "COL4A2"   "COL5A2"   "MXRA5"   
## [31] "COL6A2"   "SNAI2"    "FBLN1"    "GLIPR1"   "COL4A1"   "CCN1"    
## [37] "COL3A1"   "FBN2"     "MGP"      "LOXL2"    "TIMP1"    "SFRP1"   
## [43] "ELN"      "SERPINE1" "TPM2"     "MYL9"     "BGN"      "COL1A2"  
## [49] "COL5A1"   "COL1A1"   "TGFBI"    "COL12A1"  "FN1"      "LAMA2"   
## [55] "COL6A3"
```

```
pheatmap(assay(rld)[coregenes, ], main = "Core genes", labels_row = coresymbol)
```

<img width="863" height="623" alt="image" src="https://github.com/user-attachments/assets/42b13c13-01b5-48e9-ad5b-21b4f027649a" />

```
dev.copy(pdf, "images/GSEA Plots/ core genes symbol .pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

## Heatmap of scale data across genes

```
pheatmap(assay(rld)[coregenes, ], main = "Core genes", labels_row = coresymbol, scale = "row")
```

<img width="862" height="622" alt="image" src="https://github.com/user-attachments/assets/dbba82ab-b680-45b0-a8f8-b5b172fe0f29" />

```
dev.copy(pdf, "images/GSEA Plots/ scale data across core genes .pdf", width = 10, height = 7)
## pdf 
##   3
dev.off()
## png 
##   2
```

Choice one of the GO enrichment results and report how many categories are significant (p.adjust < 0.05) BP (biological process) enrichment result

```
GO_BP_df <- as.data.frame(GO_BP)
#View(GO_BP_df)
dim(GO_BP_df)
## [1] 672   9
sig_GO <- GO_BP_df[GO_BP_df$p.adjust < 0.05, ]
n_sig <- nrow(sig_GO)
n_sig
## [1] 672
head(sig_GO[, c("ID", "Description", "p.adjust")], 10)
##                    ID                                   Description
## GO:0030198 GO:0030198             extracellular matrix organization
## GO:0043062 GO:0043062          extracellular structure organization
## GO:0045229 GO:0045229 external encapsulating structure organization
## GO:0030199 GO:0030199                  collagen fibril organization
## GO:0001503 GO:0001503                                  ossification
## GO:0031589 GO:0031589                       cell-substrate adhesion
## GO:0048568 GO:0048568                   embryonic organ development
## GO:0051216 GO:0051216                         cartilage development
## GO:0060537 GO:0060537                     muscle tissue development
## GO:0002062 GO:0002062                   chondrocyte differentiation
##                p.adjust
## GO:0030198 3.616923e-26
## GO:0043062 3.616923e-26
## GO:0045229 3.616923e-26
## GO:0030199 2.666845e-12
## GO:0001503 5.011169e-11
## GO:0031589 7.952111e-11
## GO:0048568 8.623984e-11
## GO:0051216 8.623984e-11
## GO:0060537 8.623984e-11
## GO:0002062 3.337158e-10
```

gene present in the most enrich catogery of GO_BP enrichment result

```
GO_BP_df <- as.data.frame(GO_BP)
GO_BP_df <- GO_BP_df[order(GO_BP_df$p.adjust), ]
most_enriched_genes <- strsplit(GO_BP_df$geneID[1], "/")[[1]]
length(most_enriched_genes)
## [1] 78
```

# reads quality with FastQC:
[FastQC on disease sample](file:///C:/Users/ramis/OneDrive/Desktop/bioinformetics%201/RNA-seq%20analysis/FastQC%20and%20MultiQC%20report/FastQC_on_data_DISEASE%20sample%20__Webpage/FastQC_on_data_17__Webpage_html.html)

[FastQC on normal sample](file:///C:/Users/ramis/OneDrive/Desktop/bioinformetics%201/RNA-seq%20analysis/FastQC%20and%20MultiQC%20report/FastQC_on_data_NORMAL%20sample__Webpage/FastQC_on_data_20__Webpage_html.html)

# Generate a report with MultiQC:
[MultiQC Report](file:///C:/Users/ramis/OneDrive/Desktop/bioinformetics%201/RNA-seq%20analysis/FastQC%20and%20MultiQC%20report/MultiQC__Webpage_html.html)
