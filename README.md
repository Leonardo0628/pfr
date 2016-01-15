Across-platform imputation of DNA methylation levels using penalized functional regression
======

DNA methylation is a key epigenetic modification involved in both normal development and disease progression. Recent advances in high-throughput technologies have enabled genomewide profiling of DNA methylation. However, DNA methylation profiling often employs different designs and platforms with varying resolution, which hinders joint analysis of methylation data from multiple platforms.

Here we propose a penalized functional regression model to impute an HM27 (Illumina HumanMethylation27) dataset into an HM450 (Illumina HumanMethylation450) dataset. It is demonstrated that, by incorporating functional predictors, our model can utilize information from non-local probes and impute the missing missing methylation data accurately.

Content
-------
0. [pfr.R](pfr.R): main script used for imputation.
1. [refund_lib.R](refund_lib.R): library for penalized funuctional regression.
2. annotation/GPL8490-65.csv.gz: compressed annotation file of the HM27 array
3. annotation/GPL13534-10305.csv.gz: compressed annotation file of the HM450 array

Input Files
-----------

To run the model, you need to put two files under the same directory as the scripts:

0. meth\_450K\_QC.txt: HM450 data used to train the penalized functional regression model.
1. meth\_27K\_QC.txt: HM27 data which you want to impute into HM450 data

Notice that the methylation data in 1 and 2 should come from the same tissue. For simplicity, we assume both files are tab-delimited files where the row name is probe ID and the column name is sample ID. For example,

```
2802	2803	2804	...(more samples)...
cg03586879	0.126	0.247	0.029
cg19378133	0.294	0.578	0.037
cg03490200	0.924	0.955	0.941
...(more probes)...
```

Usage
-----

The main script is [pfr.R](pfr.R). You can either run it interactively in R or in the batch mode

```
R CMD BATCH pfr.R
```

Output Files
-----------

There will be two output files:

0. impute.txt: A tab-delimited file storing the imputed β values, where the row name is probe ID and the column name is sample ID.
1. dispersion.txt: A tab-delimited file storing the under-dispersion measure for each probe, where row name is probe ID. This can be used to remove low-quality results in 1.