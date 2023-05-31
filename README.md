# endest

<!-- badges: start -->
<!-- badges: end -->

The goal of endest is to estimate uterine cycle time of endometrium samples using expression data. This expression data is usually in the form of RNA-seq counts, but other platforms such as microarray expression values can also be used.

The estimated cycle time is given by a number between 0 and 99, where 0 is approximately the start of menstruation, ~8 is the start of the proliferative phase, and ~58 is the start of the secretory phase.

Endest is published pursuant to the terms located in the [license file](LICENSE) which permit reproduction, publication and adaptation of endest solely for non-commercial purposes.

## Installation

You can install endest using:

``` r
# First install the preprocessCore package dependency from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")

# Install the endest package
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jessicachung/endest")
```

## Usage

Using example RNA-seq data:

``` r
# Load endest package
library(endest)

# Load example RNA-seq count data
data(example_rna_exprs)

head(example_rna_exprs)
#>                      RS_1     RS_2     RS_3
#> ENSG00000135862 9.7109998 8.568146 9.722390
#> ENSG00000255690 0.6642925 2.776653 1.968152
#> ENSG00000142089 7.1937165 8.564381 8.472562
#> ENSG00000089220 6.9412745 8.118290 7.892754
#> ENSG00000108932 7.4077372 2.068414 4.304426
#> ENSG00000147533 5.2366574 5.570384 4.630944

# Estimate cycle time
results <- estimate_cycle_time(exprs=example_rna_exprs)
#> Using ensembl identifiers from rownames.
#> Using 2000 / 2000 observations from the input expression matrix.

# View estimated cycle times
results$estimated_time
#> RS_1 RS_2 RS_3 
#>    1   31   83
```

Using example microarray data:

``` r
# Load endest package
library(endest)

# Load example array data and Entrez IDs
data(example_array_exprs)
data(example_array_entrez_ids)

head(example_array_exprs)
#>                  MA_1     MA_2     MA_3
#> ILMN_1725414 4.596832 4.373430 4.463814
#> ILMN_1704163 4.319510 4.899428 4.584553
#> ILMN_1726755 9.584955 9.421494 9.377625
#> ILMN_1760593 7.498393 7.174511 7.402563
#> ILMN_1715653 4.720311 4.401360 4.212626
#> ILMN_1677636 7.773338 4.686772 5.413982

head(example_array_entrez_ids)
#> ILMN_1725414 ILMN_1704163 ILMN_1726755 ILMN_1760593 ILMN_1715653 ILMN_1677636 
#>        "368"       "6783"      "51138"       "1407"     "169693"       "1311"

# Estimate cycle time
results <- estimate_cycle_time(exprs=example_array_exprs, entrez_ids=example_array_entrez_ids)
#> Using 4669 / 5000 observations from the input expression matrix.
#> Consolidating 592 genes that have multiple observations using the mean.

# View estimated cycle times
results$estimated_time
#> MA_1 MA_2 MA_3 
#>    6   54   71
```

Example using an external dataset in the GEO database:

``` r
# Load packages
library(endest)

# Download GSE153740 RNA-seq gene counts from GREIN.
# Navigate to the count table at http://www.ilincs.org/apps/grein/?gse=GSE153740
# and download raw gene level counts.

# Load in data
file_path <- "/path/to/counts/GSE153740_GeneLevel_Raw_data.csv"
raw_counts <- read.csv(file_path, stringsAsFactors=FALSE)

# Format into matrix
counts <- as.matrix(raw_counts[,-c(1,2)])
rownames(counts) <- raw_counts[,1, drop=TRUE]
head(counts)
#>                 GSM4650386 GSM4650387 GSM4650388 GSM4650389 GSM4650390
#> ENSG00000000003        452        400        776        484        613
#> ENSG00000000005          0          0          1          2          0
#> ENSG00000000419        779        511        908        890        661
#> ENSG00000000457        715        458        487        387        468
#> ENSG00000000460        149        146        111         98        184
#> ENSG00000000938         42         75        358        458         52
#>                 GSM4650391 GSM4650392 GSM4650393
#> ENSG00000000003        802       1112        597
#> ENSG00000000005          0          0          0
#> ENSG00000000419        919        715        747
#> ENSG00000000457        629        445        423
#> ENSG00000000460        139        110        161
#> ENSG00000000938        379        638         85

# Estimate cycle time
results <- estimate_cycle_time(exprs=counts)
#> Using ensembl identifiers from rownames.
#> Using 16969 / 35413 observations from the input expression matrix.
results$estimated_time
#> GSM4650386 GSM4650387 GSM4650388 GSM4650389 GSM4650390 GSM4650391 GSM4650392 
#>         70         72         84         82         60         85         94 
#> GSM4650393 
#>         61
```
