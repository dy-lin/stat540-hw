seminar05
================
Diana Lin

## Setup

Load the packages:

``` r
library(here)
library(tidyverse)
```

## Generating a Count Table

``` r
HTseq_output <- read.table(
  here("seminars", "seminar05", "hela_count.txt"),
  sep = "\t")

HTseq_output <- HTseq_output %>%
  rename(Gene = "V1", Count = "V2")

head(HTseq_output)
```

    ##                 Gene Count
    ## 1 ENSG00000000003.14     0
    ## 2  ENSG00000000005.5     0
    ## 3 ENSG00000000419.12     0
    ## 4 ENSG00000000457.13     0
    ## 5 ENSG00000000460.16    16
    ## 6 ENSG00000000938.12     0

``` r
summary(HTseq_output)
```

    ##                      Gene           Count        
    ##  __alignment_not_unique:    1   Min.   :      0  
    ##  __ambiguous           :    1   1st Qu.:      0  
    ##  __no_feature          :    1   Median :      0  
    ##  __not_aligned         :    1   Mean   :    197  
    ##  __too_low_aQual       :    1   3rd Qu.:      0  
    ##  ENSG00000000003.14    :    1   Max.   :7775281  
    ##  (Other)               :58287

## Convert to TPM

To convert counts to TPM:

1.  **Reads per Kilobase (RPK):** Divide the read counts by the length
    of the each gene in kilobases.
2.  **Reads per Kilobase per Million (RPKM):** Count up all the RPK
    values in a sample and divide this number by 1,000,000.
3.  **Transcripts per Million (TPM):** Divide the RPK values by the “per
    million” scaling factor.

*Source:
<https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/>*
