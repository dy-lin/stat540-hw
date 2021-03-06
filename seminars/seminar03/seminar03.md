seminar03
================
Diana Lin
27/01/2020

## Deliverables

  - [Part 2 Exercise](#part-2-exercise)
  - [Part 3 Exercise](#part-3-exercise)
  - [Part 4 Exercise](#part-4-exercise)
      - [Density Plot](#density-plot)
      - [Exclusion of Transcript
        Lengths](#exclusion-of-transcript-lengths)
      - [Gene with Lowest P-value](#gene-with-lowest-p-value)

## Part 1 - Accessing data using GEOquery

Load the libraries:

``` r
library(GEOquery)
library(biomaRt)
library(tidyverse)
library(data.table)
library(reshape2)
```

Get the GEO
    dataset:

``` r
gds <- getGEO("GDS507")
```

    ## File stored at:

    ## /var/folders/n4/d75pxyjx0nn_ttn4y5fn0vgw0000gn/T//RtmpBRD5Nx/GDS507.soft.gz

    ## Parsed with column specification:
    ## cols(
    ##   ID_REF = col_character(),
    ##   IDENTIFIER = col_character(),
    ##   GSM11815 = col_double(),
    ##   GSM11832 = col_double(),
    ##   GSM12069 = col_double(),
    ##   GSM12083 = col_double(),
    ##   GSM12101 = col_double(),
    ##   GSM12106 = col_double(),
    ##   GSM12274 = col_double(),
    ##   GSM12299 = col_double(),
    ##   GSM12412 = col_double(),
    ##   GSM11810 = col_double(),
    ##   GSM11827 = col_double(),
    ##   GSM12078 = col_double(),
    ##   GSM12099 = col_double(),
    ##   GSM12269 = col_double(),
    ##   GSM12287 = col_double(),
    ##   GSM12301 = col_double(),
    ##   GSM12448 = col_double()
    ## )

Peek into the structure of the GEO object:

``` r
str(gds)
```

    ## Formal class 'GDS' [package "GEOquery"] with 3 slots
    ##   ..@ gpl      :Formal class 'GPL' [package "GEOquery"] with 2 slots
    ##   .. .. ..@ dataTable:Formal class 'GEODataTable' [package "GEOquery"] with 2 slots
    ##   .. .. .. .. ..@ columns:'data.frame':  0 obs. of  0 variables
    ##   .. .. .. .. ..@ table  :'data.frame':  0 obs. of  0 variables
    ##   .. .. ..@ header   : list()
    ##   ..@ dataTable:Formal class 'GEODataTable' [package "GEOquery"] with 2 slots
    ##   .. .. ..@ columns:'data.frame':    17 obs. of  4 variables:
    ##   .. .. .. ..$ sample       : Factor w/ 17 levels "GSM11810","GSM11815",..: 2 4 5 7 9 10 12 14 16 1 ...
    ##   .. .. .. ..$ disease.state: Factor w/ 2 levels "normal","RCC": 2 2 2 2 2 2 2 2 2 1 ...
    ##   .. .. .. ..$ individual   : Factor w/ 10 levels "001","005","011",..: 6 4 1 2 3 5 8 9 10 6 ...
    ##   .. .. .. ..$ description  : chr [1:17] "Value for GSM11815: C035 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear "| __truncated__ "Value for GSM11832: C023 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear "| __truncated__ "Value for GSM12069: C001 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear "| __truncated__ "Value for GSM12083: C005 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear "| __truncated__ ...
    ##   .. .. ..@ table  :'data.frame':    22645 obs. of  19 variables:
    ##   .. .. .. ..$ ID_REF    : chr [1:22645] "200000_s_at" "200001_at" "200002_at" "200003_s_at" ...
    ##   .. .. .. ..$ IDENTIFIER: chr [1:22645] "PRPF8" "CAPNS1" "RPL35" "MIR6805" ...
    ##   .. .. .. ..$ GSM11815  : num [1:22645] 4254 17996 41679 65391 19030 ...
    ##   .. .. .. ..$ GSM11832  : num [1:22645] 5298 12011 39117 34806 15814 ...
    ##   .. .. .. ..$ GSM12069  : num [1:22645] 4026 10284 38759 31257 16356 ...
    ##   .. .. .. ..$ GSM12083  : num [1:22645] 3498 2535 32848 28308 9580 ...
    ##   .. .. .. ..$ GSM12101  : num [1:22645] 3566 11048 39634 67448 14274 ...
    ##   .. .. .. ..$ GSM12106  : num [1:22645] 4903 13354 43511 56990 17217 ...
    ##   .. .. .. ..$ GSM12274  : num [1:22645] 6373 8564 46857 57972 19117 ...
    ##   .. .. .. ..$ GSM12299  : num [1:22645] 4829 17248 47032 57570 17488 ...
    ##   .. .. .. ..$ GSM12412  : num [1:22645] 5206 16018 22152 29062 14672 ...
    ##   .. .. .. ..$ GSM11810  : num [1:22645] 2757 6077 26661 35141 17733 ...
    ##   .. .. .. ..$ GSM11827  : num [1:22645] 3932 15704 26374 23629 18022 ...
    ##   .. .. .. ..$ GSM12078  : num [1:22645] 3730 10138 23810 22100 17957 ...
    ##   .. .. .. ..$ GSM12099  : num [1:22645] 3223 11614 24749 21651 15958 ...
    ##   .. .. .. ..$ GSM12269  : num [1:22645] 3640 8460 21937 18551 15800 ...
    ##   .. .. .. ..$ GSM12287  : num [1:22645] 4886 10283 31463 23496 16686 ...
    ##   .. .. .. ..$ GSM12301  : num [1:22645] 4070 11844 22734 21315 18817 ...
    ##   .. .. .. ..$ GSM12448  : num [1:22645] 3482 9742 25396 28631 17421 ...
    ##   .. .. .. ..- attr(*, "spec")=
    ##   .. .. .. .. .. cols(
    ##   .. .. .. .. ..   ID_REF = col_character(),
    ##   .. .. .. .. ..   IDENTIFIER = col_character(),
    ##   .. .. .. .. ..   GSM11815 = col_double(),
    ##   .. .. .. .. ..   GSM11832 = col_double(),
    ##   .. .. .. .. ..   GSM12069 = col_double(),
    ##   .. .. .. .. ..   GSM12083 = col_double(),
    ##   .. .. .. .. ..   GSM12101 = col_double(),
    ##   .. .. .. .. ..   GSM12106 = col_double(),
    ##   .. .. .. .. ..   GSM12274 = col_double(),
    ##   .. .. .. .. ..   GSM12299 = col_double(),
    ##   .. .. .. .. ..   GSM12412 = col_double(),
    ##   .. .. .. .. ..   GSM11810 = col_double(),
    ##   .. .. .. .. ..   GSM11827 = col_double(),
    ##   .. .. .. .. ..   GSM12078 = col_double(),
    ##   .. .. .. .. ..   GSM12099 = col_double(),
    ##   .. .. .. .. ..   GSM12269 = col_double(),
    ##   .. .. .. .. ..   GSM12287 = col_double(),
    ##   .. .. .. .. ..   GSM12301 = col_double(),
    ##   .. .. .. .. ..   GSM12448 = col_double()
    ##   .. .. .. .. .. )
    ##   ..@ header   :List of 23
    ##   .. ..$ channel_count           : chr "1"
    ##   .. ..$ dataset_id              : chr [1:12] "GDS507" "GDS507" "GDS507" "GDS507" ...
    ##   .. ..$ description             : chr [1:13] "Investigation into mechanisms of renal clear cell carcinogenesis (RCC). Comparison of renal clear cell tumor ti"| __truncated__ "RCC" "normal" "035" ...
    ##   .. ..$ email                   : chr "geo@ncbi.nlm.nih.gov"
    ##   .. ..$ feature_count           : chr "22645"
    ##   .. ..$ institute               : chr "NCBI NLM NIH"
    ##   .. ..$ name                    : chr "Gene Expression Omnibus (GEO)"
    ##   .. ..$ order                   : chr "none"
    ##   .. ..$ platform                : chr "GPL97"
    ##   .. ..$ platform_organism       : chr "Homo sapiens"
    ##   .. ..$ platform_technology_type: chr "in situ oligonucleotide"
    ##   .. ..$ pubmed_id               : chr "14641932"
    ##   .. ..$ ref                     : chr "Nucleic Acids Res. 2005 Jan 1;33 Database Issue:D562-6"
    ##   .. ..$ reference_series        : chr "GSE781"
    ##   .. ..$ sample_count            : chr "17"
    ##   .. ..$ sample_id               : chr [1:12] "GSM11815,GSM11832,GSM12069,GSM12083,GSM12101,GSM12106,GSM12274,GSM12299,GSM12412" "GSM11810,GSM11827,GSM12078,GSM12099,GSM12269,GSM12287,GSM12301,GSM12448" "GSM11810,GSM11815" "GSM11827,GSM11832" ...
    ##   .. ..$ sample_organism         : chr "Homo sapiens"
    ##   .. ..$ sample_type             : chr "RNA"
    ##   .. ..$ title                   : chr "Renal clear cell carcinoma (HG-U133B)"
    ##   .. ..$ type                    : chr [1:13] "Expression profiling by array" "disease state" "disease state" "individual" ...
    ##   .. ..$ update_date             : chr "Mar 04 2004"
    ##   .. ..$ value_type              : chr "count"
    ##   .. ..$ web_link                : chr "http://www.ncbi.nlm.nih.gov/geo"

Extract the metadata table and gene expression table:

``` r
meta_data <- data.frame(
  Sample = gds@dataTable@columns$sample,
  disease = gds@dataTable@columns$disease.state
)

gds_data <- gds@dataTable@table
```

## Part 2 - Exploring a gene expression dataset

Peek into the
    data:

``` r
head(gds_data)
```

    ##        ID_REF IDENTIFIER GSM11815 GSM11832 GSM12069 GSM12083 GSM12101 GSM12106
    ## 1 200000_s_at      PRPF8   4254.0   5298.2   4026.5   3498.4   3566.4   4903.1
    ## 2   200001_at     CAPNS1  17996.2  12010.7  10283.5   2534.7  11048.4  13354.0
    ## 3   200002_at      RPL35  41678.8  39116.9  38758.9  32847.7  39633.9  43511.2
    ## 4 200003_s_at    MIR6805  65390.9  34806.2  31257.2  28308.5  67447.5  56989.9
    ## 5   200004_at     EIF4G2  19030.1  15813.6  16355.7   9579.7  14273.5  17217.0
    ## 6   200005_at      EIF3D   8824.5   9706.2  10590.0   6986.7   9400.4  12835.2
    ##   GSM12274 GSM12299 GSM12412 GSM11810 GSM11827 GSM12078 GSM12099 GSM12269
    ## 1   6372.6   4829.1   5205.8   2756.8   3932.0   3729.9   3223.4   3640.5
    ## 2   8563.8  17247.6  16018.5   6077.0  15703.8  10138.5  11614.4   8460.5
    ## 3  46856.7  47032.4  22152.2  26660.7  26373.6  23809.6  24749.3  21936.8
    ## 4  57972.5  57570.5  29062.2  35140.9  23629.3  22100.5  21651.0  18550.7
    ## 5  19116.9  17487.6  14671.6  17733.1  18022.4  17957.4  15958.0  15799.8
    ## 6  10299.0  12375.2   7645.4   8661.5   7355.7   6973.4   6855.9   7949.2
    ##   GSM12287 GSM12301 GSM12448
    ## 1   4886.3   4070.2   3482.1
    ## 2  10282.6  11844.3   9741.6
    ## 3  31462.8  22733.7  25395.5
    ## 4  23496.5  21315.4  28631.4
    ## 5  16685.8  18817.3  17421.1
    ## 6   9486.5   7494.5   7252.1

How many rows?

``` r
nrow(gds_data)
```

    ## [1] 22645

How many columns?

``` r
ncol(gds_data)
```

    ## [1] 19

Computing the average count in each
    sample:

``` r
apply(gds_data[,-c(1, 2)], 2, median)
```

    ## GSM11815 GSM11832 GSM12069 GSM12083 GSM12101 GSM12106 GSM12274 GSM12299 
    ##    265.6    250.3    218.5    309.7    281.9    240.1    280.2    217.0 
    ## GSM12412 GSM11810 GSM11827 GSM12078 GSM12099 GSM12269 GSM12287 GSM12301 
    ##    264.4    273.8    264.6    266.5    269.3    288.6    238.7    244.5 
    ## GSM12448 
    ##    264.3

The first and second columns are excluded because they hold the probe
and gene names.

Tidy our
data:

``` r
melted_data <- melt(gds_data, id.vars = c("ID_REF", "IDENTIFIER"), var = "Sample")
head(melted_data)
```

    ##        ID_REF IDENTIFIER   Sample   value
    ## 1 200000_s_at      PRPF8 GSM11815  4254.0
    ## 2   200001_at     CAPNS1 GSM11815 17996.2
    ## 3   200002_at      RPL35 GSM11815 41678.8
    ## 4 200003_s_at    MIR6805 GSM11815 65390.9
    ## 5   200004_at     EIF4G2 GSM11815 19030.1
    ## 6   200005_at      EIF3D GSM11815  8824.5

Calcualte the mean gene expression per sample:

``` r
melted_data %>%
  group_by(Sample) %>%
  summarize(mean = mean(value))
```

    ## # A tibble: 17 x 2
    ##    Sample    mean
    ##    <fct>    <dbl>
    ##  1 GSM11815  751.
    ##  2 GSM11832  742.
    ##  3 GSM12069  748.
    ##  4 GSM12083  735.
    ##  5 GSM12101  803.
    ##  6 GSM12106  744.
    ##  7 GSM12274  761.
    ##  8 GSM12299  802.
    ##  9 GSM12412  685.
    ## 10 GSM11810  765.
    ## 11 GSM11827  780.
    ## 12 GSM12078  774.
    ## 13 GSM12099  766.
    ## 14 GSM12269  710.
    ## 15 GSM12287  791.
    ## 16 GSM12301  770.
    ## 17 GSM12448  757.

Calculate mean of each probe’s expression:

``` r
(new_melted_data <- melted_data %>%
  group_by(Sample, IDENTIFIER) %>%
  summarize(Count = mean(value)))
```

    ## # A tibble: 279,905 x 3
    ## # Groups:   Sample [17]
    ##    Sample   IDENTIFIER   Count
    ##    <fct>    <chr>        <dbl>
    ##  1 GSM11815 --Control   8139. 
    ##  2 GSM11815 222968_at    102. 
    ##  3 GSM11815 223641_at    200. 
    ##  4 GSM11815 224429_x_at 2385. 
    ##  5 GSM11815 224438_at     32.1
    ##  6 GSM11815 225714_s_at  291. 
    ##  7 GSM11815 225934_at    284. 
    ##  8 GSM11815 226014_at     66.3
    ##  9 GSM11815 226061_s_at   45.1
    ## 10 GSM11815 226138_s_at   23.3
    ## # … with 279,895 more rows

Use `biomaRt` to look up the chromosomal location of each gene:

``` r
# open connection between biomaRt and R.
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# function that takes in data frame, and outputs same data frame with associated chromosome annotations.
identify_gene_names <- function(df) {
  names(df) <- c("Sample", "hgnc_symbol", "Count")
  names <-
    getBM(
      attributes = c("hgnc_symbol", "chromosome_name") ,
      filters = "hgnc_symbol",
      values = unique(df$hgnc_symbol),
      mart = human
    )
  left_join(df, names, by = "hgnc_symbol")
}

# filter out all genes with annotations that are not numeric numbers between 1 and 23, X or Y.
data_with_chromosome <- identify_gene_names(new_melted_data) %>%
  filter(chromosome_name %in% c(1:23, "X", "Y"))
```

### Part 2 Exercise

**Instructions:** Modify the above code to also identify the length of
each gene captured in the dataset we have been working with in the above
exercises. This can be done by adding “transcript\_length” as attribute
in getBM function. You should end up with an extra column for
“transcript length”. We will use this number later.

I have changed the function name `identify_gene_names` to
`identify_gene_names_tx` and variable name `data_with_chromosome` to
`data_with_chromosome_tx` as to not affect the dataset for the rest of
the workflow:

``` r
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# function that takes in data frame, and outputs same data frame with associated chromosome annotations.
identify_gene_names_tx <- function(df) {
  names(df) <- c("Sample", "hgnc_symbol", "Count")
  names <-
    getBM(
      attributes = c("hgnc_symbol", "chromosome_name", "transcript_length") ,
      filters = "hgnc_symbol",
      values = unique(df$hgnc_symbol),
      mart = human
    )
  left_join(df, names, by = "hgnc_symbol")
}

# filter out all genes with annotations that are not numeric numbers between 1 and 23, X or Y.
data_with_chromosome_tx <- identify_gene_names_tx(new_melted_data) %>%
  filter(chromosome_name %in% c(1:23, "X", "Y"))

full_data_tx <- left_join(data_with_chromosome_tx, meta_data, by = "Sample")
```

    ## Warning: Column `Sample` joining factors with different levels, coercing to
    ## character vector

Here you can see the new column `transcript_length`:

``` r
data_with_chromosome_tx
```

    ## # A tibble: 1,349,171 x 5
    ## # Groups:   Sample [17]
    ##    Sample   hgnc_symbol Count chromosome_name transcript_length
    ##    <fct>    <chr>       <dbl> <chr>                       <int>
    ##  1 GSM11815 A1BG         191. 19                           2134
    ##  2 GSM11815 A1BG         191. 19                           3382
    ##  3 GSM11815 A1BG         191. 19                           2301
    ##  4 GSM11815 A1BG         191. 19                            475
    ##  5 GSM11815 A1BG         191. 19                            917
    ##  6 GSM11815 A1BG-AS1      53  19                           1718
    ##  7 GSM11815 A1BG-AS1      53  19                            959
    ##  8 GSM11815 A1BG-AS1      53  19                            559
    ##  9 GSM11815 A1BG-AS1      53  19                            977
    ## 10 GSM11815 A1BG-AS1      53  19                            712
    ## # … with 1,349,161 more rows

Let’s look at how the average expression of genes on the X chromosome
changes between RCC and normal cells.

First, combining the metadata file (`meta_data`) with expression table
(`data_with_chromosome`)
    data:

``` r
full_data <- left_join(data_with_chromosome, meta_data, by = "Sample")
```

    ## Warning: Column `Sample` joining factors with different levels, coercing to
    ## character vector

Group all samples by disease status, filter out non-X-chromosome genes,
and then calculate the mean:

``` r
full_data %>%
  group_by(disease) %>%
  filter(chromosome_name == "X") %>%
  summarize(mean = mean(Count))
```

    ## # A tibble: 2 x 2
    ##   disease  mean
    ##   <fct>   <dbl>
    ## 1 normal   686.
    ## 2 RCC      658.

## Part 3: Graphing expression data

Randomly sample 100 probes:

``` r
set.seed(5747540)
sample_to_choose <- sample(1:length(unique(full_data$hgnc_symbol)), size = 100)
names_to_choose <- as.character(unique(full_data$hgnc_symbol)[sample_to_choose])
```

Plot:

``` r
full_data %>% 
    filter(hgnc_symbol %in% names_to_choose) %>% 
    group_by(Sample) %>% 
    ggplot(aes(x = as.factor(chromosome_name), y = Count)) + geom_point()
```

![](seminar03_files/figure-gfm/plot-1.png)<!-- -->

### Part 3 Exercise

**Instructions:** By adding one additional function to the code above,
calculate the sum of all counts in each sample and divide each
expression value by that sum (hint: use mutate). Remember, you can add
multiple new columns using mutate by separating each column with a comma
(i.e mutate(x = c(“a”, “b”), y = c(“d”, “c”))). Plot this new
transformed column.

``` r
full_data %>%
  filter(hgnc_symbol %in% names_to_choose) %>%
  group_by(Sample) %>%
  mutate(sum_counts = sum(Count), normalized_count = Count/sum_counts) %>%
  ggplot(aes(x = as.factor(chromosome_name), y = normalized_count)) + geom_point()
```

![](seminar03_files/figure-gfm/mutate-1.png)<!-- -->

## Part 4 - Analyzing the results of statistical tests

Perform a T-test for each gene:

``` r
full_data %>% 
  group_by(hgnc_symbol) %>% 
  summarize( pvalue = t.test(Count ~ disease)$p.value)
```

    ## # A tibble: 9,236 x 2
    ##    hgnc_symbol pvalue
    ##    <chr>        <dbl>
    ##  1 A1BG        0.708 
    ##  2 A1BG-AS1    0.0366
    ##  3 A1CF        0.132 
    ##  4 A2MP1       0.0245
    ##  5 AADACL2     0.143 
    ##  6 AADAT       0.0304
    ##  7 AAGAB       0.469 
    ##  8 AAK1        0.0229
    ##  9 AARS2       0.0416
    ## 10 AASDH       0.0743
    ## # … with 9,226 more rows

### Part 4 Exercise

#### Density Plot

**Instructions:** Make a density plot using geom\_density() graph of the
p-value distributions of the above t-test.

``` r
full_data_p <- full_data %>% 
  group_by(hgnc_symbol) %>% 
  summarize( pvalue = t.test(Count ~ disease)$p.value)

full_data_p %>%
  ggplot(aes(x = pvalue)) +
  geom_density() +
  ggtitle("The RIGHT Plot")
```

![](seminar03_files/figure-gfm/p4-1.png)<!-- -->

#### Exclusion of Transcript Lengths

**Instructions:** Note that if you acquired transcript lengths, you
should NOT be using that data frame for this task. Can you see why?

In the data frame with `transcript_length`, there will be same rows of
gene-counts pairs repeated for each isoform of each gene. This will
wrongfully “punish” genes with higher number of isoforms if you
calculate `mean(Count)`.

For example, here is one sample and gene that has multiple entries in
the `transcript_length` data frame.

Here is the sample `GSM11815` and gene `A1BG` without
`transcript_length`:

``` r
full_data %>%
  filter(Sample == "GSM11815", hgnc_symbol == "A1BG")
```

    ## # A tibble: 1 x 5
    ## # Groups:   Sample [1]
    ##   Sample   hgnc_symbol Count chromosome_name disease
    ##   <chr>    <chr>       <dbl> <chr>           <fct>  
    ## 1 GSM11815 A1BG         191. 19              RCC

Here is the sample `GSM11815` and gene `A1BG` **with**
`transcript_length`:

``` r
full_data_tx %>%
  filter(Sample == "GSM11815", hgnc_symbol == "A1BG")
```

    ## # A tibble: 5 x 6
    ## # Groups:   Sample [1]
    ##   Sample   hgnc_symbol Count chromosome_name transcript_length disease
    ##   <chr>    <chr>       <dbl> <chr>                       <int> <fct>  
    ## 1 GSM11815 A1BG         191. 19                           2134 RCC    
    ## 2 GSM11815 A1BG         191. 19                           3382 RCC    
    ## 3 GSM11815 A1BG         191. 19                           2301 RCC    
    ## 4 GSM11815 A1BG         191. 19                            475 RCC    
    ## 5 GSM11815 A1BG         191. 19                            917 RCC

That one sample and gene combination now has 5 rows instead of 1. With
this ‘duplication’, the `mean(Count)` function is now biased against
those with more isoforms.

I will demonstrate by plotting the **WRONG** plot, using the
`transcript_length` plot:

``` r
full_data_tx %>% 
  group_by(hgnc_symbol) %>% 
  summarize( pvalue = t.test(Count ~ disease)$p.value) %>%
  ggplot(aes(x = pvalue)) +
  geom_density() +
  ggtitle("The WRONG Plot")
```

![](seminar03_files/figure-gfm/why-1.png)<!-- -->

#### Gene with Lowest P-value

**Instructions:** Also, extract a data frame of all genes with p-values
lower than 0.05. Finally, extract the name of the gene with the lowest
p-value.

Here is the data frame with all genes with p-values lower than 0.05:

``` r
(full_data_p_sig <- full_data_p %>%
  filter(pvalue < 0.05))
```

    ## # A tibble: 2,468 x 2
    ##    hgnc_symbol   pvalue
    ##    <chr>          <dbl>
    ##  1 A1BG-AS1    0.0366  
    ##  2 A2MP1       0.0245  
    ##  3 AADAT       0.0304  
    ##  4 AAK1        0.0229  
    ##  5 AARS2       0.0416  
    ##  6 ABCB1       0.00351 
    ##  7 ABCB10      0.0302  
    ##  8 ABCC3       0.0342  
    ##  9 ABCC6P1     0.000528
    ## 10 ABCG1       0.00795 
    ## # … with 2,458 more rows

To find the gene with the lowest p-value, we first sort p-value from
lowest to highest using `arrange()`, and take the first line using
`head()`, and then select the gene name field using `select()`:

``` r
# get lowest pvalue gene name
lowest_gene_p <- full_data_p_sig %>%
  arrange(pvalue) %>%
  head(n = 1L) %>%
  select(hgnc_symbol) %>%
  as.character()

# get lowest pvalue
lowest_p <- full_data_p_sig %>%
  select(pvalue) %>%
  min()

# get lowest p-value row
full_data_p_sig %>%
  arrange(pvalue) %>%
  head(n = 1L)
```

    ## # A tibble: 1 x 2
    ##   hgnc_symbol        pvalue
    ##   <chr>               <dbl>
    ## 1 CLDN2       0.00000000591

**Conclusion**: This gene with the lowest p-value is CLDN2, with a
p-value of 5.911385410^{-9}.
