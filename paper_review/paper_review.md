# Paper Review
Paper under review: [Gene Expression by Mouse Inner Ear Hair Cells during Development](paper.pdf) [[NCBI]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4405555/)

## Goal, Findings, and Conclusion

The biological problem presented in this paper is that hearing loss occurs as hair cells (HCs), responsible for hearing, lose their function. This can happen due to naturally occurring processes such as ageing or pathogenic mutations. To solve this problem, the authors propose that hearing function can be restored if the expression of genes responsible for hearing can be induced in surrounding cells (SCs). In order to identify these genes responsible for hearing, the study focuses on finding differentially expressed genes in HCs using RNA-seq data. 

This experiment begins with the separation of HCs from SCs in the cochlea and utricle, accomplished with a mouse strain with eGFP under the control of _Pou4f3_, a transcription factor of HCs. This ensures that HCs have a strong GFP signal while the SCs have a weak GFP signal when FACS sorting. RNA was extracted from each cell in each tissue at 4 different developmental stages and then sequenced. 

The authors found that there are differentially expressed genes in the HCs with known and unknown functions. Those of known function could potentially restore hearing if expression is induced in SCs. Those with  unknown function that are potential candidates for new deafness genes in HCs, especially those genes uniquely expressed in HCs.

## Datasets: Data Matrix and Metadata

This dataset consisted of 16 samples:
* 4 HCs from the utricle (1 from E16, 1 from P0, 1 from P4, 1 from P7)
* 4 SCs from the utricle (1 from E16, 1 from P0, 1 from P4, 1 from P7)
* 4 HCs from the cochlea (1 from E16, 1 from P0, 1 from P4, 1 from P7)
* 4 SCs from the cochlea (1 from E16, 1 from P0, 1 from P4, 1 from P7)

The combined metadata/data matrix should have looked something like this (4 samples shown in long format for example):

Cell Type| Tissue Type | Developmental Stage | Gene | `baseMean` |  `log2FoldChange` | `lfcSE` | `stat` | `pval` | `padj`
---------|-------------|---------------------|-------|----------|--------------------|---------|--------|--------|--------
HC|cochlea|E16|geneX|X|X|X|X|X|X|X
HC|cochlea|P0|geneX|X|X|X|X|X|X|X
HC|cochlea|P4|geneX|X|X|X|X|X|X|X
HC|cochlea|P7|geneX|X|X|X|X|X|X|X
...|...|...|...|...|...|...|...|...|...|...

#### Legend
* `baseMean`: normalized expression values
* `log2FoldChange`: log<sub>2</sub> fold change
* `lfcSE`: log<sub>2</sub> fold change standard error
* `stat`: test statistic
* `pval`: p-value
* `padj`: adjusted p-value (FDR)

After RNA extraction and sequencing, then the RNA-seq reads were aligned to reference genome mm9, quantified into ‘read counts’, normalized, and statistically tested with `DESeq` (outputs the matrix above, but without first 3 columns of metadata).

In this case, the dataset represents the gene expression of HCs and SCs (cell types) from the utricle and cochlea (tissue type) at 4 different developmental stages (time points): E16, P0, P4, P7, which are representative of before, during, and after the HC acquisition of mechanosensitivity.

## Analytical Steps

The main types of analyses that were conducted on this dataset include:
* Differential Expression Analysis (DEA) using the `DESeq` R packages
* Principal Components Analysis (PCA) using the `stats` R package 
* Hierarchical Clustering using `dChip`
* Functional Enrichment Analysis using `DAVID`

`DESeq` normalized the data by calculating the geometric mean across all samples for each gene, then dividing each count by this geometric mean. Then multiple test correction was performed using the Benjamini-Hochberg method of ranking p-values to determine the appropriate FDR for each gene.

PCA was performed to identify inherent patterns or clustering in the dataset. To reduce complexity, genes that had < 10 total reads across all developmental stages were excluded, and numbers < 1 at a single developmental stage were rounded to 1. TThis reduced matrix containing normalized expression values from `DESeq` was then log<sub>2</sub> transformed before running `prcomp()`.

Hierarchical clustering was conducted to group and visualize genes according to their spatial and temporal expression patterns, especially genes enriched in either cell type or tissue type. In this process, the normalized expression values from `DESeq` were log transformed and standardized for each gene across all samples. These values for a single gene cluster were then averaged, and used to calculate distance between gene clusters.

Functional Enrichment Analysis was done to see if there was a functional enrichment of certain biological processes within a target list of genes, using biological process ontologies. A biological process was deemed 'enriched' if its genes had a total read count > 50 and ratio of the difference in cell types over total number of cells was > 0.67.

Control for variation between each group is usually done by having biological or technical replicates which were absent. However, lots of validation steps were conducted, such as the tdTomato experiment (with _Gfi1_ promoter) as well as the RT-PCR, qPCR, _in situ_ hybridization, etc. To simulate biological replicates, samples of different levels in other factors were combined and treated as a replicate to facilitate statistical analysis with strict FDR cutoffs.

## Comments and Suggestions

### Analytical Steps 

While I agree with most of their statistical methods, I do not understand why biological and technical replicates weren’t done in the first place, rather than trying to simulate replicates later. It is universally known that replicates are very important in statistical analysis. If replicates were skipped due to financial reasons, they should have allocated more resources for conducting replicates and less for validation.

The gene function/network analyses conducted were minimal at best. The gene network analysis was conducted using the hierarchical analysis only, where it is assumed they are part of the same network if they are in the same cluster. Functional Enrichment Analysis was very broadly conducted using the biological process ontologies. However, in Gene Ontology [(GO)](http://geneontology.org/docs/ontology-documentation/), two other big categories exist: Molecular Function and Cellular Component, both of which were completely ignored in this study.

Addtionally, `DESeq2` was published in _BMC Genome Biology_ on December 5, 2014, whereas this paper was submitted to the _Journal of Neuroscience_ on December 17, 2014. If the work is as reproducible as they say it is, they should have rerun their analysis using `DESeq2`, replacing  the outdated `DESeq` and seeing if they get comparable results, and if not, why.

### Interpretation of Findings (optional)

I somewhat agree with how the authors interpreted their findings. There isn’t anything ‘wrong’ with their interpretation, but to me, it seems like a very broad conclusion, that isn’t entirely impactful to the field: 

> In summary, we performed the first comprehensive gene expression study specific for vestibular and cochlear HCs during mouse development. The public database generated is a powerful tool for discovering and understanding proteins involved in HC function. We anticipate these data will lead to the discovery of new genes important for HC mechanotransduction and new deafness genes in mouse and human.

In my opinion, these results _alone_ should not have been published-- some downstream work testing some of the specific new and unknown function HC genes that they found to be differentially expressed would be more interesting and impactful. For example, if they had conducted an experiment where they knocked out one of those unknown function genes _in vivo_, to observe some phenotypic changes it would have been more interesting. To me, this paper seems to lack _novelty_, and most of the paper is validation and confirmation of what is already known.

### Data Accessibility (optional)
Overall, I think this study was very well documented, especially with the many _wet-lab_ quality control and validation steps. At each step in the experiment, another side-experiment was performed to verify and validate, for example, checking for differential expression in genes already known to be differentially expressed. However, the quality control and validation in the dry-lab portion of the study was lacking. The authors state that their work is reproducible as demonstrated by their _wet-lab_ validations, but on the _bioinformatics_ side, almost nothing is provided. Generally with R statistical analysis, R scripts and other custom scripts should be provided as supplementary materials, as well as parameter specifications for all programs executed. 

The data provided is said to be provided in GEO accession [GSE60019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60019), but actually isn't present when you use `getGEO("GSE60019")` of the `GEOquery` R package. In fact the raw RNA-seq reads are provided, but in the [SRA](https://www.ncbi.nlm.nih.gov/sra?term=SRP045182) not GEO. However, in these studies, the matrix with the normalized read counts from `DESeq` is usually provided as well. There is a [supplementary text file](GSE60019_Processed_Data.txt) provided in GEO, but it is poorly annotated, with no indication of what some of the column names represent, and no indication of what the numeric values are. 

This experiment may be highly reproducible in the _wet-lab_, but given their provided resources, it is nearly impossible to do so in the _dry-lab_.


