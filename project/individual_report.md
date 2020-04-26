# Individual Report

## Contributions of Each Group Member

> Describes the tasks and contributions of each of your group members.

### Nairuz Elazzabi
Originally, in the `project_proposal.md`, Nairuz was assigned to do the statistical analysis using `limma`, but in the end she could not get it to work. Eventually, we learned that this course was the first time that Nairuz was using GitHub _and_ `R`. In that case, she was reassigned into the quality control portion of the work, doing the p-value histogram and density plot to check for anomalies in our data. Unfortunately we did not learn about Nairuz's lack of R and GitHub experience until too late, which is why we had minimal `limma` results in our `progress_report.md`. She failed to communicate that she was having issues until the day the report was due and we did not have enough time to re-run `limma`. 

### Almas Khan
Originally, it was intended that Almas, being the one with the most experience with methylation data, was going to 'QC' the data with the raw files, to see if the GEO beta-values were normalized appropriately, and to check certain metadata that would be reflected in the raw `.idat` files. In the end, we did not have the computational resources to work from the raw data and had to just work with the GEO data. It turned out, that Almas had the most stats background out of all of us, so she ended up doing the linear regression with `limma`, and the hierarchical clustering that generated the heatmap.

### Denitsa Vasileva
Denitsa generated the beta density plot, as well as the pathway analysis after linear regression. Given Denitsa's experience, it was apt that she was assigned gene set enrichment analysis. She attempted the gene set enrichment analysis, but due to issues with `ermineR`, which even Paul failed to resolve (was an issue with `ermineJ`, the program in which `ermineR` is an R wrapper for, so she was unable to continue with the functional gene set enrichment analysis, and simply did pathway analysis only.

----

> Do you think the task assignments were fairly assigned and appropriate given each member’s background and skills? If no, how would you change it?

As you can see from our original `project_proposal.md` and our `progress_report.md`, the division of labour got shuffled around and changed quite a bit, based on course load at the time of the deliverable, as well as the skills each person had. We tried our best to reassign the tasks as fit to each of our skillsets as soon as we realized something was out of the scope of a group member. If I were to change anything, I would like to have at least someone in the statistics program in each group. Given our backgrounds, I do not think we had a strong enough grasp of certain statistical concepts to make certain decisions (e.g. additive vs interactive model). Tasks were sort of unfairly assigned at the beginning (where not all members were involved in every step of the project), but we slowly got the hang of things and redistributed the tasks in a way that made sense, where each member had a part in each deliverable. Every member contributed to each group-level deliverable.

To see detailed contributions of each member, we have tracked it in this [issue](https://github.com/STAT540-UBC/Repo_team_The-Splice-Girls_W2020/issues/37). We updated the main `README` of our repository to reflect the changes in division of labour.

## My Specific Contributions and Comments

### My Contributions

> Explain what are your contributions to the project.

I performed the PCA analysis to find an underlying signal, the beta-value distribution boxplots, the strip plot to visualize beta values of our differentially methylated genes, and the chromosome plot (most of the content covered in `seminar06` and `seminar07`). This was within my scope given my experience with R, and lack of background experience in statistics. 

On top of these statistical analyses, I took on most of the administrative tasks. I created GitHub issues to track our progress, updated the `README` files in each directory, organized the GitHub repository, submitted the group-level deliverables. I also initiated the use of GitHub branches for us to work collaboratively with minimal merge conflicts. 

I also ran some scripts on the GSC cluster that were too much for our laptops to handle.

> What worked well and what did not? What was the most challenging or rewarding moment during your group project?

Using GitHub branches and pull requests worked well for such a large collaborative project. The GitHub filesize limit was very annoying as our data could not be pushed to the repository, meaning that everyone had to run the script that generated the data. Most of our work was done in Rmarkdown due to the ease of displaying the plots inline in a 'notebook'-like environment, but near the end of the project, we learned that some of the Rmarkdowns took a _long_ time to knit, and converted all those Rmarkdown files to Rscripts. The most challenging part of the group project had to be dealing with such a large dataset, and the uncertainty of our model. We were uncertain about our model, and was constantly re-running it, which would take a while due to the data size, and then it would affect all downstream analysis, which had to be re-done each time we switched models. RStudio would freeze a lot, and some Rmarkdowns would knit overnight and still remain unfinished. The use of the GSC cluster was a lifesaver, but not all students would've had access to such a computational resource. The most rewarding moment during the project would be when an Rscript would finish running on the cluster, especially the `limma` results!

> How did you find having members of different backgrounds for a scientific project?

Almas, Denitsa, and I are actually in the same cohort in the Bioinformatics program. This posed a challenge that I did not see coming. We were in the same program so it was very easy to contact one another, especially in person. However, the fact that we are in the same cohort, also means we have the same course load, and were always busy at the exact same times. If we had different courseloads, perhaps another group member could have picked up the slack when the Bioinformatics students were busy. In that sense, the three of us have the same background, but having Nairuz from Medical Genetics gave us a different perspective when we were discussing `limma` models.

> What have you learned from this group project?

I learned it is not always a good idea to work with friends/cohort members on group projects, and that communication is key! (Lots of avoidable mishaps from miscommunication.) It is also a good idea to start analyses early as you never know how long analyses may take with such large datasets that are used in high-dimensional biology. Additionally, I learned that during initial dataset selection, many of the QC steps (PCA analysis, etc) should be conducted so as to not select a dataset with little-to-no signal or too many confounding factors, and that the metadata should be thoroughly checked to ensure that you have a balanced dataset, and a good enough sample size to have statistical power. All in all, I learned many statistical methods, how to use them, and what to use them for, etc. as this was my first statistics project!

> Any other comments on how the group project could have been structured differently (i.e. delivery requirements, assessment“)

I think it could be a good idea to tag a release instead of copying the SHA number for deliverable submissions. Additionally, I think it would be a good idea to make use of the GSC's [ORCA](https://github.com/bcgsc/orca) cluster to give high-performance computing resources to students. It is normally set up in term 1 for MICB405, and I think it could be adapted for this class. I also think that as (new) grad students, it is very easy to put off this group project, as the deliverables seem very far apart. It would be nice to mention the upcoming deliverables in class, and remind us, or give us more minor deliverables to keep us on track. My group definitely fell behind at some point, and just didn't prioritize this project due to other course work and thesis work. Finally, I think students could benefit from some pre-approved datasets (e.g. datasets that worked well in past projects or seminars). Students could be given the option to pick from a pre-approved list or choose their own.

## Question From Presentation

Question 1 from issue [#38](https://github.com/STAT540-UBC/Repo_team_The-Splice-Girls_W2020/issues/38):

> Since you have found that the first PCs are not related to the cancer_type, did you check if the first PCs are related to other covariates in data? If not, do you think if removing the first PCs from the data can improve the results?

I did generate the PC1 vs PC2 plots for all the other covariates, in which I did not observe any clustering (see [`PCA.md`](https://github.com/STAT540-UBC/Repo_team_The-Splice-Girls_W2020/blob/master/docs/PCA.md) from the group repository). However, to answer this question, _I re-did my PCA_ using the `prcomp(center = TRUE, and scale = TRUE)` instead of doing manual scaling to rule that out as a variable for error, for PC1 to PC10, every combination and for every covariate, with 540 plots total (45 per covariate). I generated one cohesive plot for each of the 12 covariates, which you can view [here](PCA.md) in full. 

I believe that removing the first PCs from the data can improve the results, especially if PC1 is heavily skewed by outliers (e.g. confounding factors/batch effects), or if the first three PCs contribute almost equally to each other. As Sara had said as well, another reason for not seeing any clustering by `cancer_type`:

> It may not be surprising that the samples don’t cluster by cancer_type based on the entire probeset if there are only a smaller subset of probes that are differentially methylated.

If I could have gone back, we could have subsampled our data to include only probes commonly used for cancer studies. That would've resolved challenges of working with a large dataset as well.

Additionally, in the paper that Sara suggested I read, about [surrogate variable analysis](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0030161), it is possible that there are unaccounted variables in our data, such as batch effects, technical variation, etc. As stated in the paper, there are are variables that may explain the variance in PC1 but are not associated with any given gene's methylation variation, such as 'unmodeled factors' and 'unmeasured factors' and 'gene-specific noise'. Given this information, I have tried plotting other PCs in addition to the first three (which were already plotted).

Since in our presentation, we found that the first PCs are not related to `cancer_type`, here are the first 10 PCs for `cancer_type`. I have generated the same plot for other covariates, which you can see in [PCA.md](PCA.md). Looking at the plot below, it seems that within the first 10 PCs, that it is PC7 that accounts for `cancer_type`:

![](images/PCA_cancer_type.png)

Checking `summary(pcs)`, PC7 only accounts for 0.638% of the variance:

```
Importance of components:
                          PC1     PC2     PC3     PC4    PC5    PC6     PC7     PC8     PC9    PC10
Standard deviation     7.2355 1.08636 0.87153 0.71411 0.6489 0.62943 0.62370 0.55663 0.49530 0.48852
Proportion of Variance 0.8582 0.01935 0.01245 0.00836 0.0069 0.00649 0.00638 0.00508 0.00402 0.00391
Cumulative Proportion  0.8582 0.87757 0.89002 0.89838 0.9053 0.91178 0.91816 0.92324 0.92726 0.93117
```

Looking through the plots for the other covariates in [`PCA.md`](PCA.md), it does indeed look like PC1/PC2/PC3 may be confounding factors/batch effects, as they do not seem to associate with any of the other covariates as well.

**Conclusion:** PC1 vs PC2 were checked for all covariates, where no clustering was observed. However, removing the first few PCs from the data _did in fact_ help, as I did find the PC associating with `cancer_type` in the remaining PCs (i.e. PC7!), as it seems the first few PCs may be associated with sources of variations that were unmodeled or unmeasured, such as confounding factors, technical variation, etc., skewing the first few PCs.
