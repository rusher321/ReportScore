
# ReportScore

<!-- badges: start -->
<!-- badges: end -->

The goal of ReportScore is to compute the reportscore of pathway or module

## Installation

This package can be installed using [devtools](http://cran.r-project.org/web/packages/devtools/index.html).

```r
devtools::install_github('rusher321/ReportScore@master')
```

## Example

This is a basic example which shows you how to solve a common problem:

```{r example}
library(ReportScore)

## basic example code
pr <- read.table("testdata/all.KEGG.abun.txt",row.names = 1,header = T,sep = "\t")
grp <- read.table("testdata/sampleinfo.txt",row.names = 1,header = T,sep = "\t")

grp_sub <- grp[which(grp$TimePoint == "Day2"), 1, drop=F]
grp_sub$Protein <- ifelse(grp_sub$Protein == "C-Pork", "B-Beef", grp_sub$Protein)
pr_sub <- pr[, rownames(grp_sub)]

res_D2 <- ReporterScore(pr_sub, grp_sub, paired = F, database = "./database", occ = 0.1)

fig <- ReportVis(res_D2, color = c("#2470a0", "#DE3C3C"), exclude = T)

## if you want update the datebase 
download_data(db_dir = "database/")

## if you want plot the output
fig <- ReportVis(res, color = c("#2470a0", "#DE3C3C"), exclude = T)

```
## Figure 

the ***test.R*** have real data can to test,the output like 

<img src="fig/test.png" width="500">

We welcome comments, criticisms, and especially contributions! GitHub
issues are the preferred way to report bugs, ask questions, or request
new features. You can submit issues here:

<https://github.com/rusher321/ReportScore/issues>

ToDo 
-----

- [√]Visualization of reportscore result 
- [ ]parallel compute of multi group 


Meta
----

-   Please [report any issues or
    bugs](https://github.com/rusher321/microbiotaPair/issues).
-   License: MIT
-   Get citation information for `ReportScore` in R doing
    `citation(package = 'ReportScore')`
-   Please note that this project is released with a [Contributor Code
    of Conduct](CONDUCT.md). By participating in this project you agree
    to abide by its terms.
