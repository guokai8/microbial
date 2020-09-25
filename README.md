# microbial
## microbial <a href="https://travis-ci.org/guokai8/microbial"><img src="https://travis-ci.org/guokai8/microbial.svg" alt="Build status"></a>  [![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)  [![](https://img.shields.io/badge/devel%20version-0.0.3-green.svg)](https://github.com/guokai8/microbial)  ![Code Size:](https://img.shields.io/github/languages/code-size/guokai8/microbial)
An R package for microbial community analysis with phyloseq

This package is developed to enhance the available statistical analysis procedures in R by providing simple functions to analysis and visulazing the 16S rRNA data. 

## Installation
```
library(devtools)
install_github("guokai8/microbial")
``` 
## Quick tour
```{r} 
library(microbial)
```   
## Functions
```
# calcuate the alpha diversity 
richness(phyloseq,method=c("Simpson", "Shannon"))
# plot alpha diversity
plotalpha(phyloseq,method=c("Simpson", "Shannon"))
# make barplot for relative abundance
phy <- normalize(phyloseq)
plotbar(phy,level="Phylum")
# plot beta diversity(PCoA)
plotbeta(phy,distance="bray",method="PCoA")
# perform PERMANOVA test
beta_test(phy,group,distance="bray")
# do differential analysis with DESeq2
res <- diff_test(phyloseq,group)
# plot the differential results
plotdiff(res,level="Genus")
# do some test
?do_ttest
?do_wilcox
?do_aov
```
## Note
The _microbial_ package was bulit based on the phyloseq. The package is still under development. New functions will be provided soon.

## Contact information

For any questions please contact guokai8@gmail.com
