# microbial
## microbial <a href="https://travis-ci.org/guokai8/microbial"><img src="https://travis-ci.org/guokai8/microbial.svg" alt="Build status"></a>  [![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)  [![](https://img.shields.io/badge/devel%20version-0.0.8-green.svg)](https://github.com/guokai8/microbial)  ![Code Size:](https://img.shields.io/github/languages/code-size/guokai8/microbial)
An R package for microbial community analysis with phyloseq

This package is developed to enhance the available statistical analysis procedures in R by providing simple functions to analysis and visulize the 16S rRNA data. 

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
data("GlobalPatterns",package="phyloseq")
physeq <- GlobalPatterns
richness(physeq,method=c("Simpson", "Shannon"))
# plot alpha diversity
plotalpha(physeq,method=c("Simpson", "Shannon"),group="SampleType")
# make barplot for relative abundance
phy <- normalize(physeq)
plotbar(phy,level="Phylum")
# plot beta diversity(PCoA)
plotbeta(phy,group="SampleType",distance="bray",method="PCoA")
# perform PERMANOVA test
beta_test(phy,group="SampleType",distance="bray")
# do differential analysis with DESeq2
require(phyloseq)
samdf<-as(sample_data(physeq),"data.frame")
samdf$group<-c(rep("A",14),rep("B",12))
sample_data(physeq)<-samdf
res <- diff_test(physeq,group="group")
# plot the differential results
plotdiff(res,level="Genus")
# do LEfse analysis
res <- ldamarker(physeq,group="group")
# plot LEfse results
plotLDA(res,group=c("A","B"),lda=5,pvalue=0.05)
# use RandomForest to select markers
res <- biomarker(physeq,group="group")
# do some test
?do_ttest
?do_wilcox
?do_aov
```
## Note
The _microbial_ package was bulit based on the phyloseq. The package is still under development. New functions will be provided soon.

## Contact information

For any questions please contact guokai8@gmail.com
