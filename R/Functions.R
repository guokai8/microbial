#' @title calculat the richness for the phyloseq object
#' @importFrom phyloseq estimate_richness
#' @importFrom vegan rarefy
#' @importFrom vegan diversity
#' @importFrom vegan specnumber
#' @importFrom phyloseq otu_table
#' @param physeq phyloseq object
#' @param method A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "Observed","Chao1","ACE","Richness", "Fisher", "Simpson", "Shannon", "Evenness","InvSimpson".
#' @examples
#'  \dontrun{
#' data("GlobalPatterns",package="phyloseq")
#' physeq <- GlobalPatterns
#' rich <-richness(physeq,method=c("Simpson", "Shannon"))
#' }
#' @return data.frame of alpha diversity
#' @export
#' @author Kai Guo
richness<-function(physeq,method=c("Simpson", "Shannon")){
    method<-as.character(sapply(method,function(x)simpleCap(x),simplify = T))
    method<- match.arg(method,c("Observed","Chao1","ACE","Richness", "Fisher", "Simpson", "Shannon", "Evenness","InvSimpson"), several.ok = TRUE)
    df <- estimate_richness(physeq)
    if(!isTRUE(taxa_are_rows(physeq))){
        tab<-t(otu_table(physeq))
    }else{
        tab<-otu_table(physeq)
    }
    if("Evenness"%in%method){
        H<-diversity(tab)
        S <- specnumber(tab)
        J <- H/log(S)
        df_J<-data.frame(Evenness=J)
        df<-cbind(df,df_J)
    }
    if("Richness"%in%method){
        R<- rarefy(tab,min(rowSums(tab)))
        df_R<-data.frame(Richness=R)
        df<-cbind(df,df_R)
    }
    df<-df[,method]
    return(df)
}
#' @title calcaute beta diversity
#' @importFrom phyloseq ordinate
#' @param physeq A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method. A character string specifying ordination method. All methods available to the \code{ordinate} function
#'        of \code{phyloseq} are acceptable here as well.
#' @param distance. A string character specifying dissimilarity index to be used in calculating pairwise distances (Default index is "bray".).
#'                       "unifrac","wunifrac","manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower",
#'                       "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @author Kai Guo
.betadiv<-function(physeq,distance="bray",method="PCoA"){
    beta<-ordinate(physeq,method = method,distance = distance)
    pcs<-beta$values[,2]
    df<-beta$vectors
    return(list(beta=df,PCs=pcs))
}
#' @title PERMANOVA test for phyloseq
#' @importFrom phyloseq distance
#' @importFrom vegan adonis
#' @importFrom tidyr gather
#' @importFrom dplyr group_by
#' @importFrom dplyr do
#' @importFrom magrittr %>%
#' @param physeq A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param distance. A string character specifying dissimilarity index to be used in calculating pairwise distances (Default index is "bray".).
#'                       "unifrac","wunifrac","manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower",
#'                       "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param group (Required). Character string specifying name of a categorical variable that is preferred for grouping the information.
#'        information.
#' @examples
#'  \dontrun{
#' data("GlobalPatterns",package="phyloseq")
#' physeq <- GlobalPatterns
#' phy<-normalize(physeq)
#' beta <-beta_test(phy,group="SampleType")
#' }
#' @return PERMANOVA test result
#' @export
#' @author Kai Guo
beta_test<-function(physeq,group,distance="bray"){
    cat("Do PERMANOVA for: ",group,"\n")
    dist<-distance(physeq,method = distance)
    tab <- as(sample_data(physeq),"data.frame")
    tab<-tab[,group,drop=F]
    res<-NULL
    if(length(group)>1){
        res<- tab%>%gather(Group,val)%>%group_by(Group)%>%do(as.data.frame(adonis(dist~val,.)$aov.tab))
    }else{
        tab$Group <- tab[,group]
        res<-as.data.frame(adonis(dist~Group,tab)$aov.tab)
    }
    return(as.data.frame(res))
}

#' Normalize the phyloseq object with different methods
#' @importFrom phyloseq transform_sample_counts sample_data
#' @importFrom phyloseq taxa_are_rows nsamples otu_table
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors
#' @importFrom DESeq2 estimateDispersions varianceStabilizingTransformation
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom edgeR calcNormFactors
#' @param physeq A \code{phyloseq} object containing merged information of abundance,
#'      taxonomic assignment, sample data including the measured variables and categorical information
#'      of the samples, and / or phylogenetic tree if available.
#' @param method A list of character strings specifying \code{method} to be used to normalize the phyloseq object
#'      Available methods are: "relative","TMM","vst","log2".
#' @param group group (DESeq2). A character string specifying the name of a categorical variable containing  grouping information.
#' @examples
#'  \dontrun{
#' data("GlobalPatterns",package="phyloseq")
#' physeq <- GlobalPatterns
#' phy<-normalize(physeq)
#' }
#' @return phyloseq object with normalized data
#' @author Kai Guo
#' @export
normalize<-function(physeq,group,method="relative"){
    if(!taxa_are_rows(physeq)){
        physeq<-t(physeq)
    }
    otu<-as(otu_table(physeq),"matrix")
    tab<-as(sample_data(physeq),"data.frame")
    group<-tab[,group]
    if(method=="vst"){
        cat("Normalization using DESeq2 varianceStabilizingTransformation method \n")
        otu <- otu+1
        condition=group
        dds = DESeqDataSetFromMatrix(otu, DataFrame(condition), ~ condition)
        dds = estimateSizeFactors(dds)
        dds = estimateDispersions(dds)
        vst <- varianceStabilizingTransformation(dds)
        otu_table(physeq) <- otu_table(assay(vst), taxa_are_rows=TRUE)
    }
    if(method=="relative"){
        cat("Normalization using relative method \n")
        physeq<-transform_sample_counts(physeq,function(x)x/sum(x))
    }
    if(method=="TMM"){
        # modified from https://github.com/aametwally/MetaLonDA/blob/master/R/Normalization.
        cat("Normalization using TMM method \n")
        otu = otu + 1
        # Check `group` argument
        factors = calcNormFactors(otu, method="TMM")
        eff.lib.size = colSums(otu) * factors
        ref.lib.size = mean(eff.lib.size) #Use the mean of the effective library sizes as a reference library size
        count = sweep(otu, MARGIN = 2, eff.lib.size, "/") * ref.lib.size
        otu_table(physeq) <- otu_table(count, taxa_are_rows=TRUE)
    }
    if(method=="log2"){
        cat("Normalization using log2 of the RA method \n")
        physeq<-transform_sample_counts(physeq,function(x)log2(x/sum(x)+1))
    }
    return(physeq)
}
#' @title Calculate differential bacteria with DESeq2
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom phyloseq sample_data
#' @importFrom DESeq2 results DESeq
#' @param physeq A \code{phyloseq} object containing merged information of abundance,
#'      taxonomic assignment, sample data including the measured variables and categorical information
#'      of the samples, and / or phylogenetic tree if available.
#' @param group group (DESeq2). A character string specifying the name of a categorical variable containing  grouping information.
#' @param pvalue pvalue threshold for significant  results
#' @param padj adjust p value threshold for significant  results
#' @param log2FC log2 Fold Change threshold
#' @param gm_mean TRUE/FALSE calculate geometric means prior to estimate size factors
#' @param fitType either "parametric", "local", or "mean" for the type of fitting of dispersions to the mean intensity.
#' @examples
#'  \dontrun{
#' data("GlobalPatterns",package="phyloseq")
#' require(phyloseq)
#' samdf<-as(sample_data(physeq),"data.frame")
#' samdf$group<-c(rep("A",14),rep("B",12))
#' sample_data(physeq)<-samdf
#' res <- diff_test(physeq,group="group")
#' }
#' @return datafame with differential test with DESeq2
#' @author Kai Guo
#' @export
#'
diff_test<-function(physeq,group,pvalue=0.05,padj=NULL,log2FC=0,gm_mean=FALSE,fitType="local"){
    if(!taxa_are_rows(physeq)){
        physeq<-t(physeq)
    }
    otu <- as(otu_table(physeq),"matrix")
    tax <- as.data.frame(as.matrix(tax_table(physeq)))
    colData<-as(sample_data(physeq),"data.frame")
    colData$condition<-colData[,group]
    if(isTRUE(gm_mean)){
        countData<-round(otu, digits = 0)
    }else{
        countData<-round(otu, digits = 0)+1
    }
    dds <- DESeqDataSetFromMatrix(countData, colData, as.formula(~condition))
    if(isTRUE(gm_mean)){
        geoMeans = apply(counts(dds), 1, gm_mean)
        dds = estimateSizeFactors(dds, geoMeans = geoMeans)
    }
    dds = DESeq(dds, fitType=fitType)
    res = results(dds, cooksCutoff = FALSE)
    res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]))
    if(!is.null(padj)){
        pval<-padj
        sig <- rownames(subset(res,padj<pval&abs(log2FoldChange)>log2FC))
    }else{
        pval<-pvalue
        sig <- rownames(subset(res,pvalue<pval&abs(log2FoldChange)>log2FC))
    }
    res_tax$Significant<- "No"
    res_tax$Significant <- ifelse(rownames(res_tax) %in% sig, "Yes", "No")
    res_tax <- cbind(res_tax, tax[rownames(res),])
    return(as.data.frame(res_tax))
}


