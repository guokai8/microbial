#‘ calculat the richness for the phyloseq object
#' @importFrom phyloseq estimate_richness
#' @importFrom vegan rarefy
#' @importFrom vegan diversity
#' @importFrom vegan specnumber
#' @importFrom phyloseq otu_table
#' @param physeq phyloseq object
#' @param method A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "Observed","Chao1","ACE","Richness", "Fisher", "Simpson", "Shannon", "Evenness","InvSimpson".
#' @return data.frame of alpha diversity
#' @export
#' @author Kai Guo
richness<-function(physeq,method){
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
#' calcaute beta diversity
#' @importFrom phyloseq ordinate
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
#' PERMANOVA for phyloseq
#' @importFrom phyloseq distance
#' @importFrom vegan adonis
#' @importFrom dplyr gather
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
#' @export
#' @author Kai Guo
beta.test<-function(physeq,group,distance){
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

#' normalize the phyloseq object
#'

