#' plot beta diversity
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 xlab ylab stat_ellipse
#' @importFrom ggplot2 theme_light
#' @importFrom phyloseq taxa_are_rows
#' @param physeq A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method. A character string specifying ordination method. All methods available to the \code{ordinate} function
#'        of \code{phyloseq} are acceptable here as well.
#' @param distance. A string character specifying dissimilarity index to be used in calculating pairwise distances (Default index is "bray".).
#'                       "unifrac","wunifrac","manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower",
#'                       "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param group (Required). Character string specifying name of a categorical variable that is preferred for grouping the information.
#'        information.
#' @param color user defined color for group
#' @param size the point size
#' @param ellipse draw ellipse or not
plotbeta<-function(physeq,group,distance="bray",method="PCoA",color=NULL,size=3,ellipse=FALSE){
    if(!taxa_are_rows(physeq)){
        physeq <- t(physeq)
    }
    beta<-.betadiv(physeq,distance = distance,method=method)
    df <- as.data.frame(beta$beta)
    PCs <- beta$PCs
    tab <- as(sample_data(physeq),"data.frame")
    df$group<-tab[,group]
    if(is.null(color)){
        color<-lightcolor[1:length(unique(df$group))]
    }
    p <- ggplot(df,aes(Axis.1,Axis.2,color=group))+geom_point(size=size)+scale_color_manual(values=color)
    p <- p+theme_light(base_size=15)+xlab(paste0("Axis1 (",round(PCs[1]*100,2),"%)"))+ylab(paste0("Axis2 (",round(PCs[2]*100,2),"%)"))
    if(isTRUE(ellipse)){
        p <- p + stat_ellipse()
    }
    p
}

#' plot alpha diversity
#' @importFrom rstatix t_test
#' @importFrom rstatix wilcox_test
#' @importFrom ggpubr ggboxplot
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom ggpubr facet
#' @importFrom ggplot2 xlab ylab scale_color_manual theme
#' @param physeq A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param group group (Required). A character string specifying the name of a categorical variable containing  grouping information.
#' @param method A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "Observed","Chao1","ACE","Richness", "Fisher", "Simpson", "Shannon", "Evenness","InvSimpson".
#' @param pvalue pvalue threshold for significant dispersion results
#' @param padj adjust p value threshold for significant dispersion results
#' @param wilcox use wilcoxon test or not
#' @return Returns a ggplot object. This can further be manipulated as preferred by user.
#' @author Kai Guo
#' @export
plotalpha<-function(physeq,group,method=c("Simpson", "Shannon"),
                    pvalue=0.05,padj=NULL,wilcox=FALSE){
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    rich<-richness(physeq,method = method)
    name<-levels(factor(colnames(rich)))
    tab<-as(sample_data(physeq),"data.frame")
    rich$group<-tab[rownames(rich),group]
    if(isTRUE(wilcox)){
        res<-do_wilcox(rich,"group")
    }else{
        res<-do_ttest(rich,"group")
    }
    res <- subset(res,p < pvalue)
    if(!is.null(padj)){
        res <- subset(res,p.adj < padj)
    }
    vals<-rich%>%gather(type,val,-group)%>%group_by(type,group)%>%summarise(ma=max(val))%>%spread(group,ma)
    pos <- apply(res, 1, function(x)max(vals[vals$type==x[1],x[3:4]]))
    mpos <- apply(res, 1, function(x)min(vals[vals$type==x[1],x[3:4]]))
    p<-rich%>%gather(type,val,-group)%>%ggboxplot(x="group",y="val",color="group")
    p<-facet(p,facet.by = "type",scales = "free_y",ncol = length(method))+
        stat_pvalue_manual(res,label = "p",y.position = pos+mpos/5)+
        xlab("")+ylab("")+
        theme(legend.position = "none",axis.text.x=element_text(angle=90,vjust=0.5, hjust=1))+
        scale_color_manual(values=lightcolor[unique(rich$group)])
    p
}


#' plot bar
#' @importFrom phyloseq psmelt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 element_blank scale_y_continuous
#' @importFrom dplyr group_by_at
#' @importFrom dplyr vars pull
#' @importFrom dplyr one_of
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom rlang `!!`
#' @param physeq A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param level the level to plot
#' @param color A vector of character use specifying the color
#' @param top the number of most abundance bacteria to display
#' @param fontsize.x the size of x axis label
#' @param fontsize.y the size of y axis label
#' @return Returns a ggplot object. This can further be manipulated as preferred by user.
#' @author Kai Guo
#' @export
plotbar<-function(physeq,level="Phylum",color=NULL,top=5,fontsize.x = 5, fontsize.y = 12){
    pm <- psmelt(physeq)
    if(is.null(color)){
        len<-length(unique(pm[,level]))
        color<-lightcolor[1:len]
    }
    group_var<-c("Sample",level)
    d<-pm%>%group_by_at(vars(one_of(group_var)))%>%summarise(su=sum(Abundance))
    d <- as.data.frame(d)
    d[,level][is.na(d[,level])]<-"NA"
    dx <- pm%>%group_by_at(vars(one_of(level)))%>%summarise(su=sum(Abundance))
    dx <- dx[order(dx$su,decreasing = T),]
    sel <- dx%>%head(top)%>%select(!!level)%>%pull(1)
    d <- d[d[,level]%in%sel,]
    p<-ggplot(d,aes_string("Sample","su",fill=level))+
    geom_bar(stat = "identity",position = "fill")+scale_fill_manual(values=color)+
    theme_light()+
    scale_y_continuous(expand = c(0, 0.001)) +
    theme(axis.text.x=element_text(angle=90,size=fontsize.x, vjust=0.5, hjust=1),
              axis.text.y=element_text(size=fontsize.y),
              panel.background = element_blank(),axis.ticks.x = element_blank())+
        xlab("")+ylab("")
    p
}

#' @title plot differential results
#' @importFrom ggplot2 ggplot theme geom_point element_text
#' @importFrom ggplot2 aes_string scale_color_manual theme_light
#' @param sigtab
#' @param level the level to plot
#' @param color A vector of character use specifying the color
#' @param pvalue pvalue threshold for significant  results
#' @param padj adjust p value threshold for significant  results
#' @param log2FC log2 Fold Change threshold
#' @param size size for the point
#' @param fontsize.x the size of x axis label
#' @param fontsize.y the size of y axis label
#' @param horiz horizontal or not (TRUE/FALSE)
#' @return ggplot object
#' @author Kai Guo
#' @export
plotdiff<-function(sigtab,level="Genus",color=NULL,pvalue=0.05,padj=NULL,log2FC=0,size=3,fontsize.x=5,fontsize.y=10,horiz=TRUE){
    if(!is.null(padj)){
        pval<-padj
        sigtab <- subset(res,padj<pval&abs(log2FoldChange)>log2FC)
    }else{
        pval<-pvalue
        sigtab <- subset(res,pvalue<pval&abs(log2FoldChange)>log2FC)
    }
    x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels=names(x))
    if(is.null(color)){
        len<-length(unique(sigtab$Phylum))
        color<-lightcolor[1:len]
    }
    # Genus order
    x <- tapply(sigtab$log2FoldChange, sigtab[,level], function(x) max(x))
    x <- sort(x, TRUE)
    sigtab[,level] <- factor(as.character(sigtab[,level]), levels=names(x))
    p <- ggplot(sigtab, aes_string(x=level, y="log2FoldChange", color="Phylum"))+
        geom_point(size=3) +theme_light()
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=fontsize.x),
              axis.text.y = element_text(size=fontsize.y))+
        scale_color_manual(values=color)
    if(isTRUE(horiz)){
        p<-p+coord_flip()
    }
    p
}




