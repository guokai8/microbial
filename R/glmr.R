#' @title Do the generalized linear model regression
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom broom tidy
#' @param physeq phyloseq object
#' @param group the group factor to regression
#' @param factors  a vector to indicate adjuested factors
#' @param family binomial() or gaussian()
#' @export
#' @author Kai Guo

glmr<-function(physeq,group,factors=NULL,family=binomial(link = "logit")){
      if (!taxa_are_rows(physeq)) {
            physeq <- t(physeq)
      }
      otu <- as(otu_table(physeq), "matrix")
      otu <- as.data.frame(t(otu))
      colnames(otu)<-paste0('ASV_',colnames(otu))
      samx <- as(sample_data(physeq)[,c(group,factors)],"matrix")
      samd <- as.data.frame(samx)
      dd <- cbind(samd[rownames(otu),],otu)
      dd[,group]<-as.factor(dd[,group])
      cat('##########################################\n')
      cat('Do the generalized linear model regression with ',factors,'adjusted',"\n")
      cat(paste0(group,"~",paste(factors,"x",sep="+")),"\n")
      cat('##########################################\n')
      rr<-lapply(names(otu),function(x)tidy(glm(as.formula(paste0(group,"~",paste(factors,x,sep="+"))),data=dd,family=family)))
      names(rr)<- sub('ASV_','',names(otu))
      res <- do.call(rbind,rr)
      res <- res[grep('ASV_',res$term),]
      res$term<-sub('ASV_','',res$term)
      res$padj <- p.adjust(res$p.value,method="BH")
      return(res)
}
