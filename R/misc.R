simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}

#
lightcolor<-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175','#e6194b', '#3cb44b', '#ffe119', '#4363d8','#f58231', '#911eb4',
                '#46f0f0', '#f032e6', '#bcf60c',
                '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8',
                '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075',
                '#808080'
)
#
distcolor<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8',
         '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c',
         '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8',
         '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075',
         '#808080', '#ffffff', '#000000')
#' do anova test and return results as data.frame
#' @importFrom rstatix anova_test
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @export
#' @author Kai Guo
do_aov<-function(x,group,...){
    d<-x[,setdiff(colnames(x),group)]
    d$group<-x[,group]
    d<-d%>%gather(type,val,-group)
    res<-d%>%group_by(type)%>%anova_test(val~group,...)
    return(res)
}

#' do t.test
#' @importFrom rstatix t_test
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @export
#' @author Kai Guo
do_ttest<-function(x,group,ref=NULL,...){
    d<-x[,setdiff(colnames(x),group)]
    d$group<-x[,group]
    d<-d%>%gather(type,val,-group)
    res<-d%>%group_by(type)%>%t_test(val~group,ref.group = ref,...)
    return(res)
}

#' do wilcox test
#' @importFrom rstatix wilcox_test
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @export
#' @author Kai Guo
do_wilcox<-function(x,group,ref=NULL,...){
    d<-x[,setdiff(colnames(x),group)]
    d$group<-x[,group]
    d<-d%>%gather(type,val,-group)
    res<-d%>%group_by(type)%>%wilcox_test(val~group,ref.group = ref,...)
    return(res)
}
