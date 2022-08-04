#' Bootstrap CI for adjusted CIF
#'
#' Computes the covariate-adjusted CIF functions along with confidence intervals by the bootstrap percentile method.
#'
#' @param boot_n bootstrap sample size
#' @param ci_cut default c(0.025, 0.975) bootstrap 95\% CI
#' @param data the input dataset
#' @param time column name of time variable
#' @param status column name of event status
#' @param group grouping variable
#' @param covlist list of covariates that should be included in the model
#' @param event_code event of interests
#' @param stratified "Yes" refers to use stratified model, "No" refers to use Fine and Gray regression
#' @param reference_group NULL- unstratified FG when stratified = No; "G&B"- G&B when stratified = Yes; Otherwise, Storer's approach will be performed when using a self-defined reference
#'
#' @return Output is a dataframe with average number of adjusted CIF probabilities, as well as 2.5\% and 97.5\% percentiles.
#' @export
#'
#' @examples
#' library(KMsurv)
#' data(bmt)
#' bmt$arm <- bmt$group
#' bmt$arm = factor(as.character(bmt$arm), levels = c("2", "1", "3"))
#' bmt$z3 = as.character(bmt$z3)
#'
#' bmt$CenCI <- 0
#' for (ii in 1:137) {
#'   if (bmt$d3[ii] == 0) {
#'     bmt$CenCI[ii] <- 0
#'   } else {
#'     if (bmt$d2[ii] == 1) {
#'       bmt$CenCI[ii] <- 1
#'     } else {
#'       bmt$CenCI[ii] <- 2
#'     }
#'   }
#' }
#'
#' bmt$t2 = bmt$t2 * 12/365.25
#'
#' # unstratified Fine-Gray regression model
#' result1_1 = boot_ci_adj_cif(boot_n = 100, ci_cut = c(0.025, 0.975), data = bmt, time = "t2",
#'  status = "CenCI", group = "arm", covlist = c("z1", "z3"), event_code = 1, "No",
#'                             NULL)
#'
#' # stratified FG: Gail&Byar's approach
#' result1_2 = boot_ci_adj_cif(boot_n = 100, ci_cut = c(0.025, 0.975), data = bmt, time = "t2",
#'    status = "CenCI", group = "arm", covlist = c("z1", "z3"), event_code = 1, "Yes",
#'    "G&B")
#'
#'  # stratified Storer
#'  result1_3 = boot_ci_adj_cif(boot_n = 100, ci_cut = c(0.025, 0.975), data = bmt, time = "t2",
#'     status = "CenCI", group = "arm", covlist = c("z1", "z3"), event_code = 1, "Yes",
#'    "arm:2")
#'


boot_ci_adj_cif <- function(boot_n=100,ci_cut = c(0.025,0.975),data,time,status,group,covlist,event_code,stratified,reference_group){

  ########### Get bootstrap samples ################
  resample = lapply(1:boot_n,function(i){
    set.seed(i)
    index = sample(1:nrow(data),replace=T)
    dt_sample = data.frame(data[index,])
    dt_sample
  })

  ########## Calculate adjusted survival probability on each bootstrap sample####
  boot_adj_cif = function(boot_time = boot_n){
    adj_cif_prob = lapply(1:boot_time,function(x){
      t = .adj_cif(data = resample[[x]],time,status,group,covlist,event_code =event_code,stratified=stratified,reference_group=reference_group)
      #list(data.frame(cbind(t[,1],t[,2])),data.frame(cbind(t[,1],t[,3])))
      t1 = vector(mode="list",length = ncol(t)-1)
      for(i in 1:ncol(t)-1) {
        t1[i]= list(data.frame(cbind(t[,1],t[,i+1])))
      }
      t1
    })
  }
  #system.time(boot_adj_cif())
  t =boot_adj_cif(boot_n)

  boot_ci = lapply(1:length(t[[1]]),function(x){

    t_sub = sapply(t,"[",x)
    t_sub = suppressWarnings(Reduce(function(x,y) merge(x,y,all=TRUE,by="X1"),t_sub))
    t_sub_mean = rowMeans(t_sub[,-1],na.rm=T)
    t_sub_quantile = sapply(1:nrow(t_sub),function(t){
      quantile(t_sub[t,-1],na.rm=T,ci_cut)
    })
    t_sub_quantile = data.frame(cbind(t_sub[,1],t(t_sub_quantile)))
    t_sub_quantile = data.frame(cbind(t_sub_quantile$V1,t_sub_mean,t_sub_quantile[,-1]))
    names(t_sub_quantile)[1]="time"
    t_sub_quantile$time = signif(t_sub_quantile$time,8)
    new_time = data.frame("time"=sort(data[[time]]))
    new_time$time = signif(new_time$time,8)
    t_sub_quantile = suppressWarnings(merge(t_sub_quantile,new_time,by="time",all.y=T))

    for(i in 1:ncol(t_sub_quantile)){
      t_sub_quantile[,i] = zoo::na.locf(t_sub_quantile[,i],na.rm=F)
      t_sub_quantile[,i] = zoo::na.locf(t_sub_quantile[,i],fromLast=T)
    }
    t_sub_quantile = data.frame(t_sub_quantile[which(!duplicated(t_sub_quantile$time)),])
    names(t_sub_quantile)[c(ncol(t_sub_quantile)-1,ncol(t_sub_quantile))]=c("lower","upper")
    t_sub_quantile$class = rep(levels(data[[group]])[x],nrow(t_sub_quantile))
    t_sub_quantile


  })
  return(boot_ci)
}
