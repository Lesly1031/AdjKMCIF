#' Calculate the adjusted CIF probability
#'
#' Use input data, time, status, grouping variables, adjusted covariates,
#' events of interests, whether to use stratified model, and reference group as inputs
#'
#' @param data the input dataset
#' @param time column name of time variable
#' @param status column name of event status
#' @param group grouping variable
#' @param covlist list of covariates that should be included in the model
#' @param event_code event of interests
#' @param stratified "Yes" refers to use stratified model, "No" refers to use Fine and Gray regression
#' @param reference_group NULL- unstratified FG when stratified = No; "G&B"- G&B when stratified = Yes; Otherwise, Storer's approach will be performed when using a self-defined reference
#'
#' @return Output is a dataframe with adjusted CIF probabilities. If the PH assumption is invalid or if practitioners need a method by which the event time points of the adjusted function match those of the unadjusted function, the stratified model should be used (Gail and Byar and Storer et. al), otherwise, unstratified FG model can be used.
#' @export
#'
#' @examples
#'
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
#' # unstratified FG
#' result1 = adjusted_CIF(data = bmt, time = "t2", status = "CenCI", group = "arm",
#' covlist = c("z1", "z3"), event_code = 1, stratified = "No", reference_group = NULL)
#' # stratified G&B
#' result2 = adjusted_CIF(data = bmt, time = "t2", status = "CenCI", group = "arm",
#' covlist = c("z1", "z3"), event_code = 1, stratified = "Yes", reference_group = "G&B")
#' # stratified Storer
#' result3 = adjusted_CIF(data = bmt, time = "t2", status = "CenCI", group = "arm",
#' covlist = c("z1", "z3"), event_code = 1, stratified = "Yes", reference_group = "arm:2")

adjusted_CIF = function(data,time,status,group,covlist,event_code =1,stratified="Yes",reference_group="G&B"){
  res = .adj_cif(data = data,time=time,status=status,group=group,covlist=covlist,event_code=event_code,stratified=stratified,reference_group=reference_group)
  res$time = signif(res$time,8)
  new_time = data.frame("time"=sort(data[[time]]))
  new_time$time = signif(new_time$time,8)
  res = merge(res,new_time,by="time",all.y=T)

  for(i in 1:ncol(res)){
    res[,i] = zoo::na.locf(res[,i],na.rm=F)
    res[,i] = zoo::na.locf(res[,i],fromLast=T)
  }
  res = data.frame(res[which(!duplicated(res$time)),])
  colnames(res) = c("time",names(table(data[,group])))
  # Generate plot
  grouplist = group
  numgroup = length(levels(factor(data[[group]])))

  res_long = tidyr::gather(res,class,prob,names(table(data[,group]))[1]:names(table(data[,group]))[numgroup])
  return(res_long)
}


