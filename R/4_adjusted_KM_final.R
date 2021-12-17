#' Calculate the adjusted survival probability
#'
#' Use input data, time, status, grouping variables, adjusted covariates,
#' whether to use stratified model, and reference group as inputs
#'
#' @param data the input dataset
#' @param time column name of time variable
#' @param status column name of event status
#' @param group grouping variable
#' @param covlist list of covariates that should be included in the model
#' @param stratified_cox "Yes" refers to use stratified model, "No" refers to use Coxph
#' @param reference_group NULL- unstratified cox when stratified = No; "G&B"- G&B when stratified = Yes; Otherwise, Storer's approach will be performed when using a self-defined reference
#'
#' @return Output is a dataframe with adjusted survival probabilities. If the PH assumption is invalid or if practitioners need a method by which the event time points of the adjusted function match those of the unadjusted function, the stratified model should be used (Gail and Byar and Storer et. al), otherwise, unstratified FG model can be used.
#' @export
#'
#' @examples
#' # Data preparation
#' library(KMsurv)
#' data(bmt)
#' bmt$arm <- bmt$group
#' bmt$arm = factor(as.character(bmt$arm), levels = c("2", "1", "3"))
#' bmt$z3 = as.character(bmt$z3)
#' bmt$t2 = bmt$t2 * 12/365.25
#'
#' # Cox model
#' result1 = adjusted_KM(data = bmt, time = "t2", status = "d3", group = "arm", covlist = c("z1",
#' "z3"), stratified_cox = "No", reference_group = NULL)
#'
#' # Stratified cox: Gail&Byar's approach
#' result2 = adjusted_KM(data = bmt, time = "t2", status = "d3", group = "arm", covlist = c("z1",
#' "z3"), stratified_cox = "Yes", reference_group = "G&B")
#'
#' # Stratified cox: Storer's approach
#' result3 = adjusted_KM(data = bmt, time = "t2", status = "d3", group = "arm", covlist = c("z1",
#' "z3"), stratified_cox = "Yes", reference_group = "arm:2")
#'

adjusted_KM = function(data,time,status,group,covlist,stratified_cox="Yes",reference_group="G&B"){
  res = .adjusted_km(data = data,time=time,status=status,group=group,covlist=covlist,stratified_cox=stratified_cox,reference_group=reference_group)
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


