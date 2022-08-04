#' Calculate the covariate-adjusted survival probability based on Cox/stratified Cox model
#'
#' The function computes the covariate-adjusted KM functions based on Cox or stratified Cox model. More details is presented in
#' 'AdjKM.CIF: An R package for estimating the covariate-adjusted Kaplan-Meier and cumulative incidence functions'. Three approaches are
#' implemented. See the value section.
#'
#'
#' @param data the input dataset
#' @param time column name of time variable
#' @param status column name of event status
#' @param group grouping variable
#' @param covlist list of covariates that should be included in the Cox model
#' @param stratified_cox "Yes" refers to use stratified model, "No" refers to use Coxph
#' @param reference_group NULL- unstratified cox when stratified = No; "G&B"- G&B when stratified = Yes; Otherwise, Storer's approach will be performed when using a self-defined reference
#'
#' @return Output is a dataframe with adjusted survival probabilities.
#' The stratified Cox regression model can be applied if the argument "stratified_cox="Yes"" is specified, where the group variable serves as the stratification variable in the model. The reference group per the Gail and Byar method includes all subjects in each group by specifying "reference_group = "G&B"", while the Storer method selects a specific group as the reference group by "reference_group = "variable_name:level""
#' If "stratified_cox = "No"", and "reference_group = NULL", then the unstratified Cox is used to estimate the adjusted survival probability.
#' Unstratified Cox model is appropriate if the PH assumption is valid. As the common baseline survival function is used, there's no difference in "event time points" between groups, and the number of events does not match the actual number of events in each group.
#' If need a method by which the event time points of the adjusted function match those of the unadjusted function, the stratified Cox can be selected as an alternative.
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


