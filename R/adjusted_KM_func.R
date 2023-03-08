#' Covariate-adjusted survival probability based on Cox/stratified Cox model
#'
#' This function calculates the adjusted Kaplan-Meier functions by taking into account the covariates, using either Cox or stratified Cox model. You can find further information about this in the publication titled "AdjKM.CIF: An R package for estimating the covariate-adjusted Kaplan-Meier and cumulative incidence functions". The function offers three different approaches, which are discussed in the value section.
#'
#'
#' @param data the input dataset
#' @param time column name of time variable
#' @param status column names of event status, 0 is the censored and 1 is the event
#' @param group group variable being compared
#' @param covlist list of covariates included in the Cox model
#' @param stratified_cox "Yes" refers to use stratified model, "No" refers to use Coxph
#' @param reference_group "NULL" - No reference group required for the Cox PH model;"G&B" - the Gail and Byar method; "group:level" - the reference group for the Storer method (e.g., "arm:2" in the BMT data)

#'
#' @return The output of this function is a dataframe that contains the adjusted survival probabilities. If you specify "stratified_cox="Yes"", the function will apply a stratified Cox regression model, where the group variable serves as the stratification variable. You can choose the reference group using the Gail and Byar method by specifying "reference_group = "G&B"", which includes all subjects in each group, or using the Storer method by specifying "reference_group = "variable_name:level"", which selects a specific group as the reference group.

#' If you set "stratified_cox = "No"" and "reference_group = NULL", the function will use the unstratified Cox model to estimate the adjusted survival probability. This is appropriate if the PH assumption is valid, as it uses a common baseline survival function and there is no difference in event time points between groups. However, the number of events may not match the actual number of events in each group.

#' If you need to ensure that the event time points of the adjusted function match those of the unadjusted function, you can choose the stratified Cox model as an alternative.
#'
#' @export
#'
#' @examples
#' # Data preparation
#' install.packages("KMsurv")
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

  res_long = tidyr::gather(res,"class","prob",names(table(data[,group]))[1]:names(table(data[,group]))[numgroup])
  return(res_long)
}


