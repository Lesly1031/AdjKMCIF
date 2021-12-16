#' use R code to calculate the adjusted CIF
#'
#' Use input data, time, status,grouping variables, adjusted covariates,
#' events of interests, whether to use stratified model, and defining reference group as inputs
#'
#' @param data the input dataset
#' @param time column name of time variable
#' @param status column name of event status
#' @param group grouping variable
#' @param covlist list of covariates that are included in the model
#' @param event_code event of interests
#' @param stratified "Yes" refers to use stratified model, "No" refers to use FG
#' @param reference_group NULL- unstratified FG when stratified = No; "G&B"- G&B when stratified = Yes;
#'
#' @return a dataframe
#' @export
#'
#' @examples
#'
#' # unstratified FG
#' result1 = adjusted_CIF(data = bmt, time = "t2", status = "CenCI", group = "arm",
#' covlist = c("z1", "z3"), event_code = 1, stratified = "No", reference_group = NULL)
#' # stratified G$B
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


