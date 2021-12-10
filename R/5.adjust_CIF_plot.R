#' Generate adjusted CIF plot
#'
#' Use results from adjusted_CIF() and input dataset as inputs
#'
#' @param res results from adjusted_CIF()
#' @param data the input dataset
#'
#' @return a figure
#' @export
#'
#' @examples
#' result = adjusted_CIF(data = bmt, time = "t2", status = "CenCI", group = "arm",
#' covlist = c("z1", "z3"), event_code = 1, stratified = "Yes", reference_group = "arm:2")
#' adjCIF_plot(result, data = bmt)

adjCIF_plot = function(res,data){
  res_long = res

  p = ggplot2::ggplot(res_long)+
    ggplot2::geom_step(ggplot2::aes(x=time,y = prob, group =class,linetype=class,color= class),size=1.5)+
    ggplot2::theme_classic()+
    ggplot2:: ylim(c(0,1))
  return(p)
}


adjCIF_CI_plot = function(res,data){
  boot_ci = rbindlist(res)
  names(boot_ci)[2]="prob"
  p =
    ggplot2::ggplot(data.frame(boot_ci))+
    ggplot2::geom_step(ggplot2::aes(x=time,y = prob, group =class,linetype=class,color= class),size=1.2)+
    ggplot2::geom_ribbon(ggplot2::aes(x=time,y = prob, group =class,linetype=class,color= class,ymin=lower,ymax=upper,fill=class),alpha=0.3)+
    ggplot2::ylim(c(0,1))+
    theme_classic()
  return(p)
}


