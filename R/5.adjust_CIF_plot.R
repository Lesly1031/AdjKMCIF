#' Generate adjusted CIF plot without CI
#'
#' Use results from adjusted_CIF() and input dataset to generate adjusted CIF plot (users can also produce a figure by using ggplot2)
#'
#' @param res results from adjusted_CIF()
#' @param data the input dataset
#'
#' @return Adjusted CIF plot will be shown after running this function
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
#'# Adjusted CIF plot without CI
#' result = adjusted_CIF(data = bmt, time = "t2", status = "CenCI", group = "arm",
#' covlist = c("z1", "z3"), event_code = 1, stratified = "Yes", reference_group = "arm:2")
#' adjCIF_plot(result, data = bmt)
#'

adjCIF_plot = function(res,data){
  res_long = res

  p = ggplot2::ggplot(res_long)+
    ggplot2::geom_step(ggplot2::aes(x=time,y = prob, group =class,linetype=class,color= class),size=1.5)+
    ggplot2::theme_classic()+
    ggplot2:: ylim(c(0,1))
  return(p)
}

#' Generate adjusted CIF plot with bootstrap CI
#'
#' Use results from boot_ci_adj_cif() and input dataset to generate adjusted CIF plot with bootstrap CI (users can also produce a figure by using ggplot2)
#'
#' @param res results from adjusted_CIF()
#' @param data the input dataset
#'
#' @return Adjusted CIF plot with CI will be shown after running this function
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
#' # Adjusted CIF plot with bootstrap CI
#' result1_1 = boot_ci_adj_cif(boot_n = 100, ci_cut = c(0.025, 0.975), data = bmt, time = "t2",
#' status = "CenCI", group = "arm", covlist = c("z1", "z3"), event_code = 1, "No",
#'                             NULL)
#' adjCIF_CI_plot(result, data = bmt)
#'
#'
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


