#' Covariate-adjusted KM plot
#'
#' Use the results from adjusted_KM() and input dataset to generate adjusted survival plot. The figure is compatible with the results from adjusted_KM()
#'
#' @param res results from adjusted_KM()
#' @param data the input dataset
#'
#' @return Adjusted KM plot will be shown
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
#' adjKM_plot(result1, data = bmt)
#'



adjKM_plot = function(res,data){
    res_long = res

    p = ggplot2::ggplot(res_long)+
      ggplot2::geom_step(ggplot2::aes_string(x="time",y = "prob", group ="class",linetype="class",color= "class"),size=1.5)+
      ggplot2::theme_classic()+
      ggplot2:: ylim(c(0,1))
    return(p)
}


#' Generate adjusted KM plot with including bootstrap CI
#'
#' Use the results from boot_ci_adj_km() and input dataset to generate adjusted KM plot with bootstrap CI
#'
#' @param res results from boot_ci_adj_km()
#' @param data the input dataset
#'
#' @return Adjusted CIF plot with CI will be shown after running this function
#' @export
#'
#' @examples
#' install.packages("KMsurv")
#' library(KMsurv)
#' data(bmt)
#' bmt$arm <- bmt$group
#' bmt$arm = factor(as.character(bmt$arm), levels = c("2", "1", "3"))
#' bmt$z3 = as.character(bmt$z3)
#' bmt$t2 = bmt$t2 * 12/365.25
#'
#' result1_1 = boot_ci_adj_km(boot_n = 100, ci_cut = c(0.025, 0.975), data = bmt, time = "t2",
#' status = "d3", group = "arm", covlist = c("z1", "z3"), stratified_cox = "No",
#' reference_group = NULL)
#'
#' adjKM_CI_plot(result1_1, bmt)

adjKM_CI_plot = function(res,data){
  boot_ci = data.table::rbindlist(res)
  names(boot_ci)[2]="prob"
  p =
    ggplot2::ggplot(data.frame(boot_ci))+
    ggplot2::geom_step(ggplot2::aes_string(x="time",y = "prob", group ="class",linetype="class",color= "class"),size=1.2)+
    ggplot2::geom_ribbon(ggplot2::aes_string(x="time",y = "prob", group ="class",linetype="class",color= "class",ymin="lower",ymax="upper",fill="class"),alpha=0.3)+
    ggplot2::ylim(c(0,1))+
    theme_classic()
  return(p)

}


