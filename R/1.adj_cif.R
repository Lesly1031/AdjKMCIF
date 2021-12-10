#' use R code to calculate adjusted CIF
#'
#' Use input data, time, status,grouping variables, adjusted covariates,
#' events of interests, whether to use stratified model, and defining reference group as inputs
#'
#' @param data input data
#' @param time time
#' @param status status
#' @param group grouping variable
#' @param covlist adjusted covariates
#' @param event_code event of interests
#' @param stratified "Yes" refers to use stratified model, "No" refers to use FG
#' @param reference_group NULL- unstratified FG when stratified = No; "G&B"- G&B when stratified = Yes;
#' @return a dataframe
#' @keywords internal
#'
#'
adj_cif = function(data,time,status,group,covlist,event_code,stratified="Yes",reference_group=NULL){

##########################################################################
# This function is to get the adjusted CIF
# data: input data
# time, status
# group: grouping variable
# covlist: adjusted covariates
# event_code: Event of Interest
# beta: estimate of Fine-Gray regression
# reference_group: NULL- unstratified FG when stratified = No; "G&B"- G&B when stratified = Yes; self- defined results
##########################################################################
  ############ Delete this part once complete #########
  # data = dt
  # time = "RFS"
  # status = "CI_RFS"
  # group = "Arm"
  # covlist = c("Age","Severity")
  # event_code =1

  # data = bmt
  # time = "t2"
  # status = "CenCI"
  # group = "arm"
  # covlist = c("z1","z3")
  # event_code =1
####################################################
  grouplist = group
  covlist = covlist
  # number of group variable (number of strata)
  numgroup = length(levels(factor(data[[grouplist]])))
  # Number of covariates to adjust for
  numcov = ncol(data.frame(data[,covlist]))
  # Number of obs
  numobs = nrow(data)
  # create input data
  indata = cbind.data.frame(data[,time],data[,status],data[,grouplist],data[,covlist])
  names(indata) = c("time","status","strata",paste0("z",1:length(covlist)))
  indata$time = as.numeric(as.character(indata$time))
  indata$status = as.numeric(as.character(indata$status))

################### Get beta estimate #####################
  inputvars  <- paste0("z",1:numcov)
  indata_cox1 = indata
  if(stratified=="No"){
    if(is.null(reference_group)){

      suppressWarnings(
        for(i in 3:ncol(indata_cox1)){
          if(class(indata_cox1[,i])=="numeric"|class(indata_cox1[,i])=="integer"){ # if numeric, keep
            indata_cox1[,i] = indata_cox1[,i]
          }else if(class(indata_cox1[,i])=="factor"|class(indata_cox1[,i])=="character"){ # if character/factor, change to factor
            indata_cox1[,i] = as.factor(indata_cox1[,i])
            if(is.na(sum(as.numeric(levels(as.factor(indata_cox1[,i])))))){ # if string factor change to numeric, and minus 1
              indata_cox1[,i] = as.character(as.numeric(as.factor(indata_cox1[,i]))-1)
            }else{
              indata_cox1[,i] = as.character(as.numeric(as.factor(indata_cox1[,i])))
            }
          }
        }
      )

      fit_formula = as.formula(paste("prodlim::Hist(time,status)~strata+", paste(inputvars, collapse=" + ")))
      fit = riskRegression::FGR(fit_formula,data=indata,cause = event_code)
      beta = fit$crrFit$coef
      beta = matrix(beta,ncol=1)
    }else{
      stop("No reference group for non-stratified cox model")
    }
  }else if (stratified=="Yes"){
    if(!is.null(reference_group)){
      suppressWarnings(
      for(i in 3:ncol(indata_cox1)){
        if(class(indata_cox1[,i])=="numeric"|class(indata_cox1[,i])=="integer"){ # if numeric, keep
          indata_cox1[,i] = indata_cox1[,i]
        }else if(class(indata_cox1[,i])=="factor"|class(indata_cox1[,i])=="character"){ # if character/factor, change to factor
          indata_cox1[,i] = as.factor(indata_cox1[,i])
          if(is.na(sum(as.numeric(levels(as.factor(indata_cox1[,i])))))){ # if string factor change to numeric, and minus 1
            indata_cox1[,i] = as.character(as.numeric(as.factor(indata_cox1[,i]))-1)
          }else{
            indata_cox1[,i] = as.character(as.numeric(as.factor(indata_cox1[,i])))

          }
        }
      }
      )
      dtt = data.frame(na.omit(indata_cox1))
      cov_list = dtt[,-c(1:3)]
      fit = crrSC::crrs(ftime = dtt$time,fstatus = dtt$status,cov1 =cov_list,strata = dtt$strata ,failcode = event_code)
      beta = fit$coef
      beta = matrix(beta,ncol=1)

    }else{
      stop("Please choose 'G&B', or self-define a reference group for stratified cox model")
    }
  }

################ Baseline CIF ###############
  if(stratified=="No"){
    if(is.null(reference_group)){
      CIF0 = baseline_hazard_cif(data,time,status,group,covlist,event_code=event_code,beta)
      CIF0 = data.frame(CIF0)
      names(CIF0) = c("time","CIF0")
    }else{
      stop("No reference group for non-stratified cox model")
    }
  }else if(stratified=="Yes"){
    if(!is.null(reference_group)){
      beta_stratify = c(0,beta)
      beta_stratify = matrix(beta_stratify,ncol=1)

      cov_list = names(indata)[-c(1:3)]
      data_sub = lapply(1:length(unique(indata$strata)),function(x){
        dt_sub = subset(indata,indata$strata==unique(indata$strata)[[x]])
      })
      CIF0 = lapply(1:length(data_sub),function(x){
        data_sub[[x]]$strata = as.character(data_sub[[x]]$strata)
        res = baseline_hazard_cif(data_sub[[x]],"time","status","strata",cov_list,event_code=event_code,beta_stratify)
        res$strata = unique(indata$strata)[x]
        res
      })
      CIF0 = do.call(rbind,CIF0)
      names(CIF0) = c("time","CIF0","strata")
    }else{
      stop("Please choose 'G&B', or self-define a reference group for stratified cox model")
    }
  }

############### Linear predictor ###################
  if(stratified=="No"){
    if(is.null(reference_group)){
      indata_sort = indata[order(indata$time),]
      indata_temp_raw = data.frame(indata_sort[,-c(1:2)])
      new_matrix = lapply(1:ncol(indata_temp_raw),function(x){
        if(class(indata_temp_raw[,x])=="character"|class(indata_temp_raw[,x])=="factor"){
          t =  as.data.frame(model.matrix(~indata_temp_raw[,x],data = indata_temp_raw))
          t = data.frame(t[,-1])
          indata_temp_raw[,x] = factor(indata_temp_raw[,x])
          names(t) = paste0(names(indata_temp_raw)[x],levels(indata_temp_raw[,x])[-1])
        }else{
          t = as.data.frame(indata_temp_raw[,x])
          t = t-mean(indata_temp_raw[,x],na.rm=T)
          names(t) = names(indata_temp_raw)[x]
        }
        t
      })

      new_matrix = do.call(cbind,new_matrix)
      coef_fit = beta
      linear_pred = as.matrix(new_matrix)%*%coef_fit
    }else{
      stop("No reference group for non-stratified cox model")
    }
  }else if(stratified=="Yes"){
    if(!is.null(reference_group)){
      indata_sort = indata[order(indata$time),]
      indata_temp_raw = data.frame(indata_sort[,-c(1:3)])
      new_matrix = lapply(1:ncol(indata_temp_raw),function(x){
        if(class(indata_temp_raw[,x])=="character"|class(indata_temp_raw[,x])=="factor"){
          t =  as.data.frame(model.matrix(~indata_temp_raw[,x],data = indata_temp_raw))
          t = data.frame(t[,-1])
          indata_temp_raw[,x] = factor(indata_temp_raw[,x])
          names(t) = paste0(names(indata_temp_raw)[x],levels(indata_temp_raw[,x])[-1])
        }else{
          t = as.data.frame(indata_temp_raw[,x])
          t = t-mean(indata_temp_raw[,x],na.rm=T)
          names(t) = names(indata_temp_raw)[x]
        }
        t
      })
      new_matrix = do.call(cbind,new_matrix)
      coef_fit = beta
      linear_pred = as.matrix(new_matrix)%*%coef_fit

    }else{
      stop("Please choose 'G&B', or self-define a reference group for stratified cox model")
    }
  }

  coxout = data.frame(cbind(indata_sort,linear_pred))
  names(coxout)[ncol(coxout)]="linear_predictor"
  coxout$time = as.numeric(as.character(coxout$time))
  coxout$status = as.numeric(as.character(coxout$status))
  coxout$strata = factor(coxout$strata)
  coxout$time = signif(coxout$time,8)
  #############################################################################
  #---------------- CIF probability for each group ------------------------
  #############################################################################
  coxout$strata = factor(coxout$strata)
  strata_num = length(levels(coxout$strata))

 cif_probability = lapply(1:strata_num,function(x){
    cif_prob(data=data,coxout=coxout,base_res=CIF0,event_code=event_code,strata_group = levels(coxout$strata)[[x]],stratified=stratified,reference_group = reference_group)
  })
 ######################################################################
 #----------------survival probability-------------------------
 ######################################################################
 cif_prob_cb = do.call(cbind,cif_probability)
 coxout_event = subset(coxout,coxout$status==event_code)
 unique_eventTime = data.frame(signif(unique(coxout_event$time),8));names(unique_eventTime)="time"
 res = data.frame(cbind(unique_eventTime[order(unique_eventTime$time,decreasing = F),],cif_prob_cb))

 colnames(res) = c("time",names(table(data[,grouplist])))

 res

}


