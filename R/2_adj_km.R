#' use R code to calculate adjusted KM
#'
#' Internal function for adj KM. Use input data, time, status,grouping variables, adjusted covariates,
#' whether to use stratified model, and defining reference group as inputs
#'
#' @noRd
#'
#' @param data input data
#' @param time time
#' @param status status
#' @param group grouping variable
#' @param covlist adjusted covariates
#' @param stratified "Yes" refers to use stratified model, "No" refers to use Cox
#' @param reference_group NULL- unstratified cox when stratified = No; "G&B"- G&B when stratified = Yes;
#' @return a dataframe
#'
#'

.adjusted_km <- function(data,time,status,group,covlist,stratified_cox="Yes",reference_group="G&B"){

  # levels of group variable
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


  # save for generate coxout
  indata_cox = cbind.data.frame(data[,time],data[,status],data[,grouplist],data[,covlist])
  names(indata_cox)[1:3] = c("time","status","strata")
  indata_cox = data.frame(indata_cox,data[[grouplist]])
  names(indata_cox)[ncol(indata_cox)]=grouplist
  indata_cox1 = indata ## Save for later (basehaz survival survival_3.2-11)

########################## Unstratified KM ##########################
  if(stratified_cox=="No"){
    if(is.null(reference_group)){
      ####### get beta estimates #####
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
      inputvars  <- paste0("z",1:numcov)

      for (i in 3:ncol(indata_cox1)){
        if(class(indata_cox1[,i])=="numeric"|class(indata_cox1[,i])=="integer"){
          indata_cox1[,i] = indata_cox1[,i] - mean(indata_cox1[,i],na.rm=T)
        }else{
          indata_cox1[,i] = indata_cox1[,i]
        }
      }

      fit_formula = as.formula(paste("survival::Surv(time,status)~", paste0(paste("strata +",inputvars, collapse=" + "))))
      fit = survival::coxph(fit_formula,data=indata)

      #### baseline hazard ###
      fit1 = survival::coxph(fit_formula,data=indata_cox1)

      base_res = survival::basehaz(fit1,centered = F) ########### Issue with survival package survival_3.2-11  !!!!!!
      base_res1 = survival::basehaz(fit,centered = F)
      base_res$strata = base_res1$strata

      base_res$basesurv = exp(-base_res$hazard)

      ### Linear predictor ###
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

      #linear_mean = colMeans(indata_cox2[,-c(1:2)],na.rm=TRUE)
      coef_fit = as.matrix(fit$coefficients,ncol=1)
      linear_pred = as.matrix(new_matrix)%*%coef_fit
    }else{
      stop("No reference group for non-stratified cox model")
    }
  }else if(stratified_cox=="Yes"){
    if(!is.null(reference_group)){
      ### get beta estimates ###
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


      inputvars  <- paste0("z",1:numcov)

      for (i in 3:ncol(indata_cox1)){
        if(class(indata_cox1[,i])=="numeric"|class(indata_cox1[,i])=="integer"){
          indata_cox1[,i] = indata_cox1[,i] - mean(indata_cox1[,i],na.rm=T)
        }else{
          indata_cox1[,i] = indata_cox1[,i]
        }
      }

      fit_formula = as.formula(paste("survival::Surv(time,status)~", paste0(paste(inputvars, collapse=" + "),"+ strata(strata)")))
      fit = survival::coxph(fit_formula,data=indata)
      ###### baseline hazard ######
      fit1 = survival::coxph(fit_formula,data=indata_cox1)

      base_res = survival::basehaz(fit1,centered = F) ########### Issue with survival package survival_3.2-11  !!!!!!
      base_res1 = survival::basehaz(fit,centered = F)
      base_res$strata = base_res1$strata

      base_res$basesurv = exp(-base_res$hazard)

      ###### linear predictor #############
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

      #linear_mean = colMeans(indata_cox2[,-c(1:2)],na.rm=TRUE)
      coef_fit = as.matrix(fit$coefficients,ncol=1)
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

  #---------------- survival probability for each group ------------------------
  coxout$strata = factor(coxout$strata)
  strata_num = length(levels(coxout$strata))

  survival_prob = lapply(1:strata_num,function(x){
    .survival_prob(data=data,coxout=coxout,base_res=base_res,strata_group = levels(coxout$strata)[[x]],stratified_cox=stratified_cox,reference_group = reference_group)
  })

  #----------------survival probability-------------------------
  surv_prob_cb = do.call(cbind,survival_prob)
  coxout_event = subset(coxout,coxout$status==1)
  unique_eventTime = data.frame(signif(unique(coxout_event$time),8));names(unique_eventTime)="time"
  res = data.frame(cbind(unique_eventTime[order(unique_eventTime$time,decreasing = F),],surv_prob_cb))

  colnames(res) = c("time",names(table(data[,grouplist])))

  res
}
