#' Internal functions. Use R code to calculate adjusted CIF probability
#'
#' Use input data, time, status,grouping variables, adjusted covariates,
#' events of interests, whether to use stratified model, and defining reference group as inputs
#'
#' @noRd
#'
#' @return a dataframe
#'
#'
.cif_prob = function(data,coxout,base_res,event_code,strata_group,stratified,reference_group){
  ##############################################################
  # This function is to get the CIF probability
  # data: input data
  # coxout: combination of used variables and linear predictors
  # base_res: baseline hazard result
  # event_code
  # strata_group: strata for stratified FIne-gray regression
  # stratified: Yes- stratified FG, No - unstratified FG
  # reference_group: NULL- unstratified FG when stratified = No; "G&B"- G&B when stratified = Yes; self- defined results
  ###############################################################
  # make the time roughly the same for coxout and base_res for future merge
  coxout$time = as.numeric(as.character(coxout$time)) # coxout: combination of time, status, covariates and linear predictors
  coxout$status = as.numeric(as.character(coxout$status))
  coxout$strata = factor(coxout$strata)
  coxout$time = signif(coxout$time,8)
  base_res$time  = signif(base_res$time,8)

  coxout_event = subset(coxout,coxout$status==event_code)
  unique_eventTime = data.frame(signif(unique(coxout_event$time),8));names(unique_eventTime)="time"
  stratified=stratified
  if(stratified=="No"){
    if(is.null(reference_group)){
      size = sum(coxout$strata==strata_group)
      numtime = nrow(unique_eventTime) # length of unique event times

      coxout_sub = subset(coxout,coxout$strata==strata_group)
      coxout_sub = coxout_sub[order(coxout_sub$time,coxout_sub$status,decreasing=T),]

      # basehaz for each strata
      # base_res$strata =  factor(base_res$strata)
      # baseres_sub = subset(base_res,base_res$strata==strata_group)## change!! strata_group="TAC/MTX"
      baseres_sub=base_res
      baseres_sub = merge(unique_eventTime,baseres_sub,by="time",all.x=TRUE)
      baseres_sub = baseres_sub[order(baseres_sub$time),]
      # t$hazard= na.locf(t$hazard,na.rm=F)
      # t$hazard = na.locf(t$hazard,fromLast = T,na.rm=F)
      baseres_sub= as.numeric(baseres_sub$CIF0)



      adjcif = matrix(0,numtime,1)

      for(i in 1:size){
        expbz = exp(coxout_sub$linear_predictor[i])
        cumuhaz_res = baseres_sub
        cif = 1-exp(-cumuhaz_res)^expbz # estimate of survival function for group 1

        adjcif = adjcif + cif
      }
      adjcif = adjcif/size
    } else{
      stop("No reference group for non-stratified cox model")
    }
  }else if(stratified=="Yes"){
    if(is.null(reference_group)){
      stop("Please choose  'G&B', or self-define a reference group for stratified cox model")
    }

    else if(reference_group=="G&B"){ # Reference = "all", Gail and Byar approach

      size = sum(coxout$strata==strata_group)
      numtime = nrow(unique_eventTime) # length of unique event times

      coxout_sub = subset(coxout,coxout$strata==strata_group)
      coxout_sub = coxout_sub[order(coxout_sub$time,coxout_sub$status,decreasing=T),]

      # basehaz for each strata
      base_res$strata =  factor(base_res$strata)
      baseres_sub = subset(base_res,base_res$strata==strata_group)## change!! strata_group="TAC/MTX"

      baseres_sub = merge(unique_eventTime,baseres_sub,by="time",all.x=TRUE)
      baseres_sub = baseres_sub[order(baseres_sub$time),]
      # t$hazard= na.locf(t$hazard,na.rm=F)
      # t$hazard = na.locf(t$hazard,fromLast = T,na.rm=F)
      baseres_sub= as.numeric(baseres_sub$CIF0)



      adjcif = matrix(0,numtime,1)

      for(i in 1:size){
        expbz = exp(coxout_sub$linear_predictor[i])
        cumuhaz_res = baseres_sub
        cif = 1-exp(-cumuhaz_res)^expbz # estimate of survival function for group 1

        adjcif = adjcif + cif
      }
      adjcif = adjcif/size

      ################################################

    } else{ ## any user defined group
      reference_var = strsplit(reference_group,split=":")[[1]][1]
      reference_level = strsplit(reference_group,split=":")[[1]][2]
      max_num = nrow(subset(data,data[[reference_var]]==reference_level))
      numtime = nrow(unique_eventTime) # length of unique event times

      coxout_sub = subset(coxout,coxout[["strata"]]==reference_level)
      coxout_sub = coxout_sub[order(coxout_sub$time,coxout_sub$status,decreasing=T),]

      # basehaz for each strata
      base_res$strata =  factor(base_res$strata)
      baseres_sub = subset(base_res,base_res$strata==strata_group)## change!! strata_group="TAC/MTX"

      # baseres_sub = subset(base_res,base_res$strata==strata_group)## change!!
      baseres_sub = merge(unique_eventTime,baseres_sub,by="time",all.x=TRUE)
      baseres_sub = baseres_sub[order(baseres_sub$time),]
      # t$hazard= na.locf(t$hazard,na.rm=F)
      # t$hazard = na.locf(t$hazard,fromLast = T,na.rm=F)
      baseres_sub= as.numeric(baseres_sub$CIF0)

      adjcif = matrix(0,numtime,1)

      for(i in 1:max_num){
        expbz = exp(coxout_sub$linear_predictor[i])
        cumuhaz_res = baseres_sub
        cif = 1-exp(-cumuhaz_res)^expbz # estimate of survival function for group 1

        adjcif = adjcif + cif
      }
      adjcif = adjcif/max_num
    }
  }
  return(adjcif)
}

