#' use R code to calculate baseline CIF
#'
#' Use input data, time, status,grouping variables, adjusted covariates,
#' events of interests, the known beta as inputs
#'
#' @noRd
#'
#' @param data input data
#' @param time time
#' @param status status
#' @param group grouping variable
#' @param covlist adjusted covariates
#' @param event_code event of interests
#' @param beta beta from FG or stratified model
#'
#' @return a dataframe
#' @keywords internal
#'
baseline_hazard_cif = function(data,time,status,group,covlist,event_code,beta){

  EvTime0 = sort(data[[time]][data[[status]]==event_code])
  EvTimes = unique(EvTime0)
  nn = nrow(data);mm = ncol(data)
  inputvars = data.frame(data[,c(group,covlist)])
  inputvars_temp = inputvars # save for names
  ################ Need to confirm with Dr. Kim???? ##########################
  suppressWarnings(
    for(i in 1:ncol(inputvars)){
      if(class(inputvars[,i])=="numeric"|class(inputvars[,i])=="integer"){ # if numeric, keep
        inputvars[,i] = inputvars[,i]
      }else if(class(inputvars[,i])=="factor"|class(inputvars[,i])=="character"){ # if character/factor, change to factor
        inputvars[,i] = as.factor(inputvars[,i])
        if(is.na(sum(as.numeric(levels(as.factor(inputvars[,i])))))){ # if string factor change to numeric, and minus 1
          inputvars[,i] = as.character(as.numeric(as.factor(inputvars[,i]))-1)
        }else{
          inputvars[,i] = as.character(as.numeric(as.factor(inputvars[,i])))
        }
      }
    }
  )



  Zs = lapply(1:ncol(inputvars),function(x){
    if(class(inputvars[,x])=="character"|class(inputvars[,x])=="factor"){
      inputvars[,x] = factor(inputvars[,x])
      if(length(levels(inputvars[,x]))==1){
        inputvars[,x] = as.character(inputvars[,x])
        t =  data.frame(inputvars[,x])
      }else{
        inputvars[,x] = as.character(inputvars[,x])
        t =  as.data.frame(model.matrix(~inputvars[,x],data = inputvars))
        t = data.frame(t[,-1])
      }
      names(t) = paste0(names(inputvars_temp)[x],levels(inputvars_temp[,x])[-1])
    }else{
      t = as.data.frame(inputvars[,x])
      t = t-mean(inputvars[,x],na.rm=T)
      names(t) = names(inputvars)[x]
    }
    t
  })

  Zs = do.call(cbind,Zs)




  # Compute the KM estimate of censoring events
  CenTime = c(0,unique(sort(data[[time]][data[[status]]==0]))) # Unique censored time
  KM_est = matrix(0,length(CenTime),2)

  for (ij in 1:length(CenTime)){
    nj = sum(data[[time]] >= CenTime[ij])
    dj = sum((data[[time]] == CenTime[ij]) & (data[[status]]==0))

    if(ij == 1){
      KM_est[1,] = cbind(CenTime[ij],(nj-dj)/nj)
    }else{
      KM_est[ij,] = cbind(CenTime[ij],KM_est[ij-1,2]*((nj-dj)/nj))
    }
  }

  # Compute G(ti) and G(Xj)
  Gti = matrix(0,length(EvTimes),2)
  for (ij in 1:length(EvTimes)){
    Gtii = KM_est[KM_est[,1] <= EvTimes[ij], ]
    Gti[ij,] = cbind(EvTimes[ij],Gtii[length(Gtii)])
  }

  Gxi = matrix(0,nn,2)
  for (ij in 1:nn){
    Gxii = KM_est[KM_est[,1] <= data[[time]][ij], ]
    Gxi[ij,] = cbind(data[[time]][ij],Gxii[length(Gxii)])
  }

  # Compute the weight matrix
  WeightW = lapply(1:nn,function(ii){
    w_jj = c()
    for(jj in 1:length(EvTimes)){
      if(data[[time]][ii]<EvTimes[jj]){
        G1 = Gti[Gti[,1]==EvTimes[jj],2]
        G2 = unique(Gxi[Gxi[,1]==data[[time]][ii],2])

        if(data[[status]][ii]==1){ ########## Need to change
          w_jj = c(w_jj,G1/G2)
        }else if(data[[status]][ii]==0){
          w_jj = c(w_jj,0)
        }else{
          w_jj = c(w_jj, G1/G2)
        }
      }else{
        w_jj = c(w_jj,1)
      }
    }
    w_jj
  })

  WeightY = lapply(1:nn,function(ii){
    y_jj = c()
    for(jj in 1:length(EvTimes)){
      if(data[[time]][ii]<EvTimes[jj]){
        G1 = Gti[Gti[,1]==EvTimes[jj],2]
        G2 = unique(Gxi[Gxi[,1]==data[[time]][ii],2])

        if(data[[status]][ii]==1){ ########## Need to change
          y_jj = c(y_jj,0)
        }else if(data[[status]][ii]==0){
          y_jj = c(y_jj,1)
        }else{
          y_jj=c(y_jj,1)
        }
      }else{
        y_jj = c(y_jj,1)
      }
    }
    y_jj
  })

  WeightY = data.frame(do.call(rbind,WeightY))
  WeightW = data.frame(do.call(rbind,WeightW))

  beta = beta

  for(i in 1:ncol(Zs)){
    Zs[,i] = as.numeric(as.character(Zs[,i]))
  }

  Zs = as.matrix(Zs)
  CumHaz0 =  matrix(NA,ncol=2,nrow=length(EvTimes))

  for(jj in 1:length(EvTimes)){
    NumTies = sum(EvTime0 == EvTimes[jj])
    W1 = WeightW[,jj]*WeightY[,jj]
    S02 = sum(W1*exp(Zs%*%beta))/nn
    if(jj==1){
      CumHaz0[1,]= cbind(EvTimes[1],NumTies/S02)
    }else{
      CumHaz0[jj,] = cbind(EvTimes[jj],CumHaz0[jj-1,2]+NumTies/S02)
    }
  }

  CumHaz0[,2]= CumHaz0[,2]/nn;
  CumHaz0 = data.frame(rbind(0,CumHaz0)) # baseline hazard
  Z0 = matrix(0,ncol = ncol(Zs))
  CIF0 = cbind(CumHaz0[,1], 1-exp(-CumHaz0[,2]%*%exp((Z0)%*%beta))) # baseline CIF = 1 - exp(-H0)
  return(CumHaz0)
}



