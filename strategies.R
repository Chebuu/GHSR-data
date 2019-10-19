library(dplyr)

strategy.thresh <- function(result = pred.models(), thresh = list(activity=function(act){act=='active'}, effect=function(eff){as.numeric(eff)<10000}, affinity=function(aff){as.numeric(aff)<10000})){
  if(all(names(thresh) %in% colnames(result))){
    iterseq <- names(thresh)
  }else{
    iterseq <- 1:length(thresh)
  }
  result$hit <- apply(result, 1, function(r){
    res <- all(sapply(iterseq, function(tn){
      thresh[[tn]](r[[tn]])
    }))
    return(res)
  }) 
  return(result)
}

strategy.sort <- function(result = pred.models(), sortby = c('activity', 'effect', 'affinity')){
  index.col <- function(df, col){
    return(df[,col])
  }
  return(
    result[
      with(result, do.call(order, lapply(sortby, index.col, df=result))),
    ]
  )
}

