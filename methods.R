library(caret)
library(caretEnsemble)
library(rcdk)
library(RSNNS)
library(ModelMetrics)
library(assertthat)

maxn <- function(x, n) {
  partial <- length(x) - n + 1
  x[x >= sort(x, partial = partial)[partial]]
}

r2 <- function(t, p){
  rss <- sum((p - t) ^ 2)
  tss <- sum((t - mean(t)) ^ 2)
  return(1 - rss/tss)
}

applyDescriptors <- function(smiles, descriptors=c(), remove.sparsity=F){
  df <- as.data.frame(
    do.call(
      cbind, 
      lapply(descriptors, function(descriptor){
        message(paste0('Applying descriptor ', parent.frame()$i))
        res <- descriptor(smiles)
        return(res)
      })
    )  
  )
  rownames(df) <- rownames(smiles)
  if(remove.sparsity) df <- df[,-caret::nearZeroVar(df)]
  return(df)
}

# folds.group.cat <- function(y, p, k){
#   # @ param y <vector> The vector of response variables
#   # @ param p <double> If y is a factor - the probability of sampling each level in y for each fold
#   # @ param k <integer> The number of folds
#   p.hasnames <- !is.null(names(p))
#   if(length(p) == nlevels(y)){
#     if(p.hasnames){
#       p100 <- p
#     }else{
#       names(p) <- levels(y)
#     }
#   }else if(length(p) < nlevels(y)){
#     nm.missing <- levels(y) %in% names(p)
#     p.missing <- rep((1-sum(p))/length(nm.missing), times=length(nm.missing))
#     p100 <- c(p, p.missing)
#     if(p.hasnames){
#       names(p100) <- c(levels(y), nm.missing)
#     }else{
#       names(p.missing) <- levels(y)
#     }
#   }else{
#     warning()
#     if(p.hasnames){
#       
#     }else{
#       
#     }
#   }
#   y.idx.bylevel <- sapply(levels(y), function(l){
#     lapply(1:length(y), function(i){
#       return(y[i] == l)
#     })
#   })
#   lapply(1:k, function(f){
#     
#   })
# }
# folds.group.cat(as.factor(DF.binary.split$targetsTrain), )

prepare.dataset <- function(dataset, ratio, response.col='Activity', desc.set=NULL, method='random', remove.sparsity=T){
  # TODO:: 
  # - How should the dataset/resamples be split? random, bins, X, Y, KNN, SOM? 
  # # = Splitting will be implemented by folds.group when its fixed.
  
  df <- as.data.frame(dataset)
  
  if(method == 'random'){
    
    df <- df[sample(1:nrow(df), nrow(df)), sample(1:ncol(df), ncol(df))]  
  }
  
  # Select the response column
  y <- df[,response.col]
  
  # Calculate descriptors and split train/test
  IAtom <- parse.smiles(as.vector(df$Canonical.Smiles))
  X <- applyDescriptors(IAtom, desc.set, remove.sparsity)
  
  df.split <- splitForTrainingAndTest(X,y, ratio=ratio)
  
  # Sanity checks for train/test dimensions
  all(
    assert_that(
      nrow(df.split$inputsTrain) == length(df.split$targetsTrain)  
    ) &&
      assert_that(
        nrow(df.split$inputsTest) == length(df.split$targetsTest)  
      ) &&
      assert_that(
        ncol(df.split$inputsTrain) == ncol(df.split$inputsTest)  
      )
  )
  return(df.split)
}

train.ensemble <- function(X.train, y.train, X.test, y.test, ensemble.args, stack.args){
  ensemble <- do.call(caretList, ensemble.args)
  stack <- caretStack(ensemble, method=stack.args$method, metric=stack.args$metric, trControl=stack.args$trControl)
  return(stack)
}

performance.ensemble <- function(model){
  # predictions <- predict(model.stack, newdata=X.ic50.test)
  # (rmse <- RMSE(predictions, y.ic50.test))
  # (rsq <- R2(predictions, y.ic50.test))
  # plt <- plot(predictions, y.ic50.test, ylim = c(-1, 8), xlim=c(-1,8))
  # abline(0,1, col='blue')
  # 
  # resamps <- resamples(
  #   ensemble
  # )
  # ensemble$
  #   summary(resamps)
  # modelCor(resamps)
  # 
  # theme1 <- trellis.par.get()
  # theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
  # theme1$plot.symbol$pch = 16
  # theme1$plot.line$col = rgb(1, 0, 0, .7)
  # theme1$plot.line$lwd <- 2
  # trellis.par.set(theme1)
  # bwplot(resamps, layout = c(3, 1))
}

pred.models <- function(X, models){
  out <- data.frame(row.names = rownames(X))
  for(idx in 1:length(models)){
    m <- models[[idx]]
    n <- names(models)[idx]
    p <- predict(m, data.frame(X))
    if(is.null(n)) n <- idx
    out[,n] <- p
  }
  return(out)
}

screenCompunds <- function(smis, descriptors, remove.sparsity, models = list(activity=NULL, effect=NULL, affinity=NULL)){
  # for smi in smis 
  #   sample = extracted features for smi using descriptors
  #   estimate sample activity
  #   if sample is active
  #     predict sample effect class (i.e. agonist/antagonist)
  #     predict sample effect continuous (i.e. IC50/EC50)
  #     predict sample affinity continuous (i.e. Ki/Kb)
  #   end
  # end
  # sort active compounds in their respective classes (i.e. agonist/antagonist)
  # choose the top N compounds for 3D docking
}


