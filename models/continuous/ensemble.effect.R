library(tidyverse)
library(doParallel)

source('methods.R')
source('descriptors.R')

USE_SAVED_DSET = T
SAVE_DSET = F

DO_PARALLEL = F
N_CORES = detectCores() - 1 

if(DO_PARALLEL){
  if(!exists('cl')){
    cl <- makeCluster(N_CORES)  
  }
  registerDoParallel(cl)
}else{
  if(exists('cl')) stopCluster(cl)
  registerDoSEQ()
}

##########
## Data ##
##########
if(USE_SAVED_DSET){
  # Load the saved dataset which has descriptors already applied and has already been split for train/test
  DF.eff.split <- readRDS('data/datasets/Activity.Continuous.Set3.rds')
}else{
  # Load the parsed dataset and prepare it
  DF.cont <- read.csv('data/parsed/GHSR.ACTIVITY.CONTINUOUS.csv', stringsAsFactors = F, row.names=1)  
  # Apply descriptors and split train/test
  DF.eff.split <- prepare.dataset(
    DF.cont, 
    response.col='Standard.Value',
    desc.set=DESC_SET, 
    method='random', 
    remove.sparsity=F
  )
  if(SAVE_DSET){
    saveRDS(DF.eff.split, 'data/datasets/Activity.Continuous.Set3.rds') 
  }
  
}

###########
## Model ##
###########
ensemble.tuneList = list(
  svmPoly = caretModelSpec('svmPoly', tuneGrid=expand.grid(degree=c(1,2,3,4,5,6,7,8,9,10), scale=seq(0.0001,0.01, length.out=10), C=c(0.8,3, length.out=10))), # degree = 3 scale = 0.001 C = 1
  rf = caretModelSpec('rf', tuneGrid=expand.grid(mtry=c(2,5,10,50,100,150,200,300,400,500))),
  lars2 = caretModelSpec('lars2', tuneGrid=expand.grid(step=c(5,10,20,50,100,250,275,300,350,400))),
  relaxo = caretModelSpec('relaxo', tuneGrid=expand.grid(lambda=seq(2,50, length.out=10), phi=seq(0.1,10, length.out=10))),
  ctree2 = caretModelSpec('ctree2', tuneGrid=expand.grid(maxdepth=seq(2,20, length.out=10), mincriterion=c(0.01,0.5, length.out=10)))
)

ensemble.trainControl = trainControl(
  method = 'LGOCV',
  number=5,
  p = 0.90,
  index=createFolds(DF.eff.split$targetsTrain, 5),
  summaryFunction = defaultSummary,
  savePredictions = 'final',
  classProbs = FALSE,
  verboseIter = TRUE,
  returnData = FALSE,
  allowParallel = DO_PARALLEL
)

ensemble.args <- list(
  metric = 'RMSE',
  x = DF.eff.split$inputsTrain,
  y = DF.eff.split$targetsTrain,
  trControl = ensemble.trainControl,
  tuneList = ensemble.tuneList,
  continue_on_fail = F
)

stack.trainControl = trainControl(
  method = 'LGOCV',
  number  = 5,
  savePredictions = 'final',
  classProbs = FALSE,
  summaryFunction = defaultSummary,
  returnData  = FALSE
)

stack.args <- list(
  method = 'glm',
  metric = 'Rsquared',
  trControl = stack.trainControl
)

model.effect.cont <- train.ensemble(
  DF.eff.split$inputsTrain, 
  DF.eff.split$targetsTrain, 
  DF.eff.split$inputsTest, 
  DF.eff.split$targetsTest,
  ensemble.args,
  stack.args
)