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
  DF.aff.split <- readRDS('data/datasets/Affinity.Continuous.Set3.rds')
}else{
  # Load the parsed dataset and prepare it
  DF.cont <- read.csv('data/parsed/GHSR.AFFINITY.CONTINUOUS.csv', stringsAsFactors = F, row.names=1)  
  # Apply descriptors and split train/test
  DF.aff.split <- prepare.dataset(
    DF.cont, 
    response.col='Standard.Value',
    desc.set=DESC_SET, 
    method='random', 
    remove.sparsity=F
  )
  if(SAVE_DSET){
    saveRDS(DF.aff.split, 'data/datasets/Affinity.Continuous.Set3.rds')  
  }
}

###########
## Model ##
###########
ensemble.tuneList = list(
  rf = caretModelSpec('rf', tuneGrid=expand.grid(mtry=c(2,4,6,8,10,12,20,30,50,100))),
  lars2 = caretModelSpec('lars2', tuneGrid=expand.grid(step=c(5,10,20,50,100,250,275,300,350,400))),
  glmnet = caretModelSpec('glmnet', tuneGrid=expand.grid(alpha=seq(0.01, 0.9,length.out=10), lambda=seq(0.01, 0.9, length.out=10))),
  relaxo = caretModelSpec('relaxo', tuneGrid=expand.grid(lambda=seq(2,50, length.out=10), phi=seq(0.001,0.9, length.out=10)))
)

ensemble.trainControl = trainControl(
  method = 'LGOCV',
  number=5,
  p = 0.90,
  index=createFolds(DF.aff.split$targetsTrain, 5),
  summaryFunction = defaultSummary,
  savePredictions = 'final',
  classProbs = FALSE,
  verboseIter = TRUE,
  returnData = FALSE,
  allowParallel = DO_PARALLEL
)

ensemble.args <- list(
  metric = 'RMSE',
  x = DF.aff.split$inputsTrain,
  y = DF.aff.split$targetsTrain,
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

stack.affinity.cont <- train.ensemble(
  DF.aff.split$inputsTrain, 
  DF.aff.split$targetsTrain, 
  DF.aff.split$inputsTest, 
  DF.aff.split$targetsTest,
  ensemble.args,
  stack.args
)