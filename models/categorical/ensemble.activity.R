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
  DF.binary.split <- readRDS('data/datasets/Activity.Binary.Set3.rds')
}else{
  # Load the parsed dataset and prepare it
  DF.binary <- read.csv('data/parsed/GHSR.ACTIVITY.CATEGORICAL.BINARY.csv', stringsAsFactors = F, row.names=1)  
  # Convert numeric to character
  DF.binary <- DF.binary %>%
    mutate(
      Activity = case_when(
        Activity == 0 ~ 'inactive',
        Activity == 1 ~ 'active'
      )
    )
  # Apply descriptors and split train/test
  DF.binary.split <- prepare.dataset(
    DF.binary, 
    desc.set=DESC_SET, 
    method='random', 
    remove.sparsity=F
  )
  if(SAVE_DSET){
    saveRDS(DF.binary.split, 'data/datasets/Activity.Binary.Set3.rds')  
  }
}

###########
## Model ##
###########
esbl.tuneList = list(
  pda = caretModelSpec('pda', tuneGrid=expand.grid(lambda=seq(0.0001,0.1,length.out=10))),
  evtree = caretModelSpec('evtree', tuneLength=10),
  C5.0Rules = caretModelSpec('C5.0Rules', tuneLength=10),
  JRip = caretModelSpec('JRip', tuneLength=10)
)

esbl.trainControl = trainControl(
  method = 'LGOCV',
  number=5,
  p = 0.90,
  index=createFolds(DF.binary.split$targetsTrain, 5),
  summaryFunction = twoClassSummary,
  savePredictions = 'final',
  classProbs = TRUE,
  verboseIter = TRUE,
  returnData = FALSE,
  allowParallel = DO_PARALLEL
)

ensemble.args <- list(
  metric = 'ROC',
  x = DF.binary.split$inputsTrain,
  y = as.factor(DF.binary.split$targetsTrain),
  trControl = esbl.trainControl,
  tuneList = esbl.tuneList,
  continue_on_fail = F
)

stack.trainControl = trainControl(
  method = 'LGOCV',
  number  = 5,
  savePredictions = 'final',
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  returnData  = FALSE
)

stack.args <- list(
  method = 'gbm',
  metric = 'ROC',
  trControl = stack.trainControl
)

stack.activity.cat <- train.ensemble(
  DF.binary.split$inputsTrain, 
  DF.binary.split$targetsTrain, 
  DF.binary.split$inputsTest, 
  DF.binary.split$targetsTest,
  ensemble.args,
  stack.args
)
