library(caret)
library(caretEnsemble)
library(RSNNS)
library(rcdk)
library(ModelMetrics)
library(tidyverse)
library(assertthat)
library(doParallel)

source('methods.R')
source('descriptors.R')

DO_PARALLEL = F

if(DO_PARALLEL){
  cl <- makecluster(4)
  registerDoParallel(cl)
}else if(!is.null(cl)){
  stopCluster(cl)
}

## ## ## 
# Data #
## ## ##

DF.effect <- read.csv('./data/parsed/GHSR.ACTIVITY.CATEGORICAL.csv', stringsAsFactors = F, row.names=1)

# Convert int to char
DF.effect <- DF.effect %>%
  filter(
    Activity %in% c(1,2)
  ) %>%
  mutate(
    Activity = case_when(
      Activity == 1 ~ 'Agonist',
      Activity == 2 ~ 'Antagonist'
    )  
  )

# Apply descriptors and split train/test
DF.effect.split <- prepare.dataset(
  DF.effect, 
  desc.set=set4, 
  method='random', 
  remove.sparsity=T
)

esbl.tuneList = list(
  # gpls = caretModelSpec('gpls', tuneLength=10), # Slow
  # blackboost = caretModelSpec('blackboost', tuneLength=2), # Slow
  # adaboost = caretModelSpec('adaboost', tuneLength=2), # Slow
  # AdaBoostM1 = caretModelSpec('AdaBoost.M1', tuneLength=2), # Slow
  # hda = caretModelSpec('hda', tuneLength=2), # Slow
  # nodeHarvest = caretModelSpec('nodeHarvest', tuneLength=2), # Slow
  # ordinalNet = caretModelSpec('ordinalNet', tuneLength=2), # Slow
  # svmLinearWeights2 = caretModelSpec('svmLinearWeights2', tuneLength=2), # NO CLASS PROBS
  # C50Cost = caretModelSpec('C5.0Cost', tuneLength=2), # NO CLASS PROBS
  # rpartCost = caretModelSpec('rpartCost', tuneLength=2), # Fast, NO CLASS PROBS
  pda = caretModelSpec('pda', tuneGrid=expand.grid(lambda=seq(0.0001,0.1,length.out=10))),
  svmLinearWeights = caretModelSpec('svmLinearWeights', tuneGrid = expand.grid(cost=seq(0.01,5,length.out=10), weight=c(0.1,5,length.out=10))),
  evtree = caretModelSpec('evtree', tuneLength=10),
  svmRadialSigma = caretModelSpec('svmRadialSigma', tuneLength=10),
  svmRadialCost = caretModelSpec('svmRadialCost', tuneLength=10),
  svmPoly = caretModelSpec('svmPoly', tuneLength=10),
  svmLinear2 = caretModelSpec('svmLinear2', tuneLength=10),
  C5.0Rules = caretModelSpec('C5.0Rules', tuneLength=10),
  C5.0Tree = caretModelSpec('C5.0Tree', tuneLength=10),
  JRip = caretModelSpec('JRip', tuneLength=10),
  C5.0 = caretModelSpec('C5.0', tuneLength=10),
  rf = caretModelSpec('rf', tuneLength=10)
)

esbl.trainControl = trainControl(
  method = 'LGOCV',
  number=10,
  p = 0.90,
  index=createFolds(DF.effect.split$targetsTrain, 10),
  summaryFunction = defaultSummary,
  savePredictions = 'final',
  classProbs = TRUE,
  verboseIter = TRUE,
  returnData = FALSE,
  allowParallel = DO_PARALLEL
)

ensemble.args <- list(
  metric = 'Kappa',
  x = DF.effect.split$inputsTrain,
  y = as.factor(DF.effect.split$targetsTrain),
  trControl = esbl.trainControl,
  tuneList = esbl.tuneList,
  continue_on_fail = F
)

stack.trainControl = trainControl(
  method = 'LGOCV',
  number  = 10,
  savePredictions = 'final',
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  returnData  = FALSE
)

stack.args <- list(
  method = 'rpart',
  metric = 'Kappa',
  trControl = stack.trainControl
)

predictions <- train.ensemble(
  DF.effect.split$inputsTrain, 
  DF.effect.split$targetsTrain, 
  DF.effect.split$inputsTest, 
  DF.effect.split$targetsTest,
  ensemble.args,
  stack.args
)

resamps <- resamples(
  esbl
)

preds <- predict(predictions, data.frame(DF.effect.split$inputsTest))
caret::confusionMatrix(as.factor(preds), as.factor(DF.effect.split$targetsTest))

as.factor(DF.effect.split$inputsTest)

preds

modelCor(resamps)

theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(resamps, layout = c(3, 1))