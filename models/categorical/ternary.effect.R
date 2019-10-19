
library(caret)
library(kernalb)
library(rcdk)
library(randomForest)
library(ModelMetrics)
library(tidyverse)
library(assertthat)

source('methods.R')
source('descriptors.R')

registerDoParallel(4)


## ## ## 
# Data #
## ## ##

DF.ternary <- read.csv('data/parsed/GHSR.ACTIVITY.CATEGORICAL.TERNARY.csv', stringsAsFactors = F, row.names=1)
DF.ternary <- DF.ternary[sample(1:nrow(DF.ternary), nrow(DF.ternary)), sample(1:ncol(DF.ternary), ncol(DF.ternary))]
DF.ternary <- DF.ternary %>%
  mutate(
    Activity = case_when(
      Activity == 0 ~ 'NotActive',
      Activity == 1 ~ 'Agonist',
      Activity == 2 ~ 'Antagonist'
    )  
  )

ternary.IAtom <- parse.smiles(as.vector(DF.ternary$Canonical.Smiles))
X <- applyDescriptors(ternary.IAtom, set3)
y <- DF.ternary$Activity

assert_that(length(y) == nrow(X))
assert_that(!any(is.na(X)))

train.idx <- sample(1:length(y), floor(length(y) * 0.85), replace = FALSE)
test.idx <- setdiff(1:length(y), train.idx)
assert_that(!any(train.idx %in% test.idx))

X.train <- X[train.idx,]
y.train <- y[train.idx]

X.test <- X[test.idx,]
y.test <- y[test.idx]

model <- lssvm(x=as.matrix(X.train), y=as.factor(y.train), scaled = F, kernel = "rbfdot", kpar = "automatic",
               type = 'classification', tau = 0.001, reduced = TRUE, tol = 0.0001,
               rank = floor(dim(X)[1]/8), delta = 40, cross = 3, fit = F,
               subset, na.action = na.omit)
pred. <- predict(model, as.matrix(X.test))
caret::confusionMatrix(pred., as.factor(y.test))

lssvm(kernelMatrix(rbfdot(0.03), X.train))

( model.mae <- mae(y..test, pred.) )
( model.rmse <- rmse(y.test, pred.) )
( model.r2 <- r2(y..test, pred.) )

par(mfrow=c(1,1))
plot(pred., y..test, xlim=c(-2, 10), ylim=c(-2, 10))
abline(0,1,col='blue')

#########################
#  IC50  Random Forest  #
#########################
DF.ic50 <- DF.cont %>% 
  filter(Standard.Type == 'IC50') %>%
  select(Molecule, Canonical.Smiles, Standard.Value)

ic50.IAtom <- parse.smiles(as.vector(DF.ic50$Canonical.Smiles))

par(mfrow=c(2,2))
y.ic50 <- DF.ic50$Standard.Value
plot(y.ic50, main='Raw ic50 values')
hist(y.ic50, main='Raw ic50 values')

y.ic50 <- log10(DF.ic50$Standard.Value)
plot(y.ic50, main='log(ic50) values')
hist(y.ic50, main='log(ic50) values')

# Calculate molecular descriptors from smiles strings
X.ic50 <- applyDescriptors(
  smiles = ic50.IAtom, 
  descriptors = set3, 
  remove.sparsity = T
)

assert_that(length(y.ic50) == nrow(X.ic50))
assert_that(!any(is.na(X.ic50)))

train.idx <- sample(1:length(y.ic50), floor(length(y.ic50) * 0.90), replace = FALSE)
test.idx <- setdiff(1:length(y.ic50), train.idx)
assert_that(!any(train.idx %in% test.idx))

X.ic50.train <- X.ic50[train.idx,]
y.ic50.train <- y.ic50[train.idx]

X.ic50.test <- X.ic50[test.idx,]
y.ic50.test <- y.ic50[test.idx]

model.ic50 <- randomForest(data.frame(X.ic50.train), y.ic50.train, mtry=20, ntree=30000, nodesize = 15)
pred.ic50 <- predict(model.ic50, X.ic50.test)

( model.ic50.mae <- mae(y.ic50.test, pred.ic50) )
( model.ic50.rmse <- rmse(y.ic50.test, pred.ic50) )
( model.ic50.r2 <- r2(y.ic50.test, pred.ic50) )

par(mfrow=c(1,1))
plot(pred.ic50, y.ic50.test, xlim=c(0, 6), ylim=c(0, 6))
abline(0,1,col='blue')







