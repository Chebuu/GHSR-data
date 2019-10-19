TRY CUBIST !!!! https://arxiv.org/pdf/1803.06236.pdf


library(rcdk)
library(randomForest)
library(ModelMetrics)
library(tidyverse)
library(assertthat)

source('methods.R')
source('descriptors.R')

DF.cont <- read.csv('./data/parsed/GHSR.ACTIVITY.CONTINUOUS.csv')

##########################
#  Effect Random Forest  #
##########################
DF.eff <- DF.cont %>% 
  select(Molecule, Canonical.Smiles, Standard.Value)

eff.IAtom <- parse.smiles(as.vector(DF.eff$Canonical.Smiles))

par(mfrow=c(2,2))
y.eff <- DF.eff$Standard.Value
plot(y.eff, main='Raw EC50/IC50 values')
hist(y.eff, main='Raw EC50/IC50 values')

y.eff <- log(DF.eff$Standard.Value)
plot(y.eff, main='log(EC50/IC50) values')
hist(y.eff, main='log(EC50/IC50) values')

# Calculate molecular descriptors from smiles strings
X.eff <- applyDescriptors(
  smiles = eff.IAtom,
  descriptors = DESC_SET,
  remove.sparsity = T
)

assert_that(length(y.eff) == nrow(X.eff))
assert_that(!any(is.na(X.eff)))

train.idx <- sample(1:length(y.eff), floor(length(y.eff) * 0.9), replace = FALSE)
test.idx <- setdiff(1:length(y.eff), train.idx)
assert_that(!any(train.idx %in% test.idx))

X.eff.train <- X.eff[train.idx,]
y.eff.train <- y.eff[train.idx]

X.eff.test <- X.eff[test.idx,]
y.eff.test <- y.eff[test.idx]

model.eff.cont <- randomForest(data.frame(X.eff.train), y.eff.train, mtry=10, ntree=1000, nodesize = 10)
pred.eff <- predict(model.eff.cont, X.eff.test)

( model.eff.mae <- mae(y.eff.test, pred.eff) )
( model.eff.rmse <- rmse(y.eff.test, pred.eff) )
( model.eff.r2 <- r2(y.eff.test, pred.eff) )

par(mfrow=c(1,1))
plot(pred.eff, y.eff.test, xlim=c(-2, 10), ylim=c(-2, 10))
abline(0,1,col='blue')

# 
# #########################
# #  EC50  Random Forest  #
# #########################
# DF.ec50 <- DF.cont %>% 
#   filter(Standard.Type == 'EC50') %>%
#   select(Molecule, Canonical.Smiles, Standard.Value)
# 
# ec50.IAtom <- parse.smiles(as.vector(DF.ec50$Canonical.Smiles))
# 
# par(mfrow=c(2,2))
# y.ec50 <- DF.ec50$Standard.Value
# plot(y.ec50, main='Raw EC50 values')
# hist(y.ec50, main='Raw EC50 values')
# 
# y.ec50 <- log(DF.ec50$Standard.Value)
# plot(y.ec50, main='log(EC50) values')
# hist(y.ec50, main='log(EC50) values')
# 
# # Calculate molecular descriptors from smiles strings
# X.ec50 <- applyDescriptors(
#   smiles = ec50.IAtom, 
#   # descriptors = set3,
#   descriptors = set1,
#   remove.sparsity = F
# )
# 
# assert_that(length(y.ec50) == nrow(X.ec50))
# assert_that(!any(is.na(X.ec50)))
# 
# train.idx <- sample(1:length(y.ec50), floor(length(y.ec50) * 0.85), replace = FALSE)
# test.idx <- setdiff(1:length(y.ec50), train.idx)
# assert_that(!any(train.idx %in% test.idx))
# 
# X.ec50.train <- X.ec50[train.idx,]
# y.ec50.train <- y.ec50[train.idx]
# 
# X.ec50.test <- X.ec50[test.idx,]
# y.ec50.test <- y.ec50[test.idx]
# 
# model.ec50 <- randomForest(data.frame(X.ec50.train), y.ec50.train, mtry=10, ntree=15000, nodesize = 10)
# pred.ec50 <- predict(model.ec50, X.ec50.test)
# 
# ( model.ec50.mae <- mae(y.ec50.test, pred.ec50) )
# ( model.ec50.rmse <- rmse(y.ec50.test, pred.ec50) )
# ( model.ec50.r2 <- r2(y.ec50.test, pred.ec50) )
# 
# par(mfrow=c(1,1))
# plot(pred.ec50, y.ec50.test, xlim=c(-2, 10), ylim=c(-2, 10))
# abline(0,1,col='blue')
# 
# 
# #########################
# #  IC50  Random Forest  #
# #########################
# DF.ic50 <- DF.cont %>% 
#   filter(Standard.Type == 'IC50') %>%
#   select(Molecule, Canonical.Smiles, Standard.Value)
# 
# ic50.IAtom <- parse.smiles(as.vector(DF.ic50$Canonical.Smiles))
# 
# par(mfrow=c(2,2))
# y.ic50 <- DF.ic50$Standard.Value
# plot(y.ic50, main='Raw ic50 values')
# hist(y.ic50, main='Raw ic50 values')
# 
# y.ic50 <- log10(DF.ic50$Standard.Value)
# plot(y.ic50, main='log(ic50) values')
# hist(y.ic50, main='log(ic50) values')
# 
# # Calculate molecular descriptors from smiles strings
# X.ic50 <- applyDescriptors(
#   smiles = ic50.IAtom, 
#   descriptors = set3, 
#   remove.sparsity = T
# )
# 
# assert_that(length(y.ic50) == nrow(X.ic50))
# assert_that(!any(is.na(X.ic50)))
# 
# train.idx <- sample(1:length(y.ic50), floor(length(y.ic50) * 0.90), replace = FALSE)
# test.idx <- setdiff(1:length(y.ic50), train.idx)
# assert_that(!any(train.idx %in% test.idx))
# 
# X.ic50.train <- X.ic50[train.idx,]
# y.ic50.train <- y.ic50[train.idx]
# 
# X.ic50.test <- X.ic50[test.idx,]
# y.ic50.test <- y.ic50[test.idx]
# 
# model.ic50 <- randomForest(data.frame(X.ic50.train), y.ic50.train, mtry=20, ntree=30000, nodesize = 15)
# pred.ic50 <- predict(model.ic50, X.ic50.test)
# 
# ( model.ic50.mae <- mae(y.ic50.test, pred.ic50) )
# ( model.ic50.rmse <- rmse(y.ic50.test, pred.ic50) )
# ( model.ic50.r2 <- r2(y.ic50.test, pred.ic50) )
# 
# par(mfrow=c(1,1))
# plot(pred.ic50, y.ic50.test, xlim=c(0, 6), ylim=c(0, 6))
# abline(0,1,col='blue')
# 
