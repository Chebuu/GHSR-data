library(dplyr)

source('descriptors.R')
source('strategies.R')

DESC_SET = set3

source('models/categorical/ensemble.activity.R')
source('models/continuous/ensemble.affinity.R')
source('models/continuous/effect.R')

saveRDS(stack.activity.cat, 'Stack.Activity.Cat.Set3.rds')
saveRDS(stack.affinity.cont, 'Stck.Affinity.Cont.Set3.rds')
saveRDS(model.eff.cont, 'RF.Eff.Cont.Set3.rds')

tcmdb <- read.delim('data/datasets/tcmdb/compounds_final_4.txt', header=F, sep='\t', quote='')
tcmdb <- tcmdb[sample(1:nrow(tcmdb), size=50),]

tcmdb.subset <- tcmdb %>%
  select(
    V2, V3, V5, V6
  ) %>%
  filter(
    V2 != 'Not Available' & V3 != 'Not Available' & V5 != 'Not Available' & V6 != 'Not Available'
  )

X <- applyDescriptors(parse.smiles(as.vector(tcmdb.subset$V6)), descriptors=set3, remove.sparsity=FALSE)

result <- pred.models(
  X, 
  list(
    activity.cat = stack.activity.cat,
    affinity.cont = stack.affinity.cont,
    effect.cont = model.eff.cont
  )
)

results <- strategy.thresh(
  result
)

results <- strategy.sort(
  result, 
  c(
    'activity.cat',
    'effect.cont',
    'affinity.cont'
  )
)

