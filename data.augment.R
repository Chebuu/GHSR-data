# # TODO::
# # # - All thresholds should be converted to ranges
# # # - All thresholds create duplicates because they loop over the entire reference set. 
# # # # # - Must extract unqiue compounds before appending.
# # # - Partials are measured as IC50/Ki not EC50/Kd
# 
# library(dplyr)
# 
# source('util.R')
# 
# # Load aggregate datasets
# GHSR.aggregate <- readGHSR.aggregate()
# 
# # Load parsed datasets
# GHSR.parsed <- readGHSR.parsed(bind=TRUE)
# 
# ############################
# # Infer Categorical Values #
# ############################
# THRESHOLD <- list(
#   Inactive=list(
#     EC50=list(
#       min=10000,
#       max=Inf
#     ),
#     IC50=list(
#       min=10000,
#       max=Inf
#     )
#   ),
#   Agonist=list(
#     EC50=list(
#       min=-Inf,
#       max=1000
#     ),
#     IC50=list(
#       min=Inf,
#       max=-Inf
#     )
#   ),
#   Partial=list(
#     EC50=list(
#       min=1000,
#       max=1000
#     ),
#     IC50=list(
#       min=Inf,
#       max=-Inf
#     )
#   ),
#   Invagonist=list(
#     EC50=list(
#       min=Inf,
#       max=-Inf
#     ),
#     IC50=list(
#       min=-Inf,
#       max=1
#     )
#   ),
#   Antagonist=list(
#     EC50=list(
#       min=Inf,
#       max=-Inf
#     ),
#     IC50=list(
#       min=1000,
#       max=10000
#     )
#   )
# )
# 1000 >= THRESHOLD$Inactive[['EC50']][['min']] && 1000 >= THRESHOLD$Inactive[['EC50']][['max']]
# threshold.which <- function(val, type, threshold=THRESHOLD){
#   sapply(threshold, function(act){
#     act.min <- act[[type]][['min']]
#     act.max <- act[[type]][['max']]
#     is.act <- val > act.min && val <= act.max
#     return(is.act)
#   })
# }
# threshold.which(1000, 'IC50')
# 
# inferCategories <- function(GHSR, cat.dsets, ref.dsets, threshold=THRESHOLD){
#   cats <- GHSR[cat.dsets]
#   refs <- GHSR[ref.dsets]
# }
# ## ##
# ## Agonist
# ## ##
# df.cat.agonist <- df.cat %>% filter(Agonist==1)
# 
# df.refs.act <- DFS$GHSR.ACTIVITY.CONTINUOUS.csv
# EC50.refs.idx <- cont.act.EC50$Molecule %in% df.cat.agonist$Molecule
# EC50.refs.mols <- cont.act.EC50$Molecule[EC50.refs.idx]
# EC50.refs.vals <- cont.act.EC50$Standard.Value[EC50.refs.idx]
# 
# df.refs.aff <- DFS$GHSR.AFFINITY.CONTINUOUS.csv
# Kb.refs.idx <- cont.aff.Kb$Molecules %in% df.cat.agonist$Molecule
# Kb.refs.mols <- cont.aff.Kb$Molecule[Kb.refs.idx]
# Kb.refs.vals <- cont.aff.Kb$Standard.Value[Kb.refs.idx]
# 
# # By max EC50
# EC50.max <- max(EC50.refs.vals)
# EC50.imp.idx <- EC50.refs.vals < EC50.max
# EC50.imp.mol <- EC50.refs.mols[EC50.imp.idx]
# EC50.imp.val <- rep(1, sum(EC50.imp.idx))
# df.imp.agonist <- df.cat.agonist %>%
#   add_row(
#     Molecule=EC50.imp.mol, 
#     Agonist=EC50.imp.val
#   )
# 
# # By threshold EC50
# ## !!! Creates duplicates
# ## !!! This should be a range, not a threshold
# AGO_EC50_RANGE <- c(-Inf, 9999)
# EC50.thr.imp.mol <- cont.act.EC50 %>%
#   filter(Standard.Value < THRESH_EC50) %>%
#   filter(!(Molecule %in% EC50.imp.mol)) %>%
#   select(Molecule)
# EC50.thr.imp.val <- rep(1, length(EC50.thr.imp.mol))
# df.imp.agonist <- df.imp.agonist %>%
#   add_row(
#     Molecule=EC50.thr.imp.mol, 
#     Agonist=EC50.thr.imp.val
#   )
# 
# # By max Kb
# Kb.max <- max(Kb.refs.vals) # Apparently none exist
# Kb.imp.idx <- Kb.refs.vals < Kb.max
# Kb.imp.mol <- Kb.refs.mols[Kb.imp.idx]
# Kb.imp.val <- rep(1, sum(Kb.imp.idx))
# df.imp.agonist <- df.imp.agonist %>%
#   add_row(
#     Molecule=Kb.imp.mol, 
#     Agonist=Kb.imp.val
#   )
# 
# # By range Kb
# ## !!! Creates duplicates
# ## !!! This should be a range, not a threshold
# THRESH_EC50 <- 10000
# EC50.thr.imp.mol <- cont.act.EC50 %>%
#   filter(Standard.Value < THRESH_EC50) %>%
#   filter(!(Molecule %in% EC50.imp.mol)) %>%
#   select(Molecule)
# EC50.thr.imp.val <- rep(1, length(EC50.thr.imp.mol))
# df.imp.agonist <- df.imp.agonist %>%
#   add_row(
#     Molecule=EC50.thr.imp.mol, 
#     Agonist=EC50.thr.imp.val
#   )
# 
# 
# ## ##
# ## ## Partial
# ## ##
# df.cat.partial <- df.cat %>% filter(Partial==1)
# 
# df.refs.act <- DFS$GHSR.ACTIVITY.CONTINUOUS.csv
# EC50.refs.idx <- cont.act.EC50$Molecule %in% df.cat.partial$Molecule
# EC50.refs.mols <- cont.act.EC50$Molecule[EC50.refs.idx]
# EC50.refs.vals <- cont.act.EC50$Standard.Value[EC50.refs.idx]
# # # !!! There are IC50 values for partial agonists
# # which(df.cat.partial$Molecule %in% df.refs.act$Molecule[df.refs.act$Standard.Type=='IC50'])
# 
# df.refs.aff <- DFS$GHSR.AFFINITY.CONTINUOUS.csv
# Kb.refs.idx <- cont.aff.Kb$Molecules %in% df.cat.partial$Molecule
# Kb.refs.mols <- cont.aff.Kb$Molecule[Kb.refs.idx]
# Kb.refs.vals <- cont.aff.Kb$Standard.Value[Kb.refs.idx]
# # # !!! There is one Ki value for a partial agonist
# # which(df.cat.partial$Molecule %in% df.refs.act$Molecule[df.refs.aff$Standard.Type=='Kb'])
# 
# # By max EC50
# EC50.max <- max(EC50.refs.vals)
# EC50.imp.idx <- EC50.refs.vals < EC50.max
# EC50.imp.mol <- EC50.refs.mols[EC50.imp.idx]
# EC50.imp.val <- rep(1, sum(EC50.imp.idx))
# df.imp.partial <- df.cat.partial %>%
#   add_row(
#     Molecule=EC50.imp.mol, 
#     Agonist=EC50.imp.val
#   )
# 
# # By max Kb
# THRESH_KB <- 10000
# Kb.max <- max(Kb.refs.vals)# Apparently none exist
# Kb.imp.idx <- Kb.refs.vals < Kb.max
# Kb.imp.mol <- Kb.refs.mols[Kb.imp.idx]
# Kb.imp.val <- rep(1, sum(Kb.imp.idx))
# df.imp.partial <- df.imp.partial %>%
#   add_row(
#     Molecule=Kb.imp.mol, 
#     Agonist=Kb.imp.val
#   )
# 
# 
# ## ##
# ## ## Antagonist
# ## ##
# df.cat.antagonist <- df.cat %>% filter(Antagonist==1)
# 
# df.refs.act <- DFS$GHSR.ACTIVITY.CONTINUOUS.csv
# IC50.refs.idx <- cont.act.IC50$Molecule %in% df.cat.antagonist$Molecule
# IC50.refs.mols <- cont.act.IC50$Molecule[IC50.refs.idx]
# IC50.refs.vals <- cont.act.IC50$Standard.Value[IC50.refs.idx]
# 
# df.refs.aff <- DFS$GHSR.AFFINITY.CONTINUOUS.csv
# Ki.refs.idx <- cont.aff.Ki$Molecules %in% df.cat.antagonist$Molecule
# Ki.refs.mols <- cont.aff.Ki$Molecule[Ki.refs.idx]
# Ki.refs.vals <- cont.aff.Ki$Standard.Value[Ki.refs.idx]
# 
# # By max IC50
# IC50.max <- max(IC50.refs.vals)
# IC50.imp.idx <- IC50.refs.vals < IC50.max
# IC50.imp.mol <- IC50.refs.mols[IC50.imp.idx]
# IC50.imp.val <- rep(1, sum(IC50.imp.idx))
# df.imp.antagonist <- df.cat.antagonist %>%
#   add_row(
#     Molecule=IC50.imp.mol, 
#     Agonist=IC50.imp.val
#   )
# 
# # By threshold IC50
# ## !!! Creates duplicates
# ## !!! This should be a range, not a threshold
# THRESH_IC50 <- 10000
# IC50.thr.imp.mol <- cont.act.IC50 %>%
#   filter(Standard.Value < THRESH_IC50) %>%
#   filter(!(Molecule %in% IC50.imp.mol)) %>%
#   select(Molecule)
# IC50.thr.imp.val <- rep(1, length(IC50.thr.imp.mol))
# df.imp.antagonist <- df.imp.antagonist %>%
#   add_row(
#     Molecule=IC50.thr.imp.mol, 
#     Agonist=IC50.thr.imp.val
#   )
# 
# # By max Ki
# Ki.max <- max(Ki.refs.vals)# Apparently none exist
# Ki.imp.idx <- Ki.refs.vals < Kb.max
# Ki.imp.mol <- Ki.refs.mols[Ki.imp.idx]
# Ki.imp.val <- rep(1, sum(Ki.imp.idx))
# df.imp.antagonist <- df.imp.antagonist %>%
#   add_row(
#     Molecule=Ki.imp.mol, 
#     Agonist=Ki.imp.val
#   )
# 
# # By threshold Ki
# ## !!! Creates duplicates
# ## !!! This should be a range, not a threshold
# THRESH_KI <- 10000
# Ki.thr.imp.mol <- cont.act.Ki %>%
#   filter(Standard.Value < THRESH_KI) %>%
#   filter(!(Molecule %in% Ki.imp.mol)) %>%
#   select(Molecule)
# Ki.thr.imp.val <- rep(1, length(Ki.thr.imp.mol))
# df.imp.agonist <- df.imp.antagonist %>%
#   add_row(
#     Molecule=Ki.thr.imp.mol, 
#     Agonist=Ki.thr.imp.val
#   )
# 
# 
# ## ##
# ## ## Invagonist
# ## ##
# df.cat.invagonist <- df.cat %>% filter(Invagonist==1)
# 
# df.refs.act <- DFS$GHSR.ACTIVITY.CONTINUOUS.csv
# IC50.refs.idx <- cont.act.IC50$Molecule %in% df.cat.invagonist$Molecule
# IC50.refs.mols <- cont.act.IC50$Molecule[IC50.refs.idx]
# IC50.refs.vals <- cont.act.IC50$Standard.Value[IC50.refs.idx]
# 
# df.refs.aff <- DFS$GHSR.AFFINITY.CONTINUOUS.csv
# Ki.refs.idx <- cont.aff.Ki$Molecules %in% df.cat.invagonist$Molecule
# Ki.refs.mols <- cont.aff.Ki$Molecule[Ki.refs.idx]
# Ki.refs.vals <- cont.aff.Ki$Standard.Value[Ki.refs.idx]
# 
# # By max IC50
# IC50.max <- max(IC50.refs.vals)
# IC50.imp.idx <- IC50.refs.vals < IC50.max
# IC50.imp.mol <- IC50.refs.mols[IC50.imp.idx]
# IC50.imp.val <- rep(1, sum(IC50.imp.idx))
# df.imp.invagonist <- df.cat.invagonist %>%
#   add_row(
#     Molecule=IC50.imp.mol, 
#     Agonist=IC50.imp.val
#   )
# 
# # By threshold IC50
# ## !!! Creates duplicates
# ## !!! This should be a range, not a threshold
# THRESH_IC50 <- 10000
# IC50.thr.imp.mol <- cont.act.IC50 %>%
#   filter(Standard.Value < THRESH_IC50) %>%
#   filter(!(Molecule %in% IC50.imp.mol)) %>%
#   select(Molecule)
# IC50.thr.imp.val <- rep(1, length(IC50.thr.imp.mol))
# df.imp.invagonist <- df.imp.invagonist %>%
#   add_row(
#     Molecule=IC50.thr.imp.mol, 
#     Agonist=IC50.thr.imp.val
#   )
# 
# # By max Ki
# Ki.max <- max(Ki.refs.vals)# Apparently none exist
# Ki.imp.idx <- Ki.refs.vals < Kb.max
# Ki.imp.mol <- Ki.refs.mols[Ki.imp.idx]
# Ki.imp.val <- rep(1, sum(Ki.imp.idx))
# df.imp.invagonist <- df.imp.invagonist %>%
#   add_row(
#     Molecule=Ki.imp.mol, 
#     Agonist=Ki.imp.val
#   )
# 
# # By threshold Ki
# ## !!! Creates duplicates
# ## !!! This should be a range, not a threshold
# THRESH_KI <- 10000
# Ki.thr.imp.mol <- cont.act.Ki %>%
#   filter(Standard.Value < THRESH_KI) %>%
#   filter(!(Molecule %in% Ki.imp.mol)) %>%
#   select(Molecule)
# Ki.thr.imp.val <- rep(1, length(Ki.thr.imp.mol))
# df.imp.agonist <- df.imp.invagonist %>%
#   add_row(
#     Molecule=Ki.thr.imp.mol, 
#     Agonist=Ki.thr.imp.val
#   )