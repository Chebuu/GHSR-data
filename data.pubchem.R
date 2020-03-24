# TODO:
## What are the units of acvalue? uM? nM? Because some values are on the order of 10000 and some are on the order of 0.0001.
## What is the meaning of Kbapp? 
## ## - Is Kbapp proportional to Kb?
## ## ## - Reference the linked pulication for the aid reporting Kbapp
## Convert all Ka to Kd via Kd = 1/Ka
## Impute new feature 'direction' as agonist/antagonist/inactive 
## There are many duplicate compounds with non-unique acvalues
## ## - They need to be removed and the cause should be investigated
## ## ## - The call to unq.rows is probably not doing what you think
## Reformat columns to match those in parsed datasets
## Subset the dataset for modeling

source('util.R')

###
# Load raw datasets
###
DFPC.Q92847 <- read.csv('data/raw/pubchem/GHSR.ALL.Q92847.csv', stringsAsFactors = FALSE)
DFPC.2693 <- read.csv('data/raw/pubchem/GeneID_2693_bioactivity_gene.csv', stringsAsFactors = FALSE)[,colnames(DFPC.Q92847)]

###
# Extract unique samples from DFPC.2963 and bind the two datasets
###
# Convert <missing> to ''
DFPC.2693 <- DFPC.2693 %>%
  mutate_all(~blank_if(.,'NULL'))
DFPC.Q92847 <- DFPC.Q92847 %>%
  mutate_all(~blank_if_na(.))

# Find distinct rows from DFPC.2693
vals.compare <- c('activity', 'acname', 'acvalue')
unq.rows <- apply(DFPC.2693, 1, function(grow){
  all(
    apply(DFPC.Q92847, 1, function(prow){
      !all(prow == grow)
    })
  )
})
DFPC.2693 <- DFPC.2693 %>%
  filter(unq.rows)

###
# Merge datasets and reformat
###
DFPC <- rbind(DFPC.Q92847, DFPC.2693)

# # Unify kinetic terms (Kd = 1/Kb)
# DFPC %>%
#   mutate(
#     acvalue = case_when(
#       acname == 'Kd' ~ 1/acvalue,
#       TRUE ~ acvalue
#     ),
#     acname = case_when(
#       acname == 'Kd' ~ 'Kb',
#       TRUE ~ acname
#     )
#   )

###
# Imputation utilities
###
# Utilities
activityRange <- function(DF, target.activity, refs.acname=NULL){
  #' Calculate the range of activity values represented by a given activity type. For example, calculate the range of EC50 values represented by compounds classified as active.
  #' @param DF A data.frame of the pubchem bioactivity dataset
  #' @param target.activity A character representing the activities for which acvalue ranges should be calculated
  #' @param refs.acname The acname/s to range
  if(is.null(refs.acname)) refs.acname <- levels(as.factor(DF$acname))
  targets <- DFPC %>%
    filter(activity == target.activity & acname %in% refs.acname) %>%
    filter(!is.na(acvalue))
  sapply(refs.acname, function(acn){
    target <- filter(targets, acname == acn)
    return(
      c(
        min = min(target$acvalue),
        max = max(target$acvalue)
      )
    )
  })
}

imputeActivity <- function(DF, imp.lvl=c('Unspecified', 'Inconclusive'), act.rng=NULL){
  #' Impute missing activities using existing acvalues. By default, the function attempts to impute Unspecified and Inconclusive values using the reference ranges calculated by activityRange(DF, X) for each activity X in levels(as.factor(DF$activity)) .
  #' @section Note: Max values are inclusive, min values are exclusive.
  #' @param DF A data.frame of the pubchem bioactivity dataset
  #' @param imp.lvl The activities to be replaced with imputed values
  #' @param act.rng A list of data.frames used as references to impute missing values. For example act.rng = list(Active=data.frame(EC50=c(min=-Inf,max=10000)), Inactive=Active=data.frame(EC50=c(min=10001,max=Inf)))
  if(is.null(act.rng)){
    act.lvl <- levels(as.factor(DFPC$activity))
    act.rng <- setNames(lapply(act.lvl, function(act){
      activityRange(DFPC, act)
    }), act.lvl)
  }
  imp.idx <- DF$activity %in% imp.lvl
  ref.act <- names(act.rng)[!(names(act.rng) %in% imp.lvl)]
  acts.imp <- apply(DF[imp.idx,], 1, function(drow){
    act <- drow['activity']
    acn <- drow['acname']
    acv <- drow['acvalue']
    act.imp <- act
    if(!is.na(acn)){
      for(acrn in ref.act){
        acr <- act.rng[[acrn]]
        rmin <- tryCatch(acr['min', acn], error=function(x) NA)
        rmax <- tryCatch(acr['max', acn], error=function(x) NA)
        if(!any(is.na(c(acv, rmin, rmax)))){
          if(acv > rmin && acv <= rmax){
            act.imp <- acrn
          }  
        } 
      }
    }
    return(act.imp)
  })
  DF[imp.idx,'activity'] <- acts.imp
  return(DF)
}

###
# Impute missing activities by reference
###
# Convert <missing> to NA
DFPC <- DFPC %>%
  mutate_all(~na_if(., ''))

# # Create acvalue -> activity mappings
# act.levels <- unique(DFPC$activity)
# act.range <- setNames(lapply(act.levels, function(act){
#   activityRange(DFPC, act)
# }), act.levels)
# 
# act.range$Active['min',] <- -Inf
# act.range$Active['max', '%max'] <- Inf
# 
# act.range$Inactive['max',] <- Inf
# act.range$Inactive['min',] <- act.range$Active['max',]
# 
# ref.remove <- colnames(act.range$Active) %in% c('%max', 'Activity', 'E/Emax', 'Emax', 'Intrinsic activity')
# act.range$Active <- act.range$Active[,!ref.remove]
# act.range$Inactive <- act.range$Inactive[,!ref.remove]
# 
# # Impute activities
# DFPC <- imputeActivity(DFPC, act.rng=act.range)


###
# Impute activities by thresholding acvalues
###
# EC50_MAX_ACTIVE <- 10000
# IC50_MAX_ACTIVE <- 10000
# KB_MAX_ACTIVE <- 10
# KD_MAX_ACTIVE <- 10
# KI_MAX_ACTIVE <- 10
# 
# DFPC <- DFPC %>%
#   mutate(
#     activity = case_when(
#       acname == 'EC50' & acvalue > EC50_MAX_ACTIVE ~ 'Inactive',
#       acname == 'IC50' & acvalue > IC50_MAX_ACTIVE ~ 'Inactive',
#       acname == 'Kb' & acvalue > KB_MAX_ACTIVE ~ 'Inactive',
#       acname == 'Ki' & acvalue > KI_MAX_ACTIVE ~ 'Inactive',
#       TRUE ~ activity
#     )
#   )

# ###
# # Impute agonist/antagonist
# ###

  
# ###
# # Fetch SMILES
# ###
# DFPC <- DFPC %>%
#   left_join(
#     getSMILES.cid(DFPC$cid),
#     by=c('cid'='CID')
#   ) 

###
# Subset
###

DFPC <- DFPC %>%
  mutate_all(function(x) ifelse(x=='', NA_character_, x))

# INACTIVE subset
DFPC.inactive <- DFPC %>% 
  filter(activity=='Inactive') %>%
  transmute(Id=cid, Id.Type='cid') %>%
  unique(by='cid')
write.csv(DFPC.inactive, './data/parsed/pubchem/GHSR.INACTIVE.csv')

# TODO: ACTIVITY subset
DFPC.activity <- data.frame(
  Id=NA,
  Id.Type=NA, 
  Activity = NA,
  Agonist = NA,
  Partial = NA,
  Antagonist = NA,
  Invagonist = NA
)
write.csv(DFPC.activity, 'data/parsed/pubchem/GHSR.ACTIVITY.csv')

# EC50 subset
DFPC.ic50 <- DFPC %>% 
  filter(acname=='EC50') %>%
  transmute(Id=cid, Id.Type='cid', Standard.Value=acvalue)
write.csv(DFPC.ic50, 'data/parsed/pubchem/GHSR.EC50.csv')

# IC50 subset
DFPC.ic50 <- DFPC %>% 
  filter(acname=='IC50') %>%
  transmute(Id=cid, Id.Type='cid', Standard.Value=acvalue)
write.csv(DFPC.ic50, 'data/parsed/pubchem/GHSR.IC50.csv')

# KI subset
DFPC.ki <- DFPC %>%
  filter(acname=='Ki') %>%
  transmute(Id=cid, Id.Type='cid', Standard.Value=acvalue)
write.csv(DFPC.ki, 'data/parsed/pubchem/GHSR.KI.csv')

# KB subset
DFPC.kb <- DFPC %>%
  filter(acname=='Kb') %>%
  transmute(Id=cid, Id.Type='cid', Standard.Value=acvalue)
write.csv(DFPC.kb, 'data/parsed/pubchem/GHSR.KB.csv')

