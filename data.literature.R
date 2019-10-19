# TODO:
## Pubchem query throws error
## Pbchem query returns NA values


# Parse the dataset collected from the literature

library(tidyverse)
library(rlang)
require(RCurl)
require(XML)
library(assertthat)

source('util.R')

# # Get all parsed datasets
# GHSR_PARSED_DIR <- './data/parsed'
# GHSR_PARSED_FILES <- dir(GHSR_PARSED_DIR)
# GHSR_PARSED_PATHS <- sapply(GHSR_PARSED_FILES, function(file) normalizePath(paste0(GHSR_PARSED_DIR,'/',file)))
# GHSR.all <- sapply(GHSR_PARSED_PATHS, function(file) read.csv(file, row.names=1))

# Get the dataset curated from the literature
GHSR_LITERATURE_DIR <- './data/literature'
GHSR_LITERATURE_PATH <- normalizePath(paste0(GHSR_LITERATURE_DIR, '/', 'literature.csv'))
GHSR.lit <- read.csv(GHSR_LITERATURE_PATH, stringsAsFactors = FALSE)

# Trim whitespace on all values
GHSR.lit <- GHSR.lit %>%
  mutate_all(trimws)

# Unify missing values as NA
GHSR.lit <- GHSR.lit %>%
  mutate_all(~na_if(., '')) %>%
  mutate_all(~na_if(., FALSE))

# Move activities from numeric columns to Activity column
inactive.m <- lapply(1:nrow(GHSR.lit), function(n){
  lrow <- GHSR.lit[n,]
  which(tolower(lrow) %in% c('inactive', 'not active'))
})
inactive.n <- unlist(lapply(inactive.m, function(m){
  length(m) > 0
}))
GHSR.lit[inactive.n, 'Activity'] <- 'Not Active'

# #############
# # One-liner #
# #############
# # This single block will replace all blocks below (except the final SID block)
# # But the verbose blocks below are more descriptive, so I kept them
# idcols <- c('InChIKey', 'CHEMBLID', 'CID')
# fetchby <- setNames(c(getSID.inchi, getSID.chembl, getSID.cid))
# GHSR.lit$SID <- fetchSIDsByColumn(GHSR.lit, idcols, fetchby)

###############
# Add Id col # 
###############
fill.id <- rep(NA, nrow(GHSR.lit))
GHSR.lit$Id <- fill.id
GHSR.lit$Id.Type <- fill.id

# CID
fill.idx <- !is.na(GHSR.lit$CID)
GHSR.lit$Id[fill.idx] <- GHSR.lit$CID[fill.idx]
GHSR.lit$Id.Type[fill.idx] <- rep('cid', sum(fill.idx))

# InChIKey
fill.idx <- !is.na(GHSR.lit$InChlKey) & is.na(GHSR.lit$Id)
GHSR.lit$Id[fill.idx] <- GHSR.lit$InChlKey[fill.idx]
GHSR.lit$Id.Type[fill.idx] <- rep('inchikey', sum(fill.idx))

# CHEMBLID
fill.idx <- !is.na(GHSR.lit$CHEMBLID) & is.na(GHSR.lit$Id)
GHSR.lit$Id[fill.idx] <- GHSR.lit$CHEMBLID[fill.idx]
GHSR.lit$Id.Type[fill.idx] <- rep('chemblid', sum(fill.idx))

# SID
fill.idx <- !is.na(GHSR.lit$SID) & is.na(GHSR.lit$Id)
GHSR.lit$Id[fill.idx] <- GHSR.lit$SID[fill.idx]
GHSR.lit$Id.Type[fill.idx] <- rep('sid', sum(fill.idx))

###################
# Inactive Subset #
###################
# TODO::
## Some compounds are listed as inactive in IC50 assay and active in EC50 assay, need to find those that are inactive in both
## Go through the literature and find an ID for these inactive compounds, it is NA right now
GHSR.lit.inactive <- GHSR.lit %>%
  filter(Activity == 'Not Active') %>%
  transmute(
    Id = case_when(
      !is.na(CID) ~ CID,
      !is.na(SID) ~ SID,
      !is.na(InChlKey) ~ InChlKey
    ),
    Id.Type = case_when(
      !is.na(SID) ~ 'cid',
      !is.na(CID) ~ 'sid',
      !is.na(CHEMBLID) ~ 'chemblid',
      TRUE ~ 'inchikey'
    )
  )
write.csv(GHSR.lit.inactive, './data/parsed/literature/GHSR.INACTIVE.csv')
  
###################
# Activity Subset #
###################
GHSR.lit.activity <- GHSR.lit %>%
  mutate(
    Agonist = case_when(
      Activity %in% c('fag','xag') ~ 1,
      TRUE ~ 0
    ),
    Partial = case_when(
      Activity == 'pag'            ~ 1,
      TRUE ~ 0
    ),
    Antagonist = case_when(
      Activity %in% c('xan','can') ~ 1,
      TRUE ~ 0
    ),
    Invagonist = case_when(
      Activity == 'iag'            ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    Activity = case_when(
      Activity %in% c('fag','xag','pag') ~ 1,
      Activity %in% c('xan','can','iag') ~ 0,
      TRUE ~ as.double(NA)
    )
  ) %>%
  select( Id, Id.Type, Activity, Agonist, Partial, Antagonist, Invagonist)

write.csv(GHSR.lit.activity, 'data/parsed/literature/GHSR.ACTIVITY.csv')


######################
# Continuous Subsets #
######################

toNanoMolar <- function(x){
  valid.units <- 'pM|nM|uM|mM|M'
  val <- gsub(valid.units, '', x)
  unit <- gsub(val, '', x)
  is.valid.unit <- grepl(valid.units, unit)
  if(!is.valid.unit) unit <- 'nM'
  multiplier <- c(
    pM = 1e-3,
    nM = 1e0,
    uM = 1e3,
    mM = 1e6,
    M = 1e9
  )[[unit]]
  num <- suppressWarnings(as.numeric(val))
  num * multiplier
}

# EC50 subset
GHSR.lit.ec50 <- GHSR.lit %>%
  filter(!is.na(EC50)) %>%
  mutate(Standard.Value = unlist(lapply(EC50, toNanoMolar))) %>%
  select(Id, Id.Type, Standard.Value)
write.csv(GHSR.lit.ec50, 'data/parsed/literature/GHSR.EC50.csv')

# IC50 subset
GHSR.lit.ic50 <- GHSR.lit %>%
  filter(!is.na(IC50)) %>%
  mutate(Standard.Value = unlist(lapply(IC50, toNanoMolar))) %>%
  select(Id, Id.Type, Standard.Value)
write.csv(GHSR.lit.ic50, 'data/parsed/literature/GHSR.IC50.csv')

# Ki
GHSR.lit.ki <- GHSR.lit %>%
  filter(!is.na(Ki)) %>%
  mutate(Standard.Value = unlist(lapply(Ki, toNanoMolar))) %>%
  select(Id, Id.Type, Standard.Value)
GHSR.lit.ki$Standard.Value
write.csv(GHSR.lit.ki, 'data/parsed/literature/GHSR.KI.csv')

# Kb 
GHSR.lit.kb <- GHSR.lit %>%
  filter(!is.na(Kb)) %>%
  mutate(Standard.Value = unlist(lapply(Kb, toNanoMolar))) %>%
  select(Id, Id.Type, Standard.Value)
write.csv(GHSR.lit.kb, 'data/parsed/literature/GHSR.KB.csv')


