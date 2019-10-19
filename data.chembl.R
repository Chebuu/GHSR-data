##############
##  TODO::  ##
##############
## CATEGORICAL DATASETS ##
# Activity dataset:
## - A Comment of 'Not active' does not mean the compound does not bind. It just means that the compound does not produce a response. In other words, these assays are not radioligand displacements. They measure the magnitude of a metabolic effect.
## ## - This begs the question: Does 0% activity mean the compound is a neutral antagonist, or is it inactive?
## ## ## - I think all assays measure % activity compared to a full agonist. So I'm guessing 0% activity means the tested compound brings the activity down to the constitutive activity of the receptor, and not actually 0% absolute receptor activity. The problem is that both a neutral antagonist and a non-binding compound would bring the activity to constitutive levels (an effective ligand is absent), and this would be measured as 0% activation compared to a full agonist. In other words, I can't distinguish antagonist from inactive compounds with a measured activity of 0%.
## ## ## - Unfortunately, there are no samlples with 0% activity that have a non-empty Comment field to reference. In other words, no 0% activity samples are labeled as antagonist or Not active; they are simply not labeled.
## ## ## - I have decided to drop these samples. If I want to add them, I will have to dig through the literature and figure out the meaning of 0% activation.
# EC50 and IC50 datasets:
## Does Not Determined mean their test was inconclusive or that the compound is ineffective at the receptor?
## ## - May be able to determine this by finding matching CHEMBLIDs in the activity dataset
# Kb dataset:
## As of now, all rows are filtered by assay and asssign a direction of 'Antagonist'
## ## - This may be problematic, because the dataset does not explicitly state that the tested compunds are known antagonists
## ## - However two points support the assignment of antagonist to all values in the Kb dataset
## ## ## - 1. A left_join on Molecules between Kb dataset and IC50 dataset shows all matched molecules have an IC50 value, whereas all molecules in the EC50 dataset have an NA value for the same left_join operation.
## ## ## - 2. My understanding is that the Kb traditionally refers to the dissociation constant of antagonist and receptor.
# Ki dataset:


## CONTINUOUS DATASET ##
# - Include temperatures and calculate deltaG as dependent variables (Ki,Kd) for samples that have kinetic constants 
# - It may be possible to add a couple continuous samples to the dataset by extracting from the categorical datasets, because many samples in these datasets have continuous variables, but I havent checked if these are just duplicates.


library(tidyverse)
library(assertthat)

#####################
##  Load Datasets  ##
#####################
GHSR_DIR <- './data/raw/chembl'
GHSR_FILES <- dir(GHSR_DIR)
cols.keep <- c('Molecule', 'Standard.Type', 'Standard.Value', 'Standard.Units', 'Comment')

# Load datasets
GHSR.list <- setNames(lapply(GHSR_FILES, function(ghsrf){
  fpath <- normalizePath(paste0(GHSR_DIR, '/', ghsrf))
  df <- read.csv(fpath, stringsAsFactors = FALSE)
  df %>% 
    select(cols.keep) %>%
    mutate(Id.Type='chemblid')
}), GHSR_FILES)


########################
##  Inactive Dataset  ##
########################
GHSR.list.inactive <- lapply(GHSR.list, function(dset){
  dset %>%
    filter(Comment == 'Not Active') %>%
    transmute(
      Id=Molecule,
      Id.Type='chemblid'
    )
})
GHSR.inactive <- do.call(rbind, GHSR.list.inactive) %>% unique(by='Id')
write.csv(GHSR.inactive, './data/parsed/chembl/GHSR.INACTIVE.csv')

############################
##  Categorical Datasets  ##
############################
###
# Activity dataset
###
GHSR.activity <- GHSR.list[[1]]

# Convert decimal to percent
GHSR.activity <- GHSR.activity %>%
  mutate(
    Standard.Value = as.numeric(Standard.Value),
    # The assay reporting decimal activations states that measurements are relative to compound 'T1' in their paper
    Standard.Value = case_when(
      Standard.Units == '' & Comment == '' ~ Standard.Value * 100,
      TRUE ~ Standard.Value
    ), 
    # I'm clipping values over 100% because no other assay in the dataset reports activations > 100%
    Standard.Value = ifelse(Standard.Value > 100, 100, Standard.Value)
  )

# # Drop 0% activations, because I can't determine antagonist vs inactive
# GHSR.activity <- GHSR.activity %>%
#   filter(!(Comment == '' & Standard.Value == 0))

# # Threshold % activity to infer direction
# GHSR.activity <- GHSR.activity %>%
#   mutate(
#     Comment = case_when(
#       (Comment == '' | Comment == 'Active') & (Standard.Value >= 50) ~ 'Agonist',
#       (Comment == '' | Comment == 'Active') & (Standard.Value > 0 & Standard.Value < 50) ~ 'Partial agonist',
#       (Comment == '' | Comment == 'Active') & (Standard.Value < 0) ~ 'Inverse agonist',
#       TRUE ~ Comment
#     )
#   )

# # Crudely estimte % activity from direction
# GHSR.activity <- GHSR.activity %>%
#   mutate(
#     Standard.Value = case_when(
#       is.na(Standard.Value) & (Comment == 'Agonist') ~ 75,
#       is.na(Standard.Value) & (Comment == 'Partial agonist') ~ 25,
#       is.na(Standard.Value) & (Comment == 'Inverse agonist') ~ -10,
#       TRUE ~ Standard.Value
#     )
#   )

###
# EC50 dataset
###
GHSR.ec50 <- GHSR.list[[2]]

GHSR.ec50 <- GHSR.ec50 %>% 
  mutate(
    Comment = 'Agonist'
  )
  # %>% mutate(
  #   Comment = case_when(
  #     Standard.Value >= 10000 ~ 'Not active'
  #     Standard.Value > 1000 & Standard.Value < 10000 ~ 'Partial agonist'
  #     TRUE ~ as.character(Comment)
  #   )
  # )
###
# IC50 dataset
###
GHSR.ic50 <- GHSR.list[[3]]

GHSR.ic50 <- GHSR.ic50 %>%
  mutate(
    Comment = 'Antagonist'
  )
# %>% mutate(
#   Comment = case_when(
#     Standard.Value >= 10000 ~ 'Not active'
#     TRUE ~ as.character(Comment)
#   )
# )

###
# Inhibition dataset
###
GHSR.inhibition <- GHSR.list[[4]]

GHSR.inhibition <- GHSR.inhibition %>%
  filter(Standard.Type=='Inhibition') %>%
  mutate(
    Comment = case_when(
      !suppressWarnings(is.na(as.numeric(Comment))) &
        (Comment == '') & (Standard.Units == '%') & (Standard.Value > 0)  ~ 'Antagonist',
      !suppressWarnings(is.na(as.numeric(Comment))) &
        (Comment == '') & (Standard.Units == '%') & (Standard.Value < 0)  ~ 'Inverse agonist',
      !suppressWarnings(is.na(as.numeric(Comment))) &
        (Comment == '') & (Standard.Units == 'nM') & (Standard.Value > 0) ~ 'Antagonist', # With this line, even weak affinities are considered to still be antagonists, because there are no inverse agonists with molar units in this dataset. 
      Comment == 'Active' ~ 'Antagonist',
      TRUE ~ as.character(Comment)
    )
  )

###
# Kb dataset
###
GHSR.kb <- GHSR.list[[5]]

# GHSR.kb <- GHSR.kb
# 
# left_join(GHSR.kb, GHSR.list[[1]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# left_join(GHSR.kb, GHSR.activity, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# 
# left_join(GHSR.kb, GHSR.list[[2]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# left_join(GHSR.kb, GHSR.ec50, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# 
# left_join(GHSR.kb, GHSR.list[[3]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# left_join(GHSR.kb, GHSR.ic50, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# 
# left_join(GHSR.kb, GHSR.list[[4]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y','Standard.Type')]
# left_join(GHSR.kb, GHSR.inhibition, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]

###
# Ki dataset
###
GHSR.ki <- GHSR.list[[6]]

# GHSR.ki <- GHSR.ki
# 
# left_join(GHSR.ki, GHSR.list[[1]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# left_join(GHSR.ki, GHSR.activity, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# 
# left_join(GHSR.ki, GHSR.list[[2]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# left_join(GHSR.ki, GHSR.ec50, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# 
# left_join(GHSR.ki, GHSR.list[[3]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# left_join(GHSR.ki, GHSR.ic50, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# 
# left_join(GHSR.ki, GHSR.list[[4]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# left_join(GHSR.ki, GHSR.inhibition, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# 
# left_join(GHSR.ki, GHSR.list[[5]], by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]
# left_join(GHSR.ki, GHSR.kb, by='Molecule')[,c('Molecule','Comment.x','Standard.Value.x','Comment.y','Standard.Value.y')]


###
# Clean and write
###
CAT.list <- list(
  GHSR.activity,
  GHSR.ec50,
  GHSR.ic50,
  GHSR.inhibition,
  GHSR.kb,
  GHSR.ki
)

# # Fill NA by reference
# FNABR <- do.call(
#   rbind,
#   lapply(GHSR.list, function(cdf){
#     do.call(
#       rbind,
#       lapply(GHSR.list, function(gdf){
#         left_join(cdf, gdf, by='Molecule')[
#           c('Molecule','Comment.x','Comment.y')
#         ] %>%
#           mutate_all(function(x) ifelse(x=='',NA_character_, x)) %>%
#           
#           filter(
#             (is.na(Comment.x) | Comment.x == 'Not Determined') &
#               !(is.na(Comment.y) | Comment.y == 'Not Determined')
#           ) %>%
#           mutate(Comment=Comment.y) %>%
#           
#           # # This returns 5 more data points, use either block...
#           # filter(!(is.na(Comment.x) & is.na(Comment.y))) %>%
#           # filter(!(is.na(Comment.x) & Comment.y == 'Not Determined')) %>%
#           # filter(!(is.na(Comment.y) & Comment.x == 'Not Determined')) %>%
#           # mutate(
#           #   Comment=case_when(
#           #     is.na(Comment.x) | Comment.x == 'Not Determined' &
#           #       !is.na(Comment.y) & (Comment.y != 'Not Determined') ~ Comment.y,
#           #     TRUE ~ Comment.x # Default to the value in CAT.list
#           #   )
#           # ) %>%
#           
#           select(Molecule, Comment) %>%
#           filter( !is.na(Comment) & Comment != 'Not Determined') %>%
#           filter(suppressWarnings(is.na(as.numeric(Comment)))) # Remove comments encoded as numbers.
#           
#       })
#     )
#   })
# )

# (yyy <- FNABR %>% # Print statements
#   group_by(Molecule, Comment) %>% 
#   summarise(count = length(Comment)))
# (zzz <- yyy %>%
#   group_by(Molecule) %>%
#   summarize(count=length(Molecule)))
# (mmm <- zzz %>%
#   filter(count > 1))
# (FNABR %>%
#   filter(Molecule %in% mmm$Molecule) %>%
#   arrange(Molecule))

# FNABR <- FNABR %>%
#   select(Molecule, Comment) %>%
#   unique(by='Molecule')
# 
# lapply(CAT.list, function(cdf){
#   # TODO::
#   ## - Replace Comment in cdf with corresponding comment in FNABR
#   cdf <- mutate_all(cdf, function(x) ifelse(x=='',NA_character_, x))
#   return(cdf)
#   
# })

# Aggregate all categorical datasets
GHSR.activity <- do.call(rbind, CAT.list) %>%
  mutate_all(function(x) ifelse(x=='', NA_character_, x)) %>%
  unique(by='Molecule')

# Binarize activity, one-hot encode direction, change ID column name
GHSR.activity <- GHSR.activity %>%
  transmute(
    Id = Molecule,
    Id.Type = Id.Type,
    Activity = ifelse(Comment %in% c('Active', 'Agonist', 'Partial agonist', 'Antagonist', 'Inverse agonist'), 1, 0),
    Agonist = ifelse(Comment == 'Agonist', 1, 0),
    Partial = ifelse(Comment == 'Partial agonist', 1, 0),
    Antagonist = ifelse(Comment == 'Antagonist', 1, 0),
    Invagonist = ifelse(Comment == 'Inverse agonist', 1, 0)
  )

# Write to file
write.csv(GHSR.activity, './data/parsed/chembl/GHSR.ACTIVITY.csv')

###########################
##  Continuous Datasets  ##
###########################
# Ki
GHSR.ki <- GHSR.list[[5]]
GHSR.ki <- GHSR.ki %>%
  mutate_all(function(x) ifelse(x=='', NA_character_, x)) %>%
  transmute(
    Id = Molecule,
    Id.Type = Id.Type,
    Standard.Value = Standard.Value
  )
write.csv(GHSR.ki, './data/parsed/chembl/GHSR.KI.csv')

# Kb
GHSR.kb <- GHSR.list[[6]]
GHSR.kb <- GHSR.kb %>%
  mutate_all(function(x) ifelse(x=='', NA_character_, x)) %>%
  transmute(
    Id = Molecule,
    Id.Type = Id.Type,
    Standard.Value = Standard.Value
  )
write.csv(GHSR.kb, './data/parsed/chembl/GHSR.KB.csv')

# EC50
GHSR.ec50 <- GHSR.list[[2]]
GHSR.ec50 <- GHSR.ec50 %>%
  mutate_all(function(x) ifelse(x=='', NA_character_, x)) %>%
  transmute(
    Id = Molecule,
    Id.Type = Id.Type,
    Standard.Value = Standard.Value
  )
write.csv(GHSR.ec50, './data/parsed/chembl/GHSR.EC50.csv')

# IC50
GHSR.ic50 <- GHSR.list[[3]]
GHSR.ic50 <- GHSR.ic50 %>%
  mutate_all(function(x) ifelse(x=='', NA_character_, x)) %>%
  transmute(
    Id = Molecule,
    Id.Type = Id.Type,
    Standard.Value = Standard.Value
  )
write.csv(GHSR.ic50, './data/parsed/chembl/GHSR.IC50.csv')




