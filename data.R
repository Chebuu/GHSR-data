library(dplyr)
library(ChemmineR)

source('util.R')

source('data.chembl.R')
source('data.pubchem.R')
source('data.literature.R')

WRITE_SDF <- TRUE

# Read parsed datasets from each soure and rbind those of the same type
GHSR <- readGHSR.parsed(bind=TRUE)

# Convert all CHEMBLIDs to InChIKey
GHSR <- lapply(GHSR, function(dset){
  dset %>%
    mutate(
      Id = case_when(
        Id.Type == 'chemblid' ~ getINCHIKEY.chembl(Id),
        TRUE ~ Id
      ),
      Id.Type = case_when(
        Id.Type == 'chemblid' ~ 'inchikey',
        TRUE ~ Id.Type
      )
    )
})

# Convert all InChIKeys to CIDs
GHSR <- lapply(GHSR, function(dset){
  dset %>%
    mutate(
      Id = case_when(
        Id.Type == 'inchikey' ~ getCID.inchi(Id),
        TRUE ~ Id
      ),
      Id.Type = case_when(
        Id.Type == 'inchikey' ~ 'cid',
        TRUE ~ Id.Type
      )
    )
})

# Fetch isomeric SMILES for all CIDs
GHSR <- lapply(GHSR, function(dset){
  dset %>%
    mutate(
      smiles = case_when(
        Id.Type == 'cid' ~ unlist(getSMILES.cid(Id, 'isomeric')),
        TRUE ~ NA_character_
      )
    )
})

# Write aggregate datasets to file
for(dsname in names(GHSR)){
  dset <- GHSR[[dsname]]
  dsfile <- sprintf('GHSR.%s.csv', dsname)
  dspath <- normalizePath(paste(AGGREGATE_DIR, dsfile, sep='/'), mustWork=F)
  write.csv(dset, dspath)  
}

# Parse SMILES as SDFset and write to file
if(WRITE_SDF)
  for(dsname in names(GHSR)){
    dset <- GHSR[[dsname]]
    dset.sdf <- smiles2sdf(dset$smiles)
    dset.sdf@ID <- paste(dset$Id.type, dset$Id, sep='=')
    dsfile <- sprintf('GHSR.%s.sdf', dsname)
    dspath <- normalizePath(paste(AGGREGATE_DIR, dsfile, sep='/'), mustWork=F)
    write.SDF(dset.sdf, dspath)
  }
