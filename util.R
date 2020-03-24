library(dplyr)
library(httr)
library(RCurl)
library(XML)

AGGREGATE_DIR <- 'data/aggregate'
AGGREGATE_DNAMES <- c('INACTIVE','ACTIVITY', 'EC50', 'IC50', 'KI', 'KB')

PARSED_DIR <- 'data/parsed'
PARSED_SOURCES <- c('chembl','pubchem','literature')
PARSED_DNAMES <- c('INACTIVE','ACTIVITY', 'EC50', 'IC50', 'KI', 'KB')

readGHSR.parsed <- function(bind=TRUE, parsed_dir=PARSED_DIR, parsed_sources=PARSED_SOURCES, parsed_dnames=PARSED_DNAMES){
  #' Read all parsed datasets and return them as items in a named list. Optionally rbind all datasets in each subdir listed in PARSED_SOURCES.
  setNames(lapply(parsed_dnames, function(dsname){
    dsets <- setNames(lapply(parsed_sources, function(dsdir){
      dsfile <- sprintf('GHSR.%s.csv', dsname)
      dspath <- normalizePath(paste(parsed_dir, paste(dsdir, dsfile, sep='/'), sep='/'))
      read.table(dspath, header=TRUE, sep=',', row.names=NULL, strip.white=TRUE, stringsAsFactors=FALSE)
    }), parsed_sources)
    if(bind) dsets <- do.call(rbind, dsets)  
    return(dsets)
  }), parsed_dnames)
}

readGHSR.aggregate <- function(aggregate_dir=AGGREGATE_DIR, aggregate_dnames=AGGREGATE_DNAMES){
  #' Read all aggregate datasets and return them as items in a named list.
  setNames(lapply(aggregate_dnames, function(adname){
    adfile <- sprintf('GHSR.%s.csv', adname)
    adpath <- normalizePath(paste(aggregate_dir, adfile, sep='/'))
    read.table(adpath, header=TRUE, sep=',', row.names=NULL, strip.white=TRUE, stringsAsFactors=FALSE)
  }), aggregate_dnames)
}

readSDF.aggregate <- function(aggregate_dir=AGGREGATE_DIR, aggregate_dnames=AGGREGATE_DNAMES){
  #' Read SDFsets for each aggregate dataset
  setNames(lapply(aggregate_dnames, function(adname){
    adfile <- sprintf('GHSR.%s.sdf', adname)
    adpath <- normalizePath(paste(aggregate_dir, adfile, sep='/'))
    read.SDFset(adpath, skipErrors=FALSE)
  }), aggregate_dnames)
}

blank_if <- function(x, y){ ifelse(x==y, '', x) }

blank_if_na <- function(x){ ifelse(is.na(x), '', x) }

replace_if_na <- function(x, y){ ifelse(is.na(x), y, x) }

splitListForURL <- function(x, x.char.sep, n.char.max=1800, na.fill=''){
  #' Separate the items in \code{x} into a list of csv (or other delimiter) strings, where each string consists of a defined number of total characters (including sep characters). 
  #' @param x A character vector of items to concatenate with \code{x.char.sep} then split into chunks of \code{n.char.max} characters. 
  #' @param x.char.sep The characters used to separate x during concatenation.
  #' @param n.char.max The maximum number of characters (including sep characters) to allow in a given csv (or other sep) string.
  #' @examples \dontrun{ 
  #' chunks <- splitListForURL(seq_len(2000), x.char.sep=',', n.char.max=1000)
  #' lapply(chunks, length) 
  #' }
  x[is.na(x)] <- na.fill
  x.p0 <- paste0(x, sep=',')
  x.p0.nchar <- nchar(x.p0)
  x.chk <- c()
  x.chk.idx <- 1
  x.chk.sum <- 0
  chunks <- c()
  for(i in seq(from=2, to=length(x.p0.nchar))){
    x.chk <- paste(x.chk, x.p0[[i-1]], sep='')
    x.chk.sum <- x.chk.sum + x.p0.nchar[[i-1]]
    if(x.chk.sum + x.p0.nchar[[i]] > n.char.max){
      chunks <- c(chunks, x.chk)
      x.chk <- c()
      x.chk.sum <- 0
      x.chk.idx <- i
    }
  }
  x.chk <- paste(x.chk, x.p0[[length(x.p0)]], sep='')
  chunks <- c(chunks, x.chk)
  return(chunks)
}

getINCHIKEY.chembl <- function(chemblids){
  # NOTE: This grabs the first matching InChiKey from the response. I'm not sure how EBI encodes CHEMBL IDs, but I'm guessing there are multiple InChIKeys for a given CHEMBL ID.
  #' Fetch InChIKeys for one or many CHEMBLIDs. 
  #' @value A named list of InChIKeys, where names are the corresponding CHEMBLID for each InChIKey. Indices of the returned vector match those of the argument. If a corresponding InChIKey for a given CHEMBLID cannot be found, its value will be NA_character_ in the returned list.
  #' @param chemblids A vector of CHEMBLIDs
  message(sprintf('Fetching InChIKeys (%s)...', length(chemblids)))
  schema <- 'https://www.ebi.ac.uk/chembl/api/data/molecule?limit=1000&molecule_chembl_id__in='
  chemblids.split <- splitListForURL(chemblids)
  mappings <- do.call(
    rbind,
    lapply(chemblids.split, function(query){
      url <- paste(schema, query, sep='')
      res <- getURL(url)
      doc <- xmlParse(res)
      root <- xmlRoot(doc)
      molecules <- xmlElementsByTagName(root, 'molecules', recursive=FALSE)[[1]]
      molecule <- xmlElementsByTagName(molecules, 'molecule', recursive=FALSE)
      split.map <- lapply(molecule, function(mol){
        mol_cbl_id <- xmlElementsByTagName(mol, 'molecule_chembl_id', recursive=FALSE)[[1]]
        mol_std_inc <- xmlElementsByTagName(mol, 'standard_inchi_key', recursive=TRUE)[[1]]
        inchikey <- xmlValue(mol_std_inc)
        chemblid <- xmlValue(mol_cbl_id)
        return(
          data.frame(
            chemblid=chemblid,
            inchikey=inchikey
          )
        )
      })
      do.call(rbind,split.map)
    })
  )
  message('Mapping values...')
  sapply(chemblids, function(id){
    midx <- mappings$chemblid == id
    if(all(!midx) || any(is.na(midx))) return(NA_character_)
    inchi.map <- mappings$inchikey[midx][[1]]
    as.character(inchi.map)
  })
}

getCID.inchi <- function(inchis){
  # NOTE: Encoding csv values in POST request body returns an error stating CIDs cannot be found, so a GET request is used here instead. 
  # NOTE: Its tempting to use the "cids" <operation specification> endpoint, but the server returns CIDs out of order and without any reference to the corresponding InChIKey in the query, so there's no way to map a CID returned by the server to its coresponding InChIKey. The solution is to return the entire compound record, which contains both InChIKey and CID for each compound.
  # NOTE: Requesting whole compound records for a long list of compounds yields large response bodies that exceed my memory and crash the script. To limit response size I've shortened the list in each query by setting n.char.max=600. This comes with a drastic increase in execution time. A method to optimize the list size by requesting the size of the response apriori is needed.
  #' Fetch CIDs for one or many InChIKeys.
  #' @value A named list of CIDs, where names are the corresponding InChIKey for each CID. Indices of the returned vector match those of the argument. If a corresponding CID for a given InChIKey cannot be found, its value will be NA in the returned list.
  #' @param chemblids A vector of InChIKeys
  #' @example \dontrun{
  #' getCID.inchi(c("CZYZPYHWBJNCKH-UHFFFAOYSA-N","FEERXCVFQAXCNL-UHFFFAOYSA-N", "33705", 62754))
  #' }
  message(sprintf('Fetching CIDs (%s)...', length(inchis)))
  idlist.split <- splitListForURL(inchis, x.char.sep=',', n.char.max=900)
  mappings <- as.data.frame(
    do.call(
      rbind,
      lapply(idlist.split, function(idls){
        idlexpand <- strsplit(idls,',')[[1]]
        idlexpand <- idlexpand[idlexpand!='']
        noneInChI <- all(suppressWarnings(!is.na(as.numeric(idlexpand))))
        if(noneInChI) return(do.call(rbind,lapply(idlexpand, function(cid){c(inchi=NA_character_, cid=cid)})))
        get.url <- sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/record/XML', idls)
        res <- GET(get.url)
        cnt <- rawToChar(res$content)
        doc <- xmlParse(cnt)
        PC_Compounds <- xmlElementsByTagName(doc, 'PC-Compounds', recursive=FALSE)[[1]]
        PC_Compound <- xmlElementsByTagName(PC_Compounds, 'PC-Compound', recursive=FALSE)
        do.call(
          rbind,
          lapply(PC_Compound, function(node){
            pcid <- xmlElementsByTagName(node, 'PC-Compound_id', recursive=FALSE)[[1]]
            pct <- xmlElementsByTagName(pcid, 'PC-CompoundType', recursive=FALSE)[[1]]
            pcti <- xmlElementsByTagName(pct, 'PC-CompoundType_id', recursive=FALSE)[[1]]
            pctic <- xmlElementsByTagName(pcti, 'PC-CompoundType_id_cid', recursive=FALSE)[[1]]
            cid <- xmlValue(pctic)
            inchi <- NULL
            pcp <- xmlElementsByTagName(node, 'PC-Compound_props', recursive=FALSE)[[1]]
            pcid <- xmlElementsByTagName(pcp, 'PC-InfoData', recursive=FALSE)
            for(p in pcid){
              pcul <- xmlElementsByTagName(p, 'PC-Urn_label', recursive=TRUE)[[1]]
              if(xmlValue(pcul) == 'InChIKey'){
                pcifd <- pcul %>% xmlParent %>% xmlParent %>% xmlParent
                pcsvl <- xmlElementsByTagName(pcifd, 'PC-InfoData_value_sval', recursive=TRUE)
                inchi <- xmlValue(pcsvl[[1]])  
              }
            }
            return(
              c(
                inchi=inchi,
                cid=cid
              )
            )
          }) 
        )
      }) 
    ), stringsAsFactors=F
  )
  message('Mapping values...')
  sapply(inchis, function(inc){
    midx <- mappings$inchi == inc
    if(all(!midx) || any(is.na(midx))) return(NA_character_)
    cid.map <- mappings$cid[midx][[1]]
    return(cid.map)
  })
}

getSMILES.cid <- function(cids){
  #' Fetch isomeric SMILES strings for one or many CIDs.
  #' @value A named list of isomeric SMILES strings, where names are the corresponding CID for each SMILES string. Indices of the returned vector match those of the argument. If a corresponding isomeric SMILES string for a given CID cannot be found, its value will be NA_character_ in the returned list.
  #' @param chemblids A vector of CHEMBLIDs
  #' @example \dontrun{
  #' getSMILES.cid(c("9890863", "9914145", "10874830"))
  #' }
  message(sprintf('Fetching SMILES (%s)...', length(cids)))
  idlist.split <- splitListForURL(inchis, x.char.sep=',', n.char.max=900)
  mappings <- as.data.frame(
    do.call(
      rbind,
      lapply(idlist.split, function(idls){
        get.url <- sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/record/XML', idls)
        res <- GET(get.url)
        cnt <- rawToChar(res$content)
        doc <- xmlParse(cnt)
        PC_Compounds <- xmlElementsByTagName(doc, 'PC-Compounds', recursive=FALSE)[[1]]
        PC_Compound <- xmlElementsByTagName(PC_Compounds, 'PC-Compound', recursive=FALSE)
        do.call(
          rbind,
          lapply(PC_Compound, function(node){
            pcid <- xmlElementsByTagName(node, 'PC-Compound_id', recursive=FALSE)[[1]]
            pct <- xmlElementsByTagName(pcid, 'PC-CompoundType', recursive=FALSE)[[1]]
            pcti <- xmlElementsByTagName(pct, 'PC-CompoundType_id', recursive=FALSE)[[1]]
            pctic <- xmlElementsByTagName(pcti, 'PC-CompoundType_id_cid', recursive=FALSE)[[1]]
            cid <- xmlValue(pctic)
            inchi <- NULL
            smiles <- NULL
            pcp <- xmlElementsByTagName(node, 'PC-Compound_props', recursive=FALSE)[[1]]
            pcid <- xmlElementsByTagName(pcp, 'PC-InfoData', recursive=FALSE)
            for(p in pcid){
              pcun <- tryCatch((xmlElementsByTagName(p, 'PC-Urn_name', recursive=TRUE)[[1]]), error=function(e){return(FALSE)})
              if(isFALSE(pcun)) return(c(cid=cid, smiles=smiles))
              if(xmlValue(pcun) == 'Isomeric'){
                pcifd <- pcun %>% xmlParent %>% xmlParent %>% xmlParent
                pcsvl <- xmlElementsByTagName(pcifd, 'PC-InfoData_value_sval', recursive=TRUE)
                smiles <- xmlValue(pcsvl[[1]])
              }
            }
            return(
              c(
                cid=cid,
                smiles=smiles
              )
            )
          }) 
        )
      }) 
    )
  )
  message('Mapping values...')
  sapply(cids, function(id){
    midx <- mappings$cid == id
    if(sum(midx) < 1 || any(is.na(midx))) return(NA_character_)
    smiles.map <- mappings$smiles[midx][[1]]
    return(smiles.map)
  })
}

getConformers <- function(){
  # A list of diverse order conformer IDs can be obtained from CID. Valid output formats are XML, JSON(P), ASNT/B, and TXT (limited to a single CID):
  #   
  #   https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/conformers/XML
  # 
  # Individual conformer records – either computed 3D coordinates for compounds or deposited/experimental 3D coordinates for some substances – can be retrieved by conformer ID:
  #   
  #   https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/000008C400000001/SDF
}
