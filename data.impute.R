getRefsByColumn <- function(GHSR){
  reference <- do.call(
    rbind,
    lapply(names(GHSR), function(sbjname){
      dsbj <- GHSR[[sbjname]]
      do.call(
        rbind,
        lapply(names(GHSR), function(refname){
          dref <- GHSR[[refname]]
          left_join(dsbj, dref, by='Id', suffix=c('.S', '.R'))[
            c('Id','Standard.Value.S', 'Standard.Value.R')
          ] %>%
            mutate(
              Subject=sbjname,
              Reference = refname
            )
        })
      )
    })
  )
}

RFBC <- getRefsByColumn(GHSR[c('EC50', 'IC50', 'KI', 'KB')])
RFBC
# This reveals a lot of malformatted data. Not sure what the cause is.
RFBC %>%
  filter(Subject == Reference & Standard.Value.S != Standard.Value.R)


RFBC.xrf.incomplete <- RFBC %>%
  filter(Subject != Reference & xor(is.na(Standard.Value.S), is.na(Standard.Value.R)))

RFBC.xrf.complete <- RFBC %>%
  filter(Subject != Reference & !is.na(Standard.Value.S) & !is.na(Standard.Value.R))

tst <- RFBC.xrf.complete %>%
  filter(Subject == 'EC50' & Reference == 'IC50')

which(GHSR$EC50$Id == 'CHEMBL500468')
plot(log(tst$Standard.Value.S), log(tst$Standard.Value.R))

isx <- intersect(GHSR$EC50$Id, GHSR$IC50$Id)

data.frame(GHSR$EC50$Standard.Value[GHSR$EC50$Id == isx], GHSR$IC50$Standard.Value[GHSR$IC50$Id == isx])
GHSR$EC50$Standard.Value[GHSR$EC50$Id == isx]
GHSR$IC50$Standard.Value[GHSR$IC50$Id == isx]


