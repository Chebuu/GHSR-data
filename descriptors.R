library(BioMedR)

## Generally speaking, the max number of descriptors that can be applied to the largest dataset is 6. 
## ## More than 6 descriptors introduces dimensionality problems, wherein m >> n.

set0 <- c(extrDrugEstateComplete)
set1 <- c(set0, extrDrugMACCSComplete, extrDrugKappaShapeIndices)
set2 <- c(set1, extrDrugHybridizationComplete)
set3 <- c(set2, extrDrugBCUT)
set4 <- c(set3, extrDrugKRComplete, extrDrugKierHallSmarts)
set5 <- c(set4, extrDrugIPMolecularLearning)



