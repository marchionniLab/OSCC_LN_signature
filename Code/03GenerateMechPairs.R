#################################################################
### TCGA head and neck data for LN statu kTSP prediction and validation
### Luigi Marchionni
### Collaboration with Yasmen Ghantous (Sidransky lab)
### Unrestricted kTSP


#################################################################
### Clean
rm(list=ls())

### Setwd
#setwd("/Volumes/Macintosh/Research/Projects/HNSCC")


### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)


#################################################################
### Create output directory for figures
dir.create("./Figs/HallmarksComb", showWarnings = FALSE)
dir.create("./Objs/HallmarksComb", showWarnings = FALSE)
dir.create("./Text/HallmarksComb", showWarnings = FALSE)


### #################################################################
### Load pairs to create the matrix of restricted TSPs: HallmarksComb genes

### ### Read HallmarksComb gene pairs
nms <- c("./Data/Pairs/CancerHallmarks/Angiogenesis/objs/angiogenesisPairs.v6.1.rda",
	 "./Data/Pairs/CancerHallmarks/Apoptosis/objs/apoptosisPairs.v6.1.rda",
	 "./Data/Pairs/CancerHallmarks/CellGrowth/objs/growthPairs.v6.1.rda",
	 "./Data/Pairs/CancerHallmarks/ImmunoSurv/objs/immuPairs.v6.1.rda",
	 "./Data/Pairs/CancerHallmarks/Invasion/objs/emtPairs.v6.1.rda",
	 "./Data/Pairs/TumorAssociatedGenes/TAGs.rda")
names(nms) <- gsub("Pairs",  "", gsub("\\..+",  "",  gsub(".+/",  "",  nms)))
	
### Load
allPairs <- lapply(nms,  function(x) get(load(x)))

### Reorder columns
allPairs <- mapply(x=allPairs,  type= names(nms),  FUN=function(x, type) {
	## Reorder: Good gene first
	if (type == "angiogenesis") x <- cbind(x[,  c("antiAngiogenesis", "proAngiogenesis")], type)
	else if (type == "apoptosis") x <- cbind(x[,  c("proApoptotic", "antiApoptotic")], type)
	else if (type == "growth") x <- cbind(x[,  c("antiGrowth", "proGrowth")], type)
	else if (type == "immu") x <- cbind(x[,  c("ImmSurvON", "ImmSurvOFF")], type)
	else if (type == "emt") x <- cbind(x[,  c("antiInvasion", "proInvasion")], type)
	else if (type == "TAGs") x <- cbind(x[,  c("TSG", "Oncogene")], type)
	x
},  SIMPLIFY=FALSE)

###  Rename
allPairs <- lapply(allPairs,  FUN=function(x) {
	colnames(x) <- c("GoodGene",  "BadGene",  "Type")
	x
})

### Combine
allTSPs <- Reduce("rbind",  allPairs)
dim(allTSPs)
save(allTSPs, file = "./Objs/allTSPs_Annot.rda")

### prepare object
myTSPs <- unique(allTSPs[, 1:2])
dim(myTSPs)
length(unique(as.vector(myTSPs)))

save(myTSPs, file = "./Objs/myTSPs.rda")
#################################################################
