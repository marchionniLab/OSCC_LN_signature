#################################################################
### TCGA head and neck data for LN statu kTSP prediction and validation
### Luigi Marchionni
### Collaboration with Yasmen Ghantous (Sidransky lab)
### Unrestricted kTSP


#################################################################
### Clean
rm(list=ls())

### Setwd
setwd("/Volumes/Macintosh/Research/Projects/HNSCC")


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


### #################################################################
### ### Source for filtering genes in switchBox
### source("~/Rvers/Rfuncs/switchBoxFilters.R")

#################################################################
### Load data
# load("./Objs/lymphnodeDataSmall.rda")
# 
# ### Normalize to GAPDH
# trainMat <- sweep(trainMat, 2, trainMat["GAPDH",],  "-")
# trainMat <- trainMat[ -grep("GAPDH",  rownames(trainMat)), ]
# 
# ### Normalize to GAPDH
# testMat <- sweep(testMat, 2, testMat["GAPDH",],  "-")
# testMat <- testMat[ -grep("GAPDH",  rownames(testMat)), ]


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
#### Diff genes
# myN <- round(nrow(trainMat)/1)
# myN
# 
# trainGns <- SWAP.Filter.Wilcoxon(trainGrp, trainMat, featureNo=myN)
# testGns <- SWAP.Filter.Wilcoxon(testGrp, testMat, featureNo=myN)
# 
# ### ### Any in common?
# gns <- intersect(trainGns, testGns)
# ### gns <- union(trainGns, testGns)
# str(gns)
# 
# 
# #################################################################
# ### Filter genes in pairs
# 
# ### Common genes
# keepGns <-   intersect(as.vector(myTSPs), gns)
# str(keepGns)
# 
# ### For the TSP
# myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]
# 
# ### How many
# nrow(myTSPs)
# rownames(myTSPs) <- apply(myTSPs, 1, paste,  collapse=",")
# 
# ### Save Resticted TSPs
# save(myTSPs, file="./objs/HallmarksComb/allTSPs.rda")
# 
# ### Filter datasets
# trainMat <- trainMat[ keepGns , ]
# testMat <- testMat[ keepGns , ]
# 
# 
# #################################################################
# ### Compute scores
# 
# ### Set Feature number and max k
# featN <-  nrow(trainMat)
# featN
# 
# ### ### Possible number of TPS
# ### choose(featN,  2)
# 
# ### Set k
# k <- nrow(myTSPs)
# 
# ### Train a classifier with the desired maximum number of TSPs
# trainTSPs <- SWAP.Train.KTSP(trainMat, trainGrp, krange=k,
# 			    FilterFunc = SWAP.Filter.Wilcoxon,
# 			    RestrictedPairs = myTSPs,
# 			    featureNo=featN,
# 			    disjoint = FALSE)
# 
# ### Drop NAs in Train
# drop <- !is.na(trainTSPs$score)
# trainTSPs$score <- trainTSPs$score[drop]
# trainTSPs$TSPs <- trainTSPs$TSPs[drop, ]
# trainTSPs$tieVote <- droplevels(trainTSPs$tieVote[drop])
# 
# ### Rename TSPs
# trainNms <- apply(trainTSPs$TSPs, 1, paste,  collapse=",")
# names(trainTSPs$score) <- trainNms
# rownames(trainTSPs$TSPs) <- trainNms
# names(trainTSPs$tieVote) <- trainNms
# 
# ### Check consistency with biology
# keepTrain <- trainTSPs$TSPs[,1] %in% unique(myTSPs[,"GoodGene"])
# table(keepTrain)
# 
# ### Subset
# trainTSPs$score <- trainTSPs$score[keepTrain]
# trainTSPs$TSPs <- trainTSPs$TSPs[keepTrain, ]
# trainTSPs$tieVote <- droplevels(trainTSPs$tieVote[keepTrain])
# 
# 
# ### Train a classifier with the desired maximum number of TSPs
# testTSPs <- SWAP.Train.KTSP(testMat, testGrp, krange=nrow(myTSPs),
# 			    FilterFunc = SWAP.Filter.Wilcoxon,
# 			    RestrictedPairs = myTSPs,
# 			    featureNo=featN,
# 			    disjoint = FALSE)
# 
# ### Drop NAs in Test
# drop <- !is.na(testTSPs$score)
# testTSPs$score <- testTSPs$score[drop]
# testTSPs$TSPs <- testTSPs$TSP[drop, ]
# testTSPs$tieVote <- droplevels(testTSPs$tieVote[drop])
# 
# ### Rename TSPs
# testNms <- apply(testTSPs$TSPs, 1, paste,  collapse=",")
# names(testTSPs$score) <- testNms
# rownames(testTSPs$TSPs) <- testNms
# names(testTSPs$tieVote) <- testNms
# 
# ### Check consistency with biology
# keepTest <- testTSPs$TSPs[,1] %in% myTSPs[,"GoodGene"]
# table(keepTest)
# 
# ### Subset
# testTSPs$score <- testTSPs$score[keepTest]
# testTSPs$TSPs <- testTSPs$TSPs[keepTest, ]
# testTSPs$tieVote <- droplevels(testTSPs$tieVote[keepTest])
# 
# 
# #################################################################
# #################################################################
# ### Any TSPs in common?
# 
# ### Intersection
# tspNms <- intersect(names(trainTSPs$score), names(testTSPs$score))
# tspNms
# 
# ### Which genes
# goodGns <- unique(strsplit2(trainTSPs$TSPs[tspNms,], split=","))
# goodGns
# 
# ### Check scores: TSP score median
# medTSPscore <- apply(cbind(trainTSPs$score[tspNms], testTSPs$score[tspNms]), 1, median)
# 
# ### Set TSP score threshold
# tspTHR <- quantile(medTSPscore, 0.75)
# tspTHR
# 
# ### Plot
# plot(trainTSPs$score[tspNms], testTSPs$score[tspNms],  pch=16,
#      col=2+(medTSPscore > tspTHR))
# 
# ### Which ones?
# tspNms <- names(medTSPscore)[medTSPscore > tspTHR]
# 
# ### Show filtered
# tspNms
# 
# ### Which genes
# goodGns <- unique(strsplit2(trainTSPs$TSPs[tspNms,], split=","))
# goodGns
# 
# ### Train
# trainTSPsFilt <- trainTSPs
# trainTSPsFilt$score <- trainTSPsFilt$score[tspNms]
# trainTSPsFilt$TSPs <- trainTSPsFilt$TSPs[tspNms, ]
# trainTSPsFilt$tieVote <- droplevels(trainTSPsFilt$tieVote[tspNms])
# trainTSPsFilt$name <- paste(nrow(trainTSPsFilt$TSPs), "TSPS")  
# 
# ### Test
# testTSPsFilt <- testTSPs
# testTSPsFilt$score <- testTSPsFilt$score[tspNms]
# testTSPsFilt$TSPs <- testTSPsFilt$TSPs[tspNms, ]
# testTSPsFilt$tieVote <- droplevels(testTSPsFilt$tieVote[tspNms])
# testTSPsFilt$name <- paste(nrow(testTSPsFilt$TSPs), "TSPS")  
# 
# ### Check scores
# cbind(
# 	round(trainTSPsFilt$score[tspNms], 3), 
# 	round(testTSPsFilt$score[tspNms],  3)
# )
# 
# ### Assemble for saving
# featN <- length(goodGns) 
# kParams <- list(featN=featN,  ktsp=length(tspNms))
# 
# 
# #################################################################
# #################################################################
# 
# ### Compute the sum and find the best threshold: ALL TRAINING SAMPLES
# ktspStatsTrain <- SWAP.KTSP.Statistics(
# 	inputMat = trainMat,
# 	classifier = trainTSPsFilt,
# 	CombineFunc = sum)
# summary(ktspStatsTrain$statistics)
# 
# ### Threshold
# thr <- coords(roc(trainGrp, ktspStatsTrain$statistics,
#                   levels=c("POS", "NEG"),
# 		  direction="<"), "best")["threshold"]
# thr
# 
# ### Print ROC curve local maximas
# auc(roc(trainGrp, ktspStatsTrain$statistics))
# 
# ### Coordinates
# coords(roc(trainGrp, ktspStatsTrain$statistics,
# 	   levels=c("POS", "NEG"),
# 	   direction="<"), "local maximas")
# 
# ### Maximixe sensitivity based on table above
# thr <- 1.5
# 
# ### Get prediction based on best threshold from ROC curve
# ### Note the use of ">"
# trainPrediction <- SWAP.KTSP.Classify(
# 	trainMat,
# 	trainTSPsFilt,
# 	DecisionFunc = function(x) sum(x) > thr )
# 
# ### Resubstitution performance in the TRAINING set
# confusionMatrix(trainPrediction, trainGrp, positive="NEG")
# 
# ### Store thrshold
# ktspThrMechanistic <- thr
# 
# 
# #################################################################
# ### TESTING
# #################################################################
# 
# ### Compute the sum and find the best threshold
# ktspStatsTest <- SWAP.KTSP.Statistics(
# 	inputMat = testMat,
# 	classifier = testTSPsFilt,
# 	CombineFunc = sum)
# summary(ktspStatsTest$statistics)
# 
# ### print ROC curve local maximas
# auc(roc(testGrp, ktspStatsTest$statistics))
# 
# ### Coordinates
# coords(roc(testGrp, ktspStatsTest$statistics,
# 	   levels=c("POS", "NEG"),
# 	   direction="<"), "local maximas")
# 
# ### Get prediction based on best threshold from ROC curve
# ### Note the use of ">"
# testPrediction <- SWAP.KTSP.Classify(
# 	testMat,
# 	testTSPsFilt,
# 	DecisionFunc = function(x) sum(x) > thr )
# 
# ### Resubstitution performance in the TRAINING set
# confusionMatrix(testPrediction, testGrp, positive="NEG")
# 
# 
# 
# #################################################################
# ### Save
# save(list=c(ls(pattern="^ktsp"), "trainTSPsFilt", "testTSPsFilt"),
#      file=paste("objs/HallmarksComb/dge",  myN, ".good",  featN,  ".mech.ktspPredictor.rda",  sep=""))
# save(kParams, myN, file=paste("objs/HallmarksComb/dge", myN,  ".good",  featN,  ".mech.kParams.rda",  sep=""))
# 
# 
# #################################################################
# ### Session information
# sessionInfo()
# 
# ### Clean working space and quit
# rm(list=ls())
# q("no")
