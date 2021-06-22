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
### load the Mechanistic TSPs
load("./Objs/myTSPs.rda")

### Load data
load("objs/lymphnodeDataSmallBatchCorr.rda")

#################################################################
### Filter genes in pairs

### Available genes
gns <- rownames(combat.mat)

### Common genes
keepGns <-   intersect(as.vector(myTSPs), gns)
str(keepGns)

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

### ### How many
### nrow(myTSPs)
### rownames(myTSPs) <- apply(myTSPs, 1, paste,  collapse=",")


#################################################################
### Best genes

### Ttest
tGns <- t(apply(combat.mat[keepGns, ], 1, function(x,  y) {
	x <- t.test(x~y)
	c(Tstat=x$statistic, Pval=x$p.value)
},  y=grp))

### Reorder
tGns <- tGns[order(tGns[, "Tstat.t"]), ]


### Number of genes
myN <- round(nrow(tGns)/2)

### Select
gnsUP <- tail(rownames(tGns[tGns[, "Tstat.t"] > 0, ]), n=myN)
gnsDW <- head(rownames(tGns[tGns[, "Tstat.t"] < 0, ]), n=myN)

### Selection
keepGns <- c(gnsUP,  gnsDW)

myN <- 2*myN

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]


#################################################################
### Subset matrices

### Get the 2 array GEO Matrices (Grouped together)
ArrayMat <- combat.mat[c("GAPDH",  keepGns), dataset != "TCGA"]
ArrayGroup <- grp[dataset != "TCGA"]

### Get the tcga Mat
tcgaMat <- combat.mat[c("GAPDH",  keepGns), dataset == "TCGA"]
tcgaGroup <- grp[dataset == "TCGA"]

### Normalize to GAPDH
ArrayMat <- sweep(ArrayMat, 2, ArrayMat["GAPDH",],  "-")
ArrayMat <- ArrayMat[ -grep("GAPDH",  rownames(ArrayMat)), ]

### Normalize to GAPDH
tcgaMat <- sweep(tcgaMat, 2, tcgaMat["GAPDH",],  "-")
tcgaMat <- tcgaMat[ -grep("GAPDH",  rownames(tcgaMat)), ]


#################################################################
### Compute scores

### Set Feature number and max k
featN <-  nrow(ArrayMat)
featN

### ### Possible number of TPS
### choose(featN,  2)

### Set k
k <- nrow(myTSPs)

### Train a classifier with the desired maximum number of TSPs
ArrayKTSP <- SWAP.Train.KTSP(ArrayMat, ArrayGroup, krange=k,
			     FilterFunc = NULL,  #SWAP.Filter.Wilcoxon,
			    RestrictedPairs = myTSPs,
			    featureNo=featN,
			    disjoint = FALSE)

### Drop NAs in Train
drop <- !is.na(ArrayKTSP$score)
ArrayKTSP$score <- ArrayKTSP$score[drop]
ArrayKTSP$TSPs <- ArrayKTSP$TSPs[drop, ]
ArrayKTSP$tieVote <- droplevels(ArrayKTSP$tieVote[drop])

### Rename TSPs
trainNms <- apply(ArrayKTSP$TSPs, 1, paste,  collapse=",")
names(ArrayKTSP$score) <- trainNms
rownames(ArrayKTSP$TSPs) <- trainNms
names(ArrayKTSP$tieVote) <- trainNms

### Check consistency with biology
keepTrain <- ArrayKTSP$TSPs[,1] %in% unique(myTSPs[,"GoodGene"])
table(keepTrain)

### Subset
ArrayKTSP$score <- ArrayKTSP$score[keepTrain]
ArrayKTSP$TSPs <- ArrayKTSP$TSPs[keepTrain, ]
ArrayKTSP$tieVote <- droplevels(ArrayKTSP$tieVote[keepTrain])


### Train a classifier with the desired maximum number of TSPs
tcgaKTSP <- SWAP.Train.KTSP(tcgaMat, tcgaGroup, krange=k,
			    FilterFunc = NULL,  #SWAP.Filter.Wilcoxon,
			    RestrictedPairs = myTSPs,
			    featureNo=featN,
			    disjoint = FALSE)

### Drop NAs in Test
drop <- !is.na(tcgaKTSP$score)
tcgaKTSP$score <- tcgaKTSP$score[drop]
tcgaKTSP$TSPs <- tcgaKTSP$TSP[drop, ]
tcgaKTSP$tieVote <- droplevels(tcgaKTSP$tieVote[drop])

### Rename TSPs
testNms <- apply(tcgaKTSP$TSPs, 1, paste,  collapse=",")
names(tcgaKTSP$score) <- testNms
rownames(tcgaKTSP$TSPs) <- testNms
names(tcgaKTSP$tieVote) <- testNms

### Check consistency with biology
keepTest <- tcgaKTSP$TSPs[,1] %in% myTSPs[,"GoodGene"]
table(keepTest)

### Subset
tcgaKTSP$score <- tcgaKTSP$score[keepTest]
tcgaKTSP$TSPs <- tcgaKTSP$TSPs[keepTest, ]
tcgaKTSP$tieVote <- droplevels(tcgaKTSP$tieVote[keepTest])


#################################################################
#################################################################
### Any TSPs in common?

### Intersection
tspNms <- intersect(names(ArrayKTSP$score), names(tcgaKTSP$score))
tspNms

### Which genes
goodGns <- unique(strsplit2(ArrayKTSP$TSPs[tspNms,], split=","))
goodGns

### Check scores: TSP score median
medTSPscore <- apply(cbind(ArrayKTSP$score[tspNms], tcgaKTSP$score[tspNms]), 1, mean)


## Order pairs based on median score
medTSPscore <- medTSPscore[order(medTSPscore, decreasing = T)]

# Take the top 6
tspNms <- names(medTSPscore)[1:6]

### Set TSP score threshold
#tspTHR <- quantile(medTSPscore, 0.91)
#tspTHR

### Plot
#plot(ArrayKTSP$score[tspNms], tcgaKTSP$score[tspNms],  pch=16,
#     col=2+(medTSPscore > tspTHR))

### Which ones?
#tspNms <- names(medTSPscore)[medTSPscore > tspTHR]

### Show filtered
tspNms

### Which genes
goodGns <- unique(strsplit2(ArrayKTSP$TSPs[tspNms,], split=","))
goodGns

### Train
ArrayKTSPFilt <- ArrayKTSP
ArrayKTSPFilt$score <- ArrayKTSPFilt$score[tspNms]
ArrayKTSPFilt$TSPs <- ArrayKTSPFilt$TSPs[tspNms, ]
ArrayKTSPFilt$tieVote <- droplevels(ArrayKTSPFilt$tieVote[tspNms])
ArrayKTSPFilt$name <- paste(nrow(ArrayKTSPFilt$TSPs), "TSPS")  

### Test
tcgaKTSPFilt <- tcgaKTSP
tcgaKTSPFilt$score <- tcgaKTSPFilt$score[tspNms]
tcgaKTSPFilt$TSPs <- tcgaKTSPFilt$TSPs[tspNms, ]
tcgaKTSPFilt$tieVote <- droplevels(tcgaKTSPFilt$tieVote[tspNms])
tcgaKTSPFilt$name <- paste(nrow(tcgaKTSPFilt$TSPs), "TSPS")  

### Check scores
cbind(
	round(ArrayKTSPFilt$score[tspNms], 3), 
	round(tcgaKTSPFilt$score[tspNms],  3)
)

### Assemble for saving
featN <- length(goodGns) 
kParams <- list(featN=featN,  ktsp=length(tspNms))


#################################################################
#################################################################

### Compute the sum and find the best threshold: ALL TRAINING SAMPLES
ktspStatsArray <- SWAP.KTSP.Statistics(
	inputMat = ArrayMat,
	classifier = ArrayKTSPFilt,
	CombineFunc = sum)
summary(ktspStatsArray$statistics)

### Threshold Array
thr_Array <- coords(roc(ArrayGroup, ktspStatsArray$statistics,
                  levels=c("POS", "NEG"),
		  direction="<"), "best")["threshold"]
thr_Array

### Print ROC curve local maximas
auc(roc(ArrayGroup, ktspStatsArray$statistics))

### Coordinates
coords(roc(ArrayGroup, ktspStatsArray$statistics,
	   levels=c("POS", "NEG"),
	   direction="<"), "local maximas")

###### Maximixe sensitivity based on table above
###thr <- 1.5

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
ArrayPrediction <- SWAP.KTSP.Classify(
	ArrayMat,
	ArrayKTSPFilt,
	DecisionFunc = function(x) sum(x) > thr_Array )

### Resubstitution performance in the TRAINING set
confusionMatrix(ArrayPrediction, ArrayGroup, positive="POS")

### Store thrshold
ktspThrMechanistic <- thr_Array


#################################################################
### TESTING
#################################################################

### Compute the sum and find the best threshold
ktspStatsTcga <- SWAP.KTSP.Statistics(
	inputMat = tcgaMat,
	classifier = tcgaKTSPFilt,
	CombineFunc = sum)
summary(ktspStatsTcga$statistics)

### print ROC curve local maximas
auc(roc(tcgaGroup, ktspStatsTcga$statistics))

### Coordinates
coords(roc(tcgaGroup, ktspStatsTcga$statistics,
	   levels=c("POS", "NEG"),
	   direction="<"), "local maximas")

## threshold tcga
thr_TCGA <- coords(roc(tcgaGroup, ktspStatsTcga$statistics,
           levels=c("POS", "NEG"),
           direction="<"), "best")["threshold"]
thr_TCGA
### Get prediction based on best threshold from ROC curve
### Note the use of ">"
tcgaPrediction <- SWAP.KTSP.Classify(
	tcgaMat,
	tcgaKTSPFilt,
	DecisionFunc = function(x) sum(x) > thr_TCGA )

### Resubstitution performance in the TRAINING set
confusionMatrix(tcgaPrediction, tcgaGroup, positive="POS")


## Rename the 2 classifiers
ArrayKTSP <- ArrayKTSPFilt
tcgaKTSP <- tcgaKTSPFilt
#################################################################
### Save
save(ArrayKTSP, tcgaKTSP, ArrayMat, tcgaMat, ArrayGroup, tcgaGroup, file = "./Objs/FinalClassifiers.rda")
#################################################################
### Session information
sessionInfo()

### Clean working space and quit
rm(list=ls())
q("no")
