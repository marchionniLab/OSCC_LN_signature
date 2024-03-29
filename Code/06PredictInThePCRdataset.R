
rm(list = ls())

library(switchBox)
library(caret)
library(dcurves)
library(readxl)

###########################################################################################
load("./Objs/FinalClassifiers.rda")

## Load the RT-PCR data
pcrData <- read.delim("./Data/RT_PCR/deltaCt.txt")

pcrData2 <- read_excel("data/RT_PCR/Table1.xlsx")

# Get the node status
pcrGroup <- pcrData$NodeStatus
table(pcrGroup)
# Rename the groups from 0,1 to neg,pos
pcrGroup[pcrGroup == 0] <- "NEG"
pcrGroup[pcrGroup == 1] <- "POS"

# Get the matrix
pcrMat <- pcrData[, -14]
rownames(pcrMat) <- pcrMat$SampleName
pcrMat$SampleName <- NULL

# Transpose the matrix
pcrMat <- t(pcrMat)

# One gene is misspelled, rename it
rownames(pcrMat)[rownames(pcrMat) == "TFGB2"] <- "TGFB2"

### Compute the sum and find the best threshold: ALL TRAINING SAMPLES
ktspStatsPCR <- SWAP.KTSP.Statistics(
  inputMat = pcrMat,
  classifier = tcgaKTSP,
  CombineFunc = sum)

summary(ktspStatsPCR$statistics)


### Print ROC curve local maximas
auc(roc(pcrGroup, ktspStatsPCR$statistics))

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
pcrPrediction <- SWAP.KTSP.Classify(
  pcrMat,
  tcgaKTSP,
  DecisionFunc = function(x) sum(x) > 2.5 )

pcrGroup <- factor(pcrGroup, levels = c("POS", "NEG"))
table(pcrGroup)

### Resubstitution performance in the TRAINING set
PCRperf <- confusionMatrix(pcrPrediction, pcrGroup, positive="POS")
PCRperf_Classes <- as.matrix(PCRperf, what = "classes")
PCRperf_Xtabs <- as.matrix(PCRperf, what = "xtabs")

write.csv(PCRperf_Xtabs, file = "./objs/Performance_All/PCR_Xtabs.csv")

write.csv(PCRperf_Classes, file = "./objs/Performance_All/PCR_Performance.csv")






