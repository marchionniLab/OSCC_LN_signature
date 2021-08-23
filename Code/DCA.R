


rm(list = ls())

library(switchBox)
library(caret)
library(dcurves)
library(readxl)
library(epiR)

####################################
## Load the classifiers
load("./Objs/FinalClassifiers.rda")

## Load the data
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

confusionMatrix(pcrPrediction, pcrGroup, positive="POS")
data <- as.table(matrix(c(19,3,0,13), nrow = 2, byrow = TRUE))
rval <- epi.tests(data, conf.level = 0.95)
print(rval)

######
dat_PCR = data.frame('LN_status' = pcrGroup, 'signature' = pcrPrediction)


dat_PCR$LN_status <- as.numeric(dat_PCR$LN_status)
dat_PCR$LN_status[dat_PCR$LN_status == '2'] <- 0
dat_PCR$signature <- as.numeric(dat_PCR$signature)
dat_PCR$signature[dat_PCR$signature == '2'] <- 0

## Plot the DCA
pdf("./Figs/DCA/DCA_PCR.pdf", width = 10, height = 8, onefile = F)
dca(LN_status ~ signature, 
    data = dat_PCR,
    thresholds = seq(0, 0.30, by = 0.01)
) %>%
  plot(smooth = TRUE)
dev.off()

## net_intervention_avoided
pdf("./Figs/DCA/NetInt_PCR.pdf", width = 10, height = 8, onefile = F)
dca(LN_status ~ signature, 
    data = dat_PCR,
    thresholds = seq(0, 0.30, by = 0.01), as_probability = "signature"
) %>%
  net_intervention_avoided(nper = 100) %>%
  #standardized_net_benefit() %>%
  plot(smooth = TRUE, type = 'net_intervention_avoided')

dev.off()

#############################################################################
## Array

# Get the prediction stats
ktspStatsArray <- SWAP.KTSP.Statistics(
  inputMat = ArrayMat,
  classifier = ArrayKTSP,
  CombineFunc = sum)

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
ArrayPrediction <- SWAP.KTSP.Classify(
  ArrayMat,
  ArrayKTSP,
  DecisionFunc = function(x) sum(x) > 2.5 )

### Print ROC curve local maximas
auc(roc(ArrayGroup, ktspStatsArray$statistics))
confusionMatrix(ArrayPrediction, ArrayGroup, positive="POS")

data <- as.table(matrix(c(23,5,9,43), nrow = 2, byrow = TRUE))
rval <- epi.tests(data, conf.level = 0.95)
print(rval)

dat_Array = data.frame('LN_status' = ArrayGroup, 'signature' = ArrayPrediction)


dat_Array$LN_status <- as.numeric(dat_Array$LN_status)
dat_Array$LN_status[dat_Array$LN_status == '2'] <- 0
dat_Array$signature <- as.numeric(dat_Array$signature)
dat_Array$signature[dat_Array$signature == '2'] <- 0

## Plot the DCA
pdf("./Figs/DCA/DCA_Array.pdf", width = 10, height = 8, onefile = F)
dca(LN_status ~ signature, 
    data = dat_Array,
    thresholds = seq(0, 0.30, by = 0.01)
) %>%
  plot(smooth = TRUE)
dev.off()

## net_intervention_avoided
pdf("./Figs/DCA/NetInt_Array.pdf", width = 10, height = 8, onefile = F)
dca(LN_status ~ signature, 
    data = dat_Array,
    thresholds = seq(0, 0.3, by = 0.01), as_probability = "signature"
) %>%
  net_intervention_avoided(nper = 100) %>%
  #standardized_net_benefit() %>%
  plot(smooth = TRUE, type = 'net_intervention_avoided')

dev.off()


############################################################################
## TCGA

# Get the prediction stats
ktspStatsTCGA <- SWAP.KTSP.Statistics(
  inputMat = tcgaMat,
  classifier = ArrayKTSP,
  CombineFunc = sum)

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
TCGAPrediction <- SWAP.KTSP.Classify(
  tcgaMat,
  ArrayKTSP,
  DecisionFunc = function(x) sum(x) > 2.5 )

### Print ROC curve local maximas
auc(roc(tcgaGroup, ktspStatsTCGA$statistics))

confusionMatrix(TCGAPrediction, tcgaGroup, positive="POS")
data <- as.table(matrix(c(21,19,6,63), nrow = 2, byrow = TRUE))
rval <- epi.tests(data, conf.level = 0.95)
print(rval)

dat_TCGA = data.frame('LN_status' = tcgaGroup, 'signature' = TCGAPrediction)

dat_TCGA$LN_status <- as.numeric(dat_TCGA$LN_status)
dat_TCGA$LN_status[dat_TCGA$LN_status == '2'] <- 0
dat_TCGA$signature <- as.numeric(dat_TCGA$signature)
dat_TCGA$signature[dat_TCGA$signature == '2'] <- 0

## Plot the DCA
pdf("./Figs/DCA/DCA_TCGA.pdf", width = 10, height = 8, onefile = F)
dca(LN_status ~ signature, 
    data = dat_TCGA,
    thresholds = seq(0, 0.30, by = 0.01)
) %>%
  plot(smooth = TRUE)
dev.off()

## net_intervention_avoided
pdf("./Figs/DCA/NetInt_TCGA.pdf", width = 10, height = 8, onefile = F)
dca(LN_status ~ signature, 
    data = dat_TCGA,
    thresholds = seq(0, 0.3, by = 0.01), as_probability = "signature"
) %>%
  net_intervention_avoided(nper = 100) %>%
  #standardized_net_benefit() %>%
  plot(smooth = TRUE, type = 'net_intervention_avoided')

dev.off()





