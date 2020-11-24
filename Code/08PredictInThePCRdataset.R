
###########################################################################################

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
confusionMatrix(pcrPrediction, pcrGroup, positive="NEG")
