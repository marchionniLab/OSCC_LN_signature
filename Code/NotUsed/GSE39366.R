

rm(list = ls())

library(GEOquery)



GSE39366 <- getGEO("GSE39366", GSEMatrix = T, AnnotGPL = T)
GSE39366 <- GSE39366$GSE39366_series_matrix.txt.gz


Expr <- exprs(GSE39366)
Pheno <- pData(GSE39366)
FeatData <- fData(GSE39366)


###########################
## Process the expression
rownames(Expr) <- FeatData$ORF
summary(is.na(rownames(Expr)))
rownames(Expr) <- gsub("-","", rownames(Expr))
rownames(Expr) <- gsub("_","",rownames(Expr))
sel <- which(apply(Expr, 1, function(x) all(is.finite(x)) ))
Expr <- Expr[sel, ]
Expr <- Expr[!is.na(rownames(Expr)),]
Expr <- Expr[!(rownames(Expr) == ""), ] 
dim(Expr)

range(Expr)

summary(duplicated(rownames(Expr)))

# ### By variance
Expr_var <- apply(Expr, 1, var)
Expr <- Expr[order(Expr_var, decreasing = T), ]
Expr <- Expr[!duplicated(rownames(Expr)), ]

### Normalize to GAPDH
Expr <- sweep(Expr, 2, Expr["GAPDH",],  "-")
Expr <- Expr[ -grep("GAPDH",  rownames(Expr)), ]

###########################################################
### Filter to classifier genes
load("./Objs/FinalClassifiers.rda")

Genes <- as.vector(ArrayKTSP$TSPs)

GenesInter <- rownames(Expr) %in% Genes

Expr_Filt <- Expr[GenesInter, ]
Expr_Filt
Expr_Filt <- t(scale(t(Expr_Filt), center = T, scale = T))
#########################
## Process the phenotype
table(Pheno$`tumor site:ch1`)
table(Pheno$`hpv status:ch1`)
table(Pheno$`tumor status:ch1`)

Pheno <- Pheno[Pheno$`tumor site:ch1` == "Oral Cavity" & Pheno$`hpv status:ch1` == "HPV-", ]
Pheno <- Pheno[Pheno$`tumor status:ch1` %in% c("T1", "T2"), ]

table(Pheno$`node status:ch1`)
#Pheno <- Pheno[!(Pheno$`node status:ch1` == "NA"), ]

Pheno$NodeStatus <- as.factor(Pheno$`node status:ch1`)
levels(Pheno$NodeStatus) <- c("NEG", "POS", "POS", "POS")


Expr_Filt <- Expr_Filt[, colnames(Expr_Filt) %in% rownames(Pheno)]
all(rownames(Pheno) == colnames(Expr_Filt))

#############################
## Test
load("./Objs/FinalClassifiers.rda")


ClassifierGenes <- as.vector(ArrayKTSP$TSPs)
InterSect <- intersect(ClassifierGenes, rownames(Expr_Filt))

keep <- names(ArrayKTSP$score)[c(1,3:6)]

ArrayKTSPFilt <- ArrayKTSP
ArrayKTSPFilt$score <- ArrayKTSPFilt$score[keep]
ArrayKTSPFilt$TSPs <- ArrayKTSPFilt$TSPs[keep, ]
ArrayKTSPFilt$tieVote <- droplevels(ArrayKTSPFilt$tieVote[keep])
ArrayKTSPFilt$name <- paste(nrow(ArrayKTSPFilt$TSPs), "TSPS")  


### Compute the sum and find the best threshold: ALL TRAINING SAMPLES
ktspStatsGSE39366 <- SWAP.KTSP.Statistics(
  inputMat = Expr_Filt,
  classifier = ArrayKTSPFilt,
  CombineFunc = sum)

summary(ktspStatsGSE39366$statistics)


### Print ROC curve local maximas
auc(roc(Pheno$NodeStatus, ktspStatsGSE39366$statistics, levels=c("POS", "NEG"),
        direction="<"))

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
GSE39366Prediction <- SWAP.KTSP.Classify(
  Expr,
  ArrayKTSPFilt,
  DecisionFunc = function(x) sum(x) > 2.5 )

Pheno$NodeStatus <- factor(Pheno$NodeStatus, levels = c("POS", "NEG"))
table(Pheno$NodeStatus)

### Resubstitution performance in the TRAINING set
confusionMatrix(GSE39366Prediction, Pheno$NodeStatus, positive="POS")
