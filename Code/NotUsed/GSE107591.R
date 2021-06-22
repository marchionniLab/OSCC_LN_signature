rm(list = ls())
library(GEOquery)



GSE107591 <- getGEO("GSE107591", GSEMatrix = T, AnnotGPL = T)
GSE107591 <- GSE107591$GSE107591_series_matrix.txt.gz


Expr <- exprs(GSE107591)
Pheno <- pData(GSE107591)
FeatData <- fData(GSE107591)


###########################
## Process the expression
rownames(Expr) <- FeatData$`Gene symbol`
summary(is.na(rownames(Expr)))
rownames(Expr) <- gsub("-","", rownames(Expr))
rownames(Expr) <- gsub("_","",rownames(Expr))
sel <- which(apply(Expr, 1, function(x) all(is.finite(x)) ))
Expr <- Expr[sel, ]
Expr <- Expr[!is.na(rownames(Expr)),]
Expr <- Expr[!(rownames(Expr) == ""), ] 
dim(Expr)

range(Expr)

###############
## Collapse duplicate genes by variance
summary(duplicated(rownames(Expr)))

Expr_var <- apply(Expr, 1, var)
Expr <- Expr[order(Expr_var, decreasing = T), ]
Expr <- Expr[!duplicated(rownames(Expr)), ]



### Normalize to GAPDH
Expr <- sweep(Expr, 2, Expr["GAPDH",],  "-")
Expr <- Expr[ -grep("GAPDH",  rownames(Expr)), ]

#####################################
### Filter to classifier genes
load("./Objs/FinalClassifiers.rda")

Genes <- as.vector(ArrayKTSP$TSPs)

GenesInter <- rownames(Expr) %in% Genes

Expr_Filt <- Expr[GenesInter, ]
Expr_Filt
Expr_Filt <- t(scale(t(Expr_Filt), center = T, scale = T))

#########################
## Process the phenotype
table(Pheno$`site:ch1`)
#table(Pheno$`hpv status:ch1`)
table(Pheno$`t:ch1`)

Pheno <- Pheno[Pheno$`site:ch1` == "Oral Cavity", ]
Pheno <- Pheno[Pheno$`t:ch1` %in% c("T1", "T2"), ]

table(Pheno$`n:ch1`)

#Pheno <- Pheno[!(Pheno$`Stage:ch1` == "NA"), ]

Pheno$NodeStatus <- Pheno$`n:ch1`


Pheno$NodeStatus <- as.factor(Pheno$NodeStatus)
levels(Pheno$NodeStatus) <- c("POS", "NEG")
Pheno$NodeStatus <- factor(Pheno$NodeStatus, levels = c("NEG", "POS"))

Expr_Filt <- Expr_Filt[, colnames(Expr_Filt) %in% rownames(Pheno)]
all(rownames(Pheno) == colnames(Expr_Filt))

#############################
## Test
load("./Objs/FinalClassifiers.rda")


ClassifierGenes <- as.vector(ArrayKTSP$TSPs)
InterSect <- intersect(ClassifierGenes, rownames(Expr_Filt))

keep <- names(ArrayKTSP$score)[c(1:5)]

ArrayKTSPFilt <- ArrayKTSP
ArrayKTSPFilt$score <- ArrayKTSPFilt$score[keep]
ArrayKTSPFilt$TSPs <- ArrayKTSPFilt$TSPs[keep, ]
ArrayKTSPFilt$tieVote <- droplevels(ArrayKTSPFilt$tieVote[keep])
ArrayKTSPFilt$name <- paste(nrow(ArrayKTSPFilt$TSPs), "TSPS")


### Compute the sum and find the best threshold: ALL TRAINING SAMPLES
ktspStatsGSE107591 <- SWAP.KTSP.Statistics(
  inputMat = Expr_Filt,
  classifier = ArrayKTSPFilt,
  CombineFunc = sum)

summary(ktspStatsGSE107591$statistics)


### Print ROC curve local maximas
auc(roc(Pheno$NodeStatus, ktspStatsGSE107591$statistics, levels=c("POS", "NEG"),
        direction="<"))

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
GSE107591Prediction <- SWAP.KTSP.Classify(
  Expr_Filt,
  ArrayKTSPFilt,
  DecisionFunc = function(x) sum(x) > 2.5 )

Pheno$NodeStatus <- factor(Pheno$NodeStatus, levels = c("POS", "NEG"))
table(Pheno$NodeStatus)

### Resubstitution performance in the TRAINING set
confusionMatrix(GSE107591Prediction, Pheno$NodeStatus, positive="POS")
