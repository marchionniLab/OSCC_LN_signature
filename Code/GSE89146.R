
rm(list = ls())
library(GEOquery)



GSE89146 <- getGEO("GSE89146", GSEMatrix = T, AnnotGPL = T)
GSE89146 <- GSE89146$GSE89146_series_matrix.txt.gz


Expr <- exprs(GSE89146)
Pheno <- pData(GSE89146)
FeatData <- fData(GSE89146)


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

#########################
## Process the phenotype
#table(Pheno$`tumor site:ch1`)
#table(Pheno$`hpv status:ch1`)
table(Pheno$`Stage:ch1`)

#Pheno <- Pheno[Pheno$`tumor site:ch1` == "Oral Cavity" & Pheno$`hpv status:ch1` == "HPV-", ]
#Pheno <- Pheno[Pheno$`tumor status:ch1` %in% c("T1", "T2"), ]


Pheno <- Pheno[!(Pheno$`Stage:ch1` == "NA"), ]

Pheno$NodeStatus <- Pheno$`Stage:ch1`
Pheno$NodeStatus <- gsub(".+N", "N", Pheno$NodeStatus)
Pheno$NodeStatus <- gsub("M.+", "", Pheno$NodeStatus)

table(Pheno$NodeStatus)

Pheno$NodeStatus <- as.factor(Pheno$NodeStatus)
levels(Pheno$NodeStatus) <- c("NEG", "POS", "POS", "POS")


Expr <- Expr[, colnames(Expr) %in% rownames(Pheno)]
all(rownames(Pheno) == colnames(Expr))

#############################
## Test
load("./Objs/FinalClassifiers.rda")


ClassifierGenes <- as.vector(ArrayKTSP$TSPs)
#InterSect <- intersect(ClassifierGenes, rownames(Expr))

#keep <- names(ArrayKTSP$score)[c(1,3:6)]

# ArrayKTSPFilt <- ArrayKTSP
# ArrayKTSPFilt$score <- ArrayKTSPFilt$score[keep]
# ArrayKTSPFilt$TSPs <- ArrayKTSPFilt$TSPs[keep, ]
# ArrayKTSPFilt$tieVote <- droplevels(ArrayKTSPFilt$tieVote[keep])
# ArrayKTSPFilt$name <- paste(nrow(ArrayKTSPFilt$TSPs), "TSPS")  


### Compute the sum and find the best threshold: ALL TRAINING SAMPLES
ktspStatsGSE89146 <- SWAP.KTSP.Statistics(
  inputMat = Expr,
  classifier = ArrayKTSP,
  CombineFunc = sum)

summary(ktspStatsGSE89146$statistics)


### Print ROC curve local maximas
auc(roc(Pheno$NodeStatus, ktspStatsGSE89146$statistics, levels=c("POS", "NEG"),
        direction="<"))

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
GSE89146Prediction <- SWAP.KTSP.Classify(
  Expr,
  ArrayKTSP,
  DecisionFunc = function(x) sum(x) > 2.5 )

Pheno$NodeStatus <- factor(Pheno$NodeStatus, levels = c("POS", "NEG"))
table(Pheno$NodeStatus)

### Resubstitution performance in the TRAINING set
confusionMatrix(GSE89146Prediction, Pheno$NodeStatus, positive="POS")
