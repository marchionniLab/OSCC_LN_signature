rm(list = ls())
library(GEOquery)



GSE30788 <- getGEO("GSE30788", GSEMatrix = T, AnnotGPL = T)
GSE30788_1 <- GSE30788$`GSE30788-GPL13952_series_matrix.txt.gz`


Expr <- exprs(GSE30788_1)
Pheno <- pData(GSE30788_1)
FeatData <- fData(GSE30788_1)


###########################
## Process the expression
rownames(Expr) <- FeatData$gene_symbol
summary(is.na(rownames(Expr)))
rownames(Expr) <- gsub("-","", rownames(Expr))
rownames(Expr) <- gsub("_","",rownames(Expr))
sel <- which(apply(Expr, 1, function(x) all(is.finite(x)) ))
Expr <- Expr[sel, ]
Expr <- Expr[!is.na(rownames(Expr)),]
Expr <- Expr[!(rownames(Expr) == ""), ] 
dim(Expr)

range(Expr)

#Expr <- log2(Expr)
#########################
## Process the phenotype
table(Pheno$`location:ch1`)
#table(Pheno$`hpv:ch1`)
table(Pheno$`pt:ch1`)

Pheno <- Pheno[Pheno$`location:ch1` == "oral cavity", ]
Pheno <- Pheno[Pheno$`pt:ch1` %in% c("1", "2"), ]

table(Pheno$`pn (pathological n-status = lymph node metastasis status):ch1`)

#Pheno <- Pheno[!(Pheno$`Stage:ch1` == "NA"), ]

Pheno$NodeStatus <- Pheno$`pn (pathological n-status = lymph node metastasis status):ch1`


Pheno$NodeStatus <- as.factor(Pheno$NodeStatus)
levels(Pheno$NodeStatus) <- c("NEG", "POS")
#Pheno$NodeStatus <- factor(Pheno$NodeStatus, levels = c("NEG", "POS"))

Expr <- Expr[, colnames(Expr) %in% rownames(Pheno)]
all(rownames(Pheno) == colnames(Expr))

#############################
## Test
load("./Objs/FinalClassifiers.rda")


ClassifierGenes <- as.vector(ArrayKTSP$TSPs)
InterSect <- intersect(ClassifierGenes, rownames(Expr))

keep <- names(ArrayKTSP$score)[c(1,3,5,6)]

ArrayKTSPFilt <- ArrayKTSP
ArrayKTSPFilt$score <- ArrayKTSPFilt$score[keep]
ArrayKTSPFilt$TSPs <- ArrayKTSPFilt$TSPs[keep, ]
ArrayKTSPFilt$tieVote <- droplevels(ArrayKTSPFilt$tieVote[keep])
ArrayKTSPFilt$name <- paste(nrow(ArrayKTSPFilt$TSPs), "TSPS")


### Compute the sum and find the best threshold: ALL TRAINING SAMPLES
ktspStatsGSE30788 <- SWAP.KTSP.Statistics(
  inputMat = Expr,
  classifier = ArrayKTSPFilt,
  CombineFunc = sum)

summary(ktspStatsGSE30788$statistics)


### Print ROC curve local maximas
auc(roc(Pheno$NodeStatus, ktspStatsGSE30788$statistics, levels=c("POS", "NEG"),
        direction="<"))

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
GSE30788Prediction <- SWAP.KTSP.Classify(
  Expr,
  ArrayKTSPFilt,
  DecisionFunc = function(x) sum(x) > 2.5 )

Pheno$NodeStatus <- factor(Pheno$NodeStatus, levels = c("POS", "NEG"))
table(Pheno$NodeStatus)

### Resubstitution performance in the TRAINING set
confusionMatrix(GSE30788Prediction, Pheno$NodeStatus, positive="POS")

###############################

