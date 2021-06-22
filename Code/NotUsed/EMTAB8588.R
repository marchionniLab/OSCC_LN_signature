rm(list = ls())
library(GEOquery)



Expr <- read.delim("./Data/EMTAB8588/normalized_matrix.txt")


Pheno <- read.delim("./Data/EMTAB8588/Pheno.txt")

ArrayDesign <- read.delim("./Data/EMTAB8588/ArrayDesign.txt")

###########################
## Process the expression
Common <- intersect(ArrayDesign$Array.Design.Name, rownames(Expr))
ArrayDesign <- ArrayDesign[ArrayDesign$Array.Design.Name %in% Common, ]

Expr <- as.matrix(Expr)

#ArrayDesign <- ArrayDesign[!duplicated(ArrayDesign$X.2), ]
#Common2 <- intersect(ArrayDesign$Array.Design.Name, rownames(Expr))
Expr <- Expr[Common, ]
all(rownames(Expr) == ArrayDesign$Array.Design.Name)


rownames(Expr) <- ArrayDesign$X.2
summary(is.na(rownames(Expr)))
#rownames(Expr) <- gsub("-","", rownames(Expr))
#rownames(Expr) <- gsub("_","",rownames(Expr))
sel <- which(apply(Expr, 1, function(x) all(is.finite(x)) ))
Expr <- Expr[sel, ]
Expr <- Expr[!is.na(rownames(Expr)),]
Expr <- Expr[!(rownames(Expr) == ""), ] 
dim(Expr)

range(Expr)

################
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
table(Pheno$Characteristics.sampling.site.)
Pheno <- Pheno[Pheno$Characteristics.sampling.site. == "tumour", ]

table(Pheno$Characteristics.disease.)
table(Pheno$Characteristics.pt.)

Pheno <- Pheno[Pheno$Characteristics.disease. %in% c("carcinoma af the buccal mucosa", "carcinoma of floor of mouth", "carcinoma of the base of tongue", "lip carcinoma", "tongue carcinoma", "tonsillar carcinoma"), ]
Pheno <- Pheno[Pheno$Characteristics.pt. %in% c("T1", "T2"), ]

table(Pheno$Characteristics.pn.)

#Pheno <- Pheno[!(Pheno$`Stage:ch1` == "NA"), ]

Pheno$NodeStatus <- Pheno$Characteristics.pn.


Pheno$NodeStatus <- as.factor(Pheno$NodeStatus)
levels(Pheno$NodeStatus) <- c("NEG", "POS", "POS", "POS", "POS")
#Pheno$NodeStatus <- factor(Pheno$NodeStatus, levels = c("NEG", "POS"))
rownames(Pheno) <- Pheno$Source.Name

colnames(Expr_Filt) <- gsub("X", "", colnames(Expr_Filt))

CommSamples <- intersect(rownames(Pheno), colnames(Expr_Filt))

Expr_Filt <- Expr_Filt[, CommSamples]
Pheno <- Pheno[CommSamples, ]

all(rownames(Pheno) == colnames(Expr_Filt))

#Expr <- as.matrix(Expr)
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
ktspStats <- SWAP.KTSP.Statistics(
  inputMat = Expr_Filt,
  classifier = ArrayKTSPFilt,
  CombineFunc = sum)

summary(ktspStats$statistics)


### Print ROC curve local maximas
auc(roc(Pheno$NodeStatus, ktspStats$statistics, levels=c("POS", "NEG"),
        direction="<"))

### Get prediction based on best threshold from ROC curve
### Note the use of ">"
Prediction <- SWAP.KTSP.Classify(
  Expr_Filt,
  ArrayKTSPFilt,
  DecisionFunc = function(x) sum(x) > 2.5 )

Pheno$NodeStatus <- factor(Pheno$NodeStatus, levels = c("POS", "NEG"))
table(Pheno$NodeStatus)

### Resubstitution performance in the TRAINING set
confusionMatrix(Prediction, Pheno$NodeStatus, positive="POS")

###############################

