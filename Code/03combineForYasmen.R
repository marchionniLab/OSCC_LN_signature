#################################################################
### Solit TCGA head and neck data for LN statu kTSP prediction and validation
### Luigi Marchionni
### Collaboration with Yasmen Ghantous (Sidransky lab)
### T1 and T2,  HPV negative,  oral cavity only


#################################################################
### Clean
rm(list=ls())

### Setwd
setwd("/Volumes/Macintosh/Research/Projects/HNSCC")

### Libraries
require(sampling)
require(limma)


#################################################################
### Load expression data and phenotypes
load("./Objs/combinedData.rda")


#################################################################
### Phenotype of interest

### Node status
pheno$NodeStatus <- pheno$NstagePath
levels(pheno$NodeStatus)[levels(pheno$NodeStatus) %in% paste("N", 1:4,  sep="")] <- "POS"
levels(pheno$NodeStatus)[levels(pheno$NodeStatus) != "POS"] <- "NEG"
table(pheno$NodeStatus)

### Modify site
levels(pheno$site)[levels(pheno$site) == "ORAL CAVITY"] <- "OC"

### SIMPLIFY Clinical N
levels(pheno$NstageClin)[levels(pheno$NstageClin) != "N0"] <- "POS"
levels(pheno$NstageClin)[levels(pheno$NstageClin) != "POS"] <- "NEG"

### Identifiy samples to be analyzed
attach(pheno)

### Explore variables
table(site)

### Keep only T1 and T2
smps <- (
	(Tstage %in% c("T1",  "T2")) &
	(HPV == "Neg" | is.na(HPV)) &
	(site %in% c("Lip", "OC")) &
	(! is.na(NodeStatus))
)

### Count
table(dataset, NodeStatus,  smps)


#################################################################
### Normalize and split
datN <- normalizeBetweenArrays(dat)[, smps]

### Pheno data
phenoN <- droplevels(pheno[ colnames(datN), ])


#################################################################
#################################################################
### Selection of variables 
sel <- c(
	"dataset",
	"Tstage",
	"HPV",
	"smoking",
	"smokingDetails",
	"NstageClin",
	"gender",
	"age",
	"OSevent",
	"OStime"
	 )

### Assemble in one data.frame and turn numeric all variables of interest
covs <- phenoN[,  sel]
summary(covs)

### Prepare vocs for sampling
covs <- sapply(covs , function(x) as.numeric(factor(paste(x))) )

### Outcome
outcome <- 1*(phenoN$NodeStatus != "NEG")


#################################################################
###SAMPLING

### Set seed
set.seed(736251)

### Balanced stratification
trainingOrTesting <- balancedstratification(
	covs[ , , drop=FALSE], strata=outcome,
	pik=inclusionprobabilities(1:nrow(covs), nrow(covs) * 0.5),
	comment=TRUE, method=1)

### Show
covsSplit <- apply(covs[ , drop=FALSE], 2, table, outcome, trainingOrTesting)


#################################################################
### Split data

### Subset Training
trainMatPooled <- datN[ , trainingOrTesting == 0]
trainPhenoPooled <- phenoN[ trainingOrTesting == 0,  ]
trainGrpPooled <- phenoN$NodeStatus[ trainingOrTesting == 0]
trainGrpPooled <- relevel(trainGrpPooled, "POS")

### Subset Testing
testMatPooled <- datN[ , trainingOrTesting == 1]
testPhenoPooled <- phenoN[ trainingOrTesting == 1,  ]
testGrpPooled <- phenoN$NodeStatus[ trainingOrTesting == 1]
testGrpPooled <- relevel(testGrpPooled, "POS")


#################################################################
### Split data by Study and platform
tcga <- phenoN$dataset == "TCGA"

### Subset Training
trainMatStudy <- datN[ , !tcga]
trainPhenoStudy <- phenoN[ !tcga,  ]
trainGrpStudy <- phenoN$NodeStatus[ !tcga]
trainGrpStudy <- relevel(trainGrpStudy, "POS")

### Save in text format
write.csv(trainMatStudy, file="./Text/arrayStudyExprs.csv")
write.csv(trainPhenoStudy, file="./Text/arrayStudyPheno.csv")

### Subset Testing
testMatStudy <- datN[ , tcga]
testPhenoStudy <- phenoN[ tcga,  ]
testGrpStudy <- phenoN$NodeStatus[ tcga]
testGrpStudy <- relevel(testGrpStudy, "POS")

### Save in text format
write.csv(testMatStudy, file="./Text/rnaSeqStudyExprs.csv")
write.csv(testPhenoStudy, file="./Text/rnaSeqStudyPheno.csv")


#################################################################
### Save

### Pooled
save(list=ls(pattern="Pooled"),  file="./Objs/lymphnodeDataPooled.rda")

### Study
save(list=ls(pattern="Study"),  file="./Objs/lymphnodeDataStudy.rda")


#################################################################
### Session Info
date()
sessionInfo()

### Date
q("no")
