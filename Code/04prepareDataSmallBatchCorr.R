################################################################
### Solit TCGA head and neck data for LN statu kTSP prediction and validation
### Luigi Marchionni
### Collaboration with Yasmen Ghantous (Sidransky lab)


#################################################################
### Clean
rm(list=ls())

### Setwd
setwd("/Volumes/Macintosh/Research/Projects/HNSCC")

### Libraries
require(sva)
require(limma)


#################################################################
### Load expression data and phenotypes
load("./Objs/tcgaData.rda")

### Edit phenotypes data.frame
phenoTCGA$PatientIDIdentifier <- gsub("-..$",  "",  phenoTCGA$sampleID)
any(duplicated(phenoTCGA$PatientIDIdentifier))

### Add tissue type
phenoTCGA$TissueType <- factor(gsub(".+-",  "",  phenoTCGA$sampleID))
levels(phenoTCGA$TissueType)[levels(phenoTCGA$TissueType) == "01"] <- "TUM"
levels(phenoTCGA$TissueType)[levels(phenoTCGA$TissueType) == "06"] <- "METS"
levels(phenoTCGA$TissueType)[levels(phenoTCGA$TissueType) == "11"] <- "NORM"
table(phenoTCGA$TissueType)

### Add Grade
phenoTCGA$Grade <- factor(phenoTCGA$neoplasm_histologic_grade)

### Tissue composition: lymphocytes
sel <- grep("percent_lymphocyte",  colnames(phenoTCGA))
phenoTCGA$LymphocytesPerc <- apply(phenoTCGA[, sel], 1, median, na.rm=TRUE)

### Tissue composition: monocyte
sel <- grep("percent_monocyte",  colnames(phenoTCGA))
phenoTCGA$MonocytesPerc <- apply(phenoTCGA[, sel], 1, median, na.rm=TRUE)

### Tissue composition: neutrophil
sel <- grep("percent_neutrophil",  colnames(phenoTCGA))
phenoTCGA$NeutrophilsPerc <- apply(phenoTCGA[, sel], 1, median, na.rm=TRUE)

### Tissue composition: normal cells
sel <- grep("percent_normal_cells",  colnames(phenoTCGA))
phenoTCGA$NormalCellsPerc <- apply(phenoTCGA[, sel], 1, median, na.rm=TRUE)

### Tissue composition: stromal cells
sel <- grep("percent_stromal_cells",  colnames(phenoTCGA))
phenoTCGA$StromalCellsPerc <- apply(phenoTCGA[, sel], 1, median, na.rm=TRUE)

### Tissue composition: tumor cells
sel <- grep("percent_tumor_cells",  colnames(phenoTCGA))
phenoTCGA$TumorCellsPerc <- apply(phenoTCGA[, sel], 1, median, na.rm=TRUE)

### Tissue composition: necrosis
sel <- grep("percent_necrosis",  colnames(phenoTCGA))
phenoTCGA$NecrosisPerc <- apply(phenoTCGA[, sel], 1, median, na.rm=TRUE)

### Tissue composition: tumor nuclei
sel <- grep("percent_tumor_nuclei",  colnames(phenoTCGA))
phenoTCGA$TumorNucleiPerc <- apply(phenoTCGA[, sel], 1, median, na.rm=TRUE)


#################################################################
### Load phenotypes from Farhoud
load("objs/lymphnodeDataStudy.rda")

### Rename   (testPhenoStudy is the tcga pheno from Farhoud)
targets <- testPhenoStudy

### Clean
rm(list=c("sel", ls(pattern="train"), ls(pattern="test")))


#################################################################
### Compare targets and UCSC Xena phenotypes

### Rename column and compare
phenoTCGA$patientID <- phenoTCGA$PatientIDIdentifier
length(intersect(phenoTCGA$patientID, targets$patientID))
all(targets$patientID %in% phenoTCGA$patientID)

### Selection fo columns
sel <- c("patientID",  "sampleID",  "TissueType", "Grade",
	 grep(".+Perc$",  colnames(phenoTCGA), value=TRUE))
sel

### Merge
pheno <- merge(targets, phenoTCGA[,  sel],  all.x=TRUE,  all.y=FALSE)

### Remove useless columns
pheno <- pheno[,  sapply(pheno, function(x) length(unique(x)) > 1)]

### Add rownames
rownames(pheno) <- gsub("-",  ".",  pheno$sampleID)

### Check
dim(pheno)
sum(duplicated(pheno$PatientIDIdentifier))


#################################################################
### Drop zeros
zeros <- rowSums(mrnaTCGA)
mrnaTCGA <- mrnaTCGA[ zeros > 0, ]

#### quantile normalize
mrnaTCGA <- normalizeBetweenArrays(1+mrnaTCGA)


#################################################################
### Subset gene expression
testMat <- mrnaTCGA[,  rownames(pheno)]
testPheno <- droplevels(pheno)
testGrp <- testPheno$NodeStatus
testGrp <- relevel(testGrp, "POS")

### Check
all(colnames(testMat) == rownames(testPheno))


########################################################