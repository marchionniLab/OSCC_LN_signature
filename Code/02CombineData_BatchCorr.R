################################################################
### Combine datasets together 
### Filter to samples of interest: HPV-ve T1-T2 OSCC with LN status 
### Correct for batch effect


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
### Load phenotypes information for the TCGA dataset from Farhoud
load("objs/tcga_Farhood.rda")

### Rename   (testPhenoStudy is the tcga pheno from Farhoud)
targets <- tcgaPhenoStudy
targets$patientID <- targets$pt_id

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
sum(duplicated(pheno$patientID))

# Remove duplicated samples
pheno <- pheno[!duplicated(pheno$patientID), ]

#################################################################
### Drop zeros
zeros <- rowSums(mrnaTCGA)
mrnaTCGA <- mrnaTCGA[ zeros > 0, ]

#### quantile normalize
mrnaTCGA <- normalizeBetweenArrays(1+mrnaTCGA)


#################################################################
### Subset gene expression
tcgaMat <- mrnaTCGA[,  rownames(pheno)]
tcgaPheno <- droplevels(pheno)
tcgaGrp <- tcgaPheno$NodeStatus
tcgaGrp <- relevel(tcgaGrp, "POS")

### Check
all(colnames(tcgaMat) == rownames(tcgaPheno))


########################################################
#################################################################
#################################################################
### Load expression data and phenotypes from the two microarray datasets
load("./Objs/ArrayData.rda")


#################################################################
### Phenotype of interest

### Node status
pheno$NodeStatus <- pheno$NstagePath
levels(pheno$NodeStatus)[levels(pheno$NodeStatus) %in% paste("N", 1:4,  sep="")] <- "POS"
levels(pheno$NodeStatus)[levels(pheno$NodeStatus) != "POS"] <- "NEG"
table(pheno$NodeStatus)

### Identifiy samples to be analyzed
attach(pheno)

### Explore variables
table(site)
table(HPV)
table(Tstage)

### Keep only T1 and T2,  ecc
smps <- (
  (Tstage %in% c("T1",  "T2")) &
    (HPV == "Negative" | is.na(HPV)) &
    (site %in% c("Lip", "ORAL CAVITY")) &
    (! is.na(NodeStatus))
)

### Count
table(dataset, NodeStatus,  smps)

### Subset
pheno <- pheno[ smps, ]


#################################################################
### Subset gene expression
arrayMat <- dat[,  rownames(pheno)]
arrayPheno <- droplevels(pheno)
arrayGrp <- arrayPheno$NodeStatus
arrayGrp <- relevel(arrayGrp, "POS")

### Check
all(colnames(tcgaMat) == rownames(tcgaPheno))
all(colnames(arrayMat) == rownames(arrayPheno))
all(colnames(arrayMat) == pheno$sampleID)


#################################################################
### Common Genes
gns <- intersect(rownames(arrayMat), rownames(tcgaMat))
str(gns)

### Subset
arrayMat <- arrayMat[gns, ]
tcgaMat <- tcgaMat[gns, ]


#################################################################
### Combine and combat and psva correct

### Matrix
mat <- cbind(arrayMat, tcgaMat)

### Factor
grp <- factor(c(arrayGrp, tcgaGrp))
levels(grp) <- c("POS", "NEG")

### Source
dataset <- c(pheno$dataset,  rep("TCGA", ncol(tcgaMat)))


#################################################################
### Scale and log transforms if necessary
maxDat <- apply(mat, 2,  max)
minDat <- apply(mat, 2,  min)
tapply(maxDat, dataset, mean)
tapply(minDat, dataset, mean)


#############################################################
### Remove batch: run pSVA
psva.mat <- psva(mat,  batch=dataset)
colnames(psva.mat) <- colnames(mat)

### Remove batch: run Combat
combat.mat <- ComBat(mat, batch=dataset)


#################################################################
### Save
save(mat, combat.mat, psva.mat, grp, dataset,  file="objs/Data_BatchCorr.rda")

## Save the pheno also for later use
save(arrayPheno, tcgaPheno, file = "./Objs/Pheno.rda")
#################################################################
### Session Info
sessionInfo()



