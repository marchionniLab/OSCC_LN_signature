### Process the TCGA data obtained from UCSC cancer genome browser


######################################################################
### Clean working space ans set working directory
rm(list=ls())

setwd("/Volumes/Macintosh/Research/Projects/HNSCC")



######################################################################
### Process phenotype data
datDir <- "./Data/TCGA"
gzFile <- list.files(".", pattern="TCGA.+gz$", recursive=TRUE)

### Decompress
untar(gzFile, exdir=datDir)


###########################################################################
### Process the phenotype file

### Find phenotype file
phenoFile <- list.files(datDir, pattern="clinical_data$", full.names = TRUE, recursive=TRUE)

### Read file
pheno <- read.table(phenoFile, sep="\t", header=TRUE, colClasses = "character", fill=TRUE)

### Set NA values
pheno[pheno == ""] <- NA

### Set numeric values
pheno <- data.frame(sapply(pheno, function(x) {
    if (all(is.na(as.numeric(x)))) x
    else as.numeric(x)
}, simplify=FALSE), stringsAsFactors=FALSE)

### Get rid of useless columns (all redundant information)
dim(pheno)
pheno <- pheno[ , apply(pheno, 2, function(x) length(unique(x)) > 1 ) ]
dim(pheno)

### Set rownames
rownames(pheno) <- gsub("-", ".", pheno$sampleID)

### Rename
phenoTCGA <- pheno


######################################################################
### Process mRNA

### List files with data
dataFile <- list.files(datDir, pattern=".+Matrix$", full.names = TRUE, recursive=TRUE)

### Read file
dat <- read.table(dataFile, sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, quote="")

### Set rownames
rownames(dat) <- dat$sample

### Drop first colunm with gene names
dat$sample <- NULL

### transform into a matrix
dat <- as.matrix(dat)

### Rename
mrnaTCGA <- dat


###########################################################################
### Make sure data and phenotypes are in the same order
commonID <- intersect(rownames(phenoTCGA), colnames(mrnaTCGA))
length(commonID)

### Subset accordingly
phenoTCGA <- phenoTCGA[ commonID , ]
mrnaTCGA <- mrnaTCGA[ , commonID ]

### Verify
all( rownames(phenoTCGA) == colnames(mrnaTCGA) )

###########################################################################
### Remove decompressed archive
#unlink(datDir, recursive=TRUE)


######################################################################
### Save objects
save(list=ls(pattern="TCGA$"), file="./Objs/tcgaData.rda")

################################################################
################################################################
### Load expression data and phenotypes from Farhoud dataset

### Base path
bPath1 <- "./Data/fromFarhoud/"

### Load Farhoud data
load(paste(bPath1,  "2016.12.28_AllData_clin_matched_pre_correction_n1135.RData",  sep=""))

### Rename
dat <- combinednew
rm(combinednew)

### Load clinical info
pheno <- read.csv(paste(bPath1,  "2017.02.08_Combined_clin_dat_v13.csv",  sep=""),
                  stringsAsFactors=FALSE)

### Check conformity and reorder if needed
all(colnames(dat) == pheno$sample_id)
all(colnames(dat) %in% pheno$sample_id)
dat <- dat[ ,  pheno$sample_id]
all(colnames(dat) == pheno$sample_id)



#################################################################
### Add gender
pheno$gender <- factor(pheno$sex)
#levels(pheno$gender) <- c("F", "F",  "M", "M")

### Add site
pheno$site <- factor(pheno$site)
levels(pheno$site)

### Add HPV
pheno$HPV <- factor(pheno$hpv_status)
levels(pheno$HPV)

### Add T stage
pheno$Tstage <- factor(paste("T", pheno$t_stage, sep=""))

levels(pheno$Tstage) <- gsub("[ab]$",  "",  levels(pheno$Tstage))
levels(pheno$Tstage)[levels(pheno$Tstage) %in% c("TNA",  "TX")] <- NA
levels(pheno$Tstage)

### Add clincal n_stage
pheno$NstageClin <- factor(rep(NA, nrow(pheno)))
levels(pheno$NstageClin)

### Add path n_stage
pheno$NstagePath <- factor(toupper(paste("N", pheno$n_stage, sep="")))
levels(pheno$NstagePath)
levels(pheno$NstagePath)[levels(pheno$NstagePath) %in% c("NNA",  "NX")] <- NA
levels(pheno$NstagePath) <- gsub("[ABC]$",  "",  levels(pheno$NstagePath))
levels(pheno$NstagePath)

### Add smoking
pheno$smoking <- factor(pheno$tobacco)
levels(pheno$smoking)[levels(pheno$smoking) %in% c("1",  "CURRENT",  "FORMER")] <- "YES"
levels(pheno$smoking)[levels(pheno$smoking) %in% c("0",  "NEVERSMOKER")] <- "NO"
levels(pheno$smoking)[levels(pheno$smoking) %in% c("888")] <- NA

### Add smoking details
pheno$smokingDetails <- factor(pheno$tobacco)
levels(pheno$smokingDetails)

### Add pack years
pheno$packyears <- as.numeric(pheno$packyears)
summary(pheno$packyears)

### Add OS time
pheno$OStime <- pheno$os_mo
summary(pheno$OStime)

### Add OS event
pheno$OSevent <- factor(pheno$os_event)
levels(pheno$OSevent)

### Add dataset
pheno$dataset <- factor(pheno$source)
levels(pheno$dataset)

### Add rownames
rownames(pheno) <- pheno$sample_id

#################################################################
### Check conformity
all(colnames(dat) == rownames(pheno))

### turn to matrix
dat <- data.matrix(dat)


#################################################################
### Scale and log transforms if necessary
maxDat <- apply(dat, 2,  max)
minDat <- apply(dat, 2,  min)
tapply(maxDat, pheno$dataset, mean)
tapply(minDat, pheno$dataset, mean)

### ### Need to add 2 to TCGA: better shift everything (rather than bringing up TCGA)
### dat[ , pheno$dataset == "TCGA"] <- 1+dat[ , pheno$dataset == "TCGA"]

### Floor all datasets
dat <- 2+dat

### Need to rescale this: GSE40774
dat[ , pheno$dataset == "GSE40774"] <- log2(dat[ , pheno$dataset == "GSE40774"])


#################################################################
#### Remove duplicated samples in  TCGA
dropTCGA <- grep("^TCGA.+06$", colnames(dat), value=TRUE)
str(dropTCGA)
dat <- dat[ ,  !colnames(dat) %in% dropTCGA]
pheno <- pheno[ !rownames(pheno) %in% dropTCGA, ]


#### Remove duplicated samples in Chung's and TCGA
### Collapsing functions
source("./Code/collapseFunc.R")

### Identify data to be collapsed
dupPt <- pheno$pt_id[duplicated(pheno$pt_id)]
dupSmp <- pheno$sample_id[ pheno$pt_id %in% dupPt]
dupPt <- pheno[ pheno$pt_id %in% dupPt, "pt_id"]

### Tmp exprs
tmp <- dat[, colnames(dat) %in% dupSmp]
dim(tmp)

### Collapse expression by median
tmp <- t(apply(tmp, 1, collapseData,  ind=dupPt, func=median))
colnames(tmp) <- collapseData(dupSmp, dupPt, func=function(z) z[[1]] )

### Remove from dat
dim(dat)
dat <- dat[ , !colnames(dat) %in% dupSmp ]

### Add median sample and subset pheno
dat <- cbind(tmp, dat)
pheno <- pheno[ colnames(dat), ]

### Check conformity
any(duplicated(pheno$patientID))
all(rownames(pheno) == colnames(dat))

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
### Split data by Study and platform
tcga <- phenoN$dataset == "TCGA"

### Subset Testing
tcgaMatStudy <- datN[ , tcga]
tcgaPhenoStudy <- phenoN[ tcga,  ]
tcgaGrpStudy <- phenoN$NodeStatus[ tcga]
tcgaGrpStudy <- relevel(tcgaGrpStudy, "POS")


#################################################################
### Save
save(list=ls(pattern="Study"),  file="./Objs/tcga_Farhood.rda")


######################################################################
### Session information
sessionInfo()

### Clean working space and quit
rm(list=ls())
q("no")
