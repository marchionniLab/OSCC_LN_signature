#################################################################
### Solit TCGA head and neck data for LN statu kTSP prediction and validation
### Luigi Marchionni
### Collaboration with Yasmen Ghantous (Sidransky lab)


#################################################################
### Clean
rm(list=ls())

### Setwd
#setwd("/Volumes/Macintosh/Research/Projects/HNSCC")

### Libraries
require(sampling)
require(Biobase)
require(limma)
library(sva)


################################################################
################################################################
### Load expression data and phenotypes from Wichmann's dataset: GSE65858

### Load GSE65858
load("./Objs/esetGSE65858.rda")

### Expression
dat1 <- exprs(esetgse)
rownames(dat1) <- fData(esetgse)$Symbol

### Pheno data
pheno1 <- pData(esetgse)


################################################################
################################################################
### Load expression data and phenotypes from Hollinger's dataset: GSE42743

### Load GSE42743
load("./Objs/GSE42743FromLuigi/esetGSE42743.rda")

### Expression
dat2 <- exprs(esetgse)
rownames(dat2) <- fData(esetgse)$SYMBOL

### Pheno data
pheno2 <- pData(esetgse)


#################################################################
### Combine
pheno <- data.frame(sampleID=c(rownames(pheno1),  rownames(pheno2)),
		    stringsAsFactors=FALSE)

### Add patient ID
pheno$age <- c(pheno1$age,  pheno2$age.dx)

### Add gender
pheno$gender <- factor(c(as.character(pheno1$gender),  as.character(pheno2$gender)))
levels(pheno$gender) <- c("F", "F",  "M", "M")

### Add site
pheno$site <- factor(c(as.character(pheno1$tumor_site),  pheno2$oc.vs..op))
levels(pheno$site) <- toupper(levels(pheno$site))
levels(pheno$site)[levels(pheno$site) == "CAVUM ORIS"] <- "ORAL CAVITY"

### Add HPV
pheno$HPV <- factor(c(as.character(pheno1$hpv_dna),  rep(NA,  nrow(pheno2))))
levels(pheno$HPV)[levels(pheno$HPV) != "Negative"] <- "Positive"

### Add T stage
pheno$Tstage <- factor(c(paste("T", pheno1$t_category, sep=""),
		       as.character( pheno2$t.stage)))
levels(pheno$Tstage) <- gsub("[ab]$",  "",  levels(pheno$Tstage))
levels(pheno$Tstage)[levels(pheno$Tstage) %in% c("TNA",  "TX")] <- NA
levels(pheno$Tstage)

### Add path n_stage
pheno$NstagePath <- factor(toupper(c(paste("N", pheno1$n_category, sep=""), 
			   as.character( pheno2$pn.stage))))
levels(pheno$NstagePath)
levels(pheno$NstagePath)[levels(pheno$NstagePath) %in% c("NNA",  "NX")] <- NA
levels(pheno$NstagePath) <- gsub("[ABC]$",  "",  levels(pheno$NstagePath))
levels(pheno$NstagePath)

### Add smoking
pheno$smoking <- factor(c(as.character(pheno1$smoking), pheno2$smoking.status))
levels(pheno$smoking)[levels(pheno$smoking) %in% c("Yes",  "CURRENT",  "FORMER")] <- "YES"
levels(pheno$smoking)[levels(pheno$smoking) %in% c("No",  "NEVERSMOKER")] <- "NO"

### Add pack years
pheno$packyears <- c(as.numeric(pheno1$packyears), rep(NA,  nrow(pheno2)))
summary(pheno$packyears)

### Add OS time
pheno$OStime <- c(round(pheno1$os/365,1), round(pheno2$futime/365,1))
summary(pheno$OStime)

### Add OS event
pheno$OSevent <- factor(c(as.character(pheno1$os_event), pheno2$survivallastfollowup))
levels(pheno$OSevent) <- gsub("^DEA.+",  "Yes",  levels(pheno$OSevent))
levels(pheno$OSevent) <- gsub("^DIED.+",  "Yes",  levels(pheno$OSevent))
levels(pheno$OSevent) <- gsub("TRUE",  "Yes",  levels(pheno$OSevent))
levels(pheno$OSevent) <- gsub("^LIVING.+",  "No",  levels(pheno$OSevent))
levels(pheno$OSevent) <- gsub("FALSE",  "No",  levels(pheno$OSevent))
levels(pheno$OSevent)

### Add dataset
pheno$dataset <- factor(c(rep("GSE65858", nrow(pheno1)),  rep("GSE42743", nrow(pheno2))))
levels(pheno$dataset)

### Add rownames
rownames(pheno) <- pheno$sampleID


#################################################################
### Select common genes
gns <- intersect(rownames(dat1),  rownames(dat2))
str(gns)

### Combine expression
dat <- cbind(dat1[gns,] ,  dat2[gns, ])

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


### Check conformity
any(duplicated(pheno$patientID))
all(rownames(pheno) == colnames(dat))


#############################################################
### Remove batch: run pSVA
psva.dat <- psva(dat,  batch=pheno$dataset)
colnames(psva.dat) <- colnames(dat)

### Remove batch: run Combat
combat.dat <- ComBat(dat, batch=pheno$dataset)


#################################################################
### Save
save(dat, pheno, combat.dat,  psva.dat,  file="./Objs/combinedDataSmall.rda")


#################################################################
### Session Info
date()
sessionInfo()

### Quit
q("no")
