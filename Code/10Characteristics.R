
### Clean
rm(list = ls())

### Setwd
setwd("/Volumes/Macintosh/Research/Projects/HNSCC")


## Load the Phenotype data of both Array and TCGA data
load("objs/Pheno.rda")

ArrayPheno <- trainPheno
tcgaPheno <- testPheno

rm(trainPheno, testPheno)

### Combine together
tcgaPheno$patientID <- NULL
tcgaPheno$smokingDetails <- NULL
ArrayDataset <- ArrayPheno$dataset
ArrayPheno$dataset <- NULL
ArrayPheno$dataset <- ArrayDataset
tcgaPheno$dataset <- rep("TCGA", nrow(tcgaPheno))

all(colnames(tcgaPheno) == colnames(ArrayPheno))

AllPheno <- rbind(ArrayPheno, tcgaPheno)


########
## Get the characteristics
NofSamples <- nrow(AllPheno) # 189 samples

MeanAge <- mean(AllPheno$age[!is.na(AllPheno$age)]) # 60.5 years
sdAge <- sd(AllPheno$age[!is.na(AllPheno$age)])


Gender <- table(AllPheno$gender)
Gender # F:62, M: 127

PerMales <- (table(AllPheno$gender)["M"]/nrow(AllPheno)) *100
PerMales


HPV_status <- summary(is.na(AllPheno$HPV))
HPV_status # 48: Negative -- 141 NAs

T_Stage <- table(AllPheno$Tstage)
T_Stage # T1: 34 - T2: 155

N_Stage <- table(AllPheno$NstagePath)
N_Stage # N0:130 - N1: 27 - N2: 32

Smoking <- table(AllPheno$smoking)
Smoking # Yes: 130 - No: 54 - NAs : 5

PercSmokers <- (table(AllPheno$smoking)["YES"]/nrow(AllPheno)) *100
PercSmokers


meanPackYear <- mean(!is.na(AllPheno$packyears))
meanPackYear

######################
## Each  dataset
GSE65858Pheno <- ArrayPheno[ArrayPheno$dataset == "GSE65858", ]
GSE42743Pheno <- ArrayPheno[ArrayPheno$dataset == "GSE42743", ]

######
## Age
GSE65858_MeanAge <- mean(GSE65858Pheno$age[!is.na(GSE65858Pheno$age)]) # 56 years
GSE42743_MeanAge <- mean(GSE42743Pheno$age[!is.na(GSE42743Pheno$age)]) # 59 years
TCGA_MeanAge <- mean(tcgaPheno$age[!is.na(tcgaPheno$age)]) # 62.5 years

#######
## % males
PercMalesGSE65858 <- (table(GSE65858Pheno$gender)["M"]/nrow(GSE65858Pheno)) *100
PercMalesGSE65858

PercMalesGSE42743 <- (table(GSE42743Pheno$gender)["M"]/nrow(GSE42743Pheno)) *100
PercMalesGSE42743

PercMalesTCGA <- (table(tcgaPheno$gender)["M"]/nrow(tcgaPheno)) *100
PercMalesTCGA

########
## Smokers
PercSmokersGSE65858 <- (table(GSE65858Pheno$smoking)["YES"]/nrow(GSE65858Pheno)) *100
PercSmokersGSE65858

PercSmokersGSE42743 <- (table(GSE42743Pheno$smoking)["YES"]/nrow(GSE42743Pheno)) *100
PercSmokersGSE42743

PercSmokersTCGA <- (table(tcgaPheno$smoking)["YES"]/nrow(tcgaPheno)) *100
PercSmokersTCGA

#########
## T-stage
Tstage_GSE65858 <- table(GSE65858Pheno$Tstage)
Tstage_GSE65858

Tstage_GSE42743 <- table(GSE42743Pheno$Tstage)
Tstage_GSE42743

Tstage_TCGA <- table(tcgaPheno$Tstage)
Tstage_TCGA


#########
## N-stage
Nstage_GSE65858 <- table(GSE65858Pheno$NstagePath)
Nstage_GSE65858

Nstage_GSE42743 <- table(GSE42743Pheno$NstagePath)
Nstage_GSE42743

Nstage_TCGA <- table(tcgaPheno$NstagePath)
Nstage_TCGA

#####################################3
## Load the RT-PCR data
library("readxl")
pcrPheno <- read_xlsx("./Data/RT_PCR/Table 1.xlsx")

#####
## Age
PCR_MeanAge <- mean(pcrPheno$age) # 60.4 years

#######
## % males
PercMalesPCR <- (table(pcrPheno$Gender)["M"]/nrow(pcrPheno)) *100
PercMalesPCR # 74.3%

#########
## T-stage
Tstage_PCR <- table(pcrPheno$`T staging`)
Tstage_PCR

#########
## N-stage
Nstage_PCR <- table(pcrPheno$`N staging`)
Nstage_PCR

Nstage_GSE42743 <- table(GSE42743Pheno$NstagePath)
Nstage_GSE42743

Nstage_TCGA <- table(tcgaPheno$NstagePath)
Nstage_TCGA
