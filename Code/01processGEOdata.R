#################################################################
## GSE65858
## GSE42743


#################################################################
### Clean
rm(list=ls())

### Setwd
#setwd("/Volumes/Macintosh/Research/Projects/HNSCC")

### Libraries
require(GEOquery)
require(sampling)
require(Biobase)
require(limma)
library(sva)

###############################################################
## GSE65858
### Get the gse
gse <- getGEO("GSE65858")#, getGPL = FALSE,  AnnotGPL = FALSE)

### Check names
names(gse)

### Get the Expression Set
esetgse <- gse[[1]]

### Check length of platform lists: if more than 1 merge
length(esetgse)

##################################################
### Process esetgse: it has length == 1

### Check phenotypes
pheno <- pData(esetgse)
dim(pheno)
str(pheno)

### Remove useless columns and factors
pheno <- sapply(pheno, as.character)
pheno <- pheno[, -grep("characteristics_ch",  colnames(pheno))]
pheno <- pheno[,  apply(pheno, 2, function(x) length(unique(x)) > 1)]

### Change values 
pheno[pheno == "NA"] <- NA

### Make into dataset
pheno <- data.frame(pheno,  stringsAsFactors=TRUE)

### Edit names
colnames(pheno) <- gsub("\\.ch1$",  "", colnames(pheno))
rownames(pheno) <- pheno$geo_accession


##################################################
### Process the pheno data.frame to assemble phenotypes correctly

### Some variable should be character  (ie IDs)
sel <- c("title", "geo_accession", "source_name_ch1",  "ID")
pheno[,sel] <- apply(pheno[,sel], 2, as.character)

### Some variable should be numeric (ie time to recurrence)
sel <- c("age", "os", "pfs",  "packyears")
pheno[,sel] <- apply(pheno[,sel], 2, function(x) as.numeric(as.character(x)))

### Some variable should be factors (ie time to recurrence)
sel <- grep("_mutation", colnames(pheno))
pheno[,sel] <- apply(pheno[,sel], 2, as.factor)


##################################################
### Check order
if ( all(sampleNames(esetgse) == rownames(pheno)) ) {
  print("All rows in phenoData are in the good order")
} else if ( all(sampleNames(esetgse) %in% rownames(pheno)) ) {
  pheno <- pheno[sampleNames(esetgse), ]
  print("Reordering samples")
}  else stop("check the row order")


##################################################
### Create a new ExpressionSet instance
pData(esetgse) <- pheno

### Check validity
validObject(esetgse)


##################################################
### writing the data to a file
write.table(pheno, file="./Objs/GSE65858.txt", sep="\t", row.names=F, col.names=T)


##################################################
### Save
save(esetgse, file="./Objs/esetGSE65858.rda")

################################################################
################################################################
## GSE42743
###function to retrive the information from the messed up pheno
parsePheno <- function(myPat, pheno) {
  out <- rep(NA, length=nrow(pheno))
  sel <- unlist(apply(pheno, 2, grep, pattern=myPat, value=FALSE))
  val <- unlist(apply(pheno, 2, grep, pattern=myPat, value=TRUE))
  out[sel] <- gsub(".+: ", "", val)
  out <- toupper(out)
}

### Get the gse
gse <- getGEO("GSE42743")#,  getGPL = F,  GSEMatrix=F)

### Check names
names(gse)

### Get the Expression Set
esetgse <- gse[[1]]

### Check length of platform lists: if more than 1 merge
length(esetgse)

#######################
### Process esetgse: it has length == 1

### Check phenotypes
pheno <- pData(esetgse)
dim(pheno)
str(pheno)

### Remove useless columns
pheno <- pheno[, apply(pheno, 2, function(x) length(unique(x)) > 1 )
               | colnames(pheno) == "source_name_ch1" ]

######################
### Process the pheno data.frame to assemble phenotypes correctly

### What kind of information is available: columns with potential
### Phenotypic information
phenoTags <- apply(pheno[, grep("characteristics", colnames(pheno)), drop=FALSE], 2, gsub,
                   pattern=":.+", replacement="")
phenoTags <- unique(as.character(phenoTags))
phenoTags <- phenoTags[!is.na(phenoTags) & phenoTags != ""]
phenoTags


### Check if the information is present also in other columns
### Not labeled by the description tag in GEO
missedTags <- apply(pheno[, -grep("characteristics", colnames(pheno)), drop=FALSE], 2, gsub,
                    pattern=":.+", replacement="")
missedTags <- unique(as.character(missedTags))
missedTags <- missedTags[!is.na(missedTags) & missedTags != ""]
str(missedTags)

### If the following is TRUE you can go ahead with phenoTags only
if (! any(missedTags %in% phenoTags)) {
  print("go ahead 'phenoTags' contains all the needed information")
} else {
  stop("Check 'missedTags' for phenotypic information")
}

### Show tags found in the pheno data.frame
sort(phenoTags)


################
### IMPORTANT: checking for metacharacters in phenotypes tags
### It is crucial to remove special characters BEFORE the semicolon!!!
########

### The list of special characters
specialChar <- c("\\.", "\\\\", "\\|", "\\(", "\\)", "\\[", "\\{",
                 "\\^", "\\$", "\\*", "\\+", "\\?")

### Check for metacharactes in the phenotypes tags
anyMetaChar <- sapply(specialChar, function(x, y) {
  length(grep(x, y)) > 0
}, y = phenoTags)
anyMetaChar <- names(anyMetaChar[anyMetaChar])

### Show matacharacters used in 'phenoTags' if any
anyMetaChar

### Replace metacharacter in 'pheno' and 'phenoTags' if these are present
if (length(anyMetaChar) > 0) {
  ## Process phenotypes tags
  phenoTags <- sapply(phenoTags, function(x, y) {
    y <- paste("[", paste(y, collapse=""), "]", sep="")
    gsub(y, replacement="", x)
  }, y=anyMetaChar)
  ## Process phenotypes data.frame
  rnms <- rownames(pheno)
  pheno <- sapply(pheno, function(x, y) {
    y <- paste("[", paste(y, collapse=""), "]", sep="")
    gsub(y, replacement="", x)
  }, y=anyMetaChar)
  ## Add ronames back to pheno
  rownames(pheno) <- rnms
}


####################
### Resume processing phenotype information
####################

### Parse pheno data.frame
parsedPheno <- sapply(phenoTags, parsePheno, pheno)
colnames(parsedPheno) <- gsub(" ", ".", colnames(parsedPheno))

### Remove factors
parsedPheno <- data.frame(apply(parsedPheno, 2, as.character), stringsAsFactors=FALSE)
parsedPheno$order1 <- 1:nrow(parsedPheno)
rownames(parsedPheno) <- rownames(pheno)

### If there are Values like N/A or empty strings "", or UNKNOWN: use NA instead
parsedPheno[parsedPheno=="NA"] <- NA
parsedPheno[parsedPheno=="UNKNOWN"] <- NA
parsedPheno[parsedPheno=="N/A"] <- NA
parsedPheno[parsedPheno==""] <- NA

### Some variable should be numeric (ie time to recurrence)
sel <- c("age.dx", "futime", "order1")
parsedPheno[,sel] <- apply(parsedPheno[,sel], 2,
                           function(x) as.numeric((gsub("\\s", "", x))))

######################
### Modidify other columns: stage
parsedPheno$t.stage <- factor(paste("T", parsedPheno$t.stage,  sep=""))

### Modidify other columns: clinical node status
parsedPheno$cn.stage <- factor(paste("N", parsedPheno$cn.stage,  sep=""))

### Modidify other columns: pathological node status
parsedPheno$pn.stage <- factor(paste("N", parsedPheno$pn.stage,  sep=""))

### Modidify other columns: gender
parsedPheno$gender <- factor(parsedPheno$gender)


#####################
### Prepare extra info about samples and other GEO metadata
extraInfo <- as.data.frame(pheno[, -grep("characteristics", colnames(pheno))],
                           stringsAsFactors=FALSE)
extraInfo <- data.frame(apply(extraInfo, 2, as.character),stringsAsFactors=FALSE)
extraInfo$order2 <- 1:nrow(extraInfo)
rownames(extraInfo) <- extraInfo[, "geo_accession"]

### Merge the metadata with rest of the phenotypic information
finalPheno <- merge(extraInfo,  parsedPheno, by=0, all=FALSE,
                    stringsAsFactors=FALSE)

### Reorder and then check order
finalPheno <- finalPheno[order(as.numeric(finalPheno$order1)),]
all(finalPheno$order1 == finalPheno$order2)
all(finalPheno[, "Row.names"] %in% rownames(pheno))
all(finalPheno[, "Row.names"] == rownames(pheno))

### Set rownames
if (all(finalPheno[, "Row.names"] == rownames(pheno))) {
  rownames(finalPheno) <- rownames(pheno)
  print("All rows in phenoData are in the good order")
} else {
  stop("check the row order")
}


#########################
### Create a new ExpressionSet instance
pData(esetgse) <- finalPheno


#########################
### Add information about genes
tmp <- esetgse

### Load annotated eset
load("./data/GSE3292/esetGSE3292.rda")

### Check order of features
feat <- intersect(featureNames(tmp), featureNames(esetgse))

### Extract
ann <- featureData(esetgse)[feat, ]

###Assign
featureData(tmp) <- ann

### Check validity
validObject(tmp)

### replace and save
esetgse <- tmp

### Check validity
validObject(esetgse)


###############################
###writinf the data to a file
write.table(finalPheno, file="./Objs/GSE42743.txt", sep="\t", row.names=F, col.names=T)

##################################################
### Save
save(esetgse, file="./Objs/esetGSE42743.rda")


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
save(dat, pheno, combat.dat,  psva.dat,  file="./Objs/ArrayData.rda")


#################################################################
### Session Info
date()
sessionInfo()

### Quit
q("no")
