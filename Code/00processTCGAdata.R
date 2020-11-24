######################################################################
######################################################################
### October 07, 2015
### Process TCGA data obtained from UCSC cancer genome browser
### Luigi Marchionni


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



######################################################################
### Session information
sessionInfo()

### Clean working space and quit
rm(list=ls())
q("no")
