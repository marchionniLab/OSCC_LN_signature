###################################################
### Luigi Marchionni
### Processed GEO data from series gse
### Clinical information from the matrix available from GEO

###################################################
### Clear Workspace
rm(list=ls())

### Date
date()

### Set WD
#setwd("~/Research/HeadNeck/pubDomain/GEOdata/GPL10588/GSE65858/")

### Load libraries
library("GEOquery")

 
##################################################
###################################################
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
###writinf the data to a file
write.table(pheno, file="./Objs/GSE65858.txt", sep="\t", row.names=F, col.names=T)


##################################################
### Save
save(esetgse, file="./Objs/esetGSE65858.rda")


###################################################
### Session information and clean quit

### Session
sessionInfo()

### Clean
rm(list=ls())

### Quit
q("no")
