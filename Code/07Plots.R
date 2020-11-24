
rm(list = ls())
setwd("/Volumes/Macintosh/Research/Projects/HNSCC")

library(ggplot2)
library(switchBox)
library(superheat)

## Load the 2 classifiers and data
load("./Objs/FinalClassifiers.rda")

############################################################################
## Scatter Plots
###########################################################################
### Plot genes in the training set (the 2 array datsets)
## Which TSPs
i <- 1:nrow(ArrayKTSP$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ArrayKTSP$TSPs)

## Assemble
dfTspArray <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), out)
}, x=ArrayMat, g=ArrayGroup)

names(dfTspArray) <- rownames(ArrayKTSP$TSPs)

# Change the names of elements inside each element in dfTspTrain (For Plotting)  
for(i in seq_along(dfTspArray)) names(dfTspArray[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")

## Reduce
datArray <- Reduce("rbind", dfTspArray)

##
## Scatter plot
png(filename = "./Figs/ScatterPlots/KTSP_Array_ScatterPlot.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datArray), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 12),
        axis.title=element_text(face="bold", size = 12),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=17.5),
        legend.title = element_text(face="bold", size=17.5),
        strip.text.x = element_text(face="bold", size=11))
sctplt
dev.off()

###########################################################################
### Plot genes in the test set (the tcga datset)
## Which TSPs
i <- 1:nrow(tcgaKTSP$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=tcgaKTSP$TSPs)

## Assemble
dfTspTcga <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair = paste("TSP:", paste(i, collapse = "-")), out)
}, x = tcgaMat, g = tcgaGroup)

names(dfTspTcga) <- rownames(tcgaKTSP$TSPs)

# Change the names of elements inside each element in dfTspTrain (For Plotting)  
for(i in seq_along(dfTspTcga)) names(dfTspTcga[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")

## Reduce
datTcga <- Reduce("rbind", dfTspTcga)

####
## Scatter plots
png(filename = "./Figs/ScatterPlots/KTSP_TCGA_ScatterPlot.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datTcga), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 12),
        axis.title=element_text(face="bold", size = 12),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=17.5),
        legend.title = element_text(face="bold", size=17.5),
        strip.text.x = element_text(face="bold", size=11))
sctplt
dev.off()

#################################
## Load the RT-PCR data
pcrData <- read.delim("./Data/RT_PCR/deltaCt.txt")

# Get the node status
pcrGroup <- pcrData$NodeStatus
table(pcrGroup)
# Rename the groups from 0,1 to neg,pos
pcrGroup[pcrGroup == 0] <- "NEG"
pcrGroup[pcrGroup == 1] <- "POS"

# Get the matrix
pcrMat <- pcrData[, -14]
rownames(pcrMat) <- pcrMat$SampleName
pcrMat$SampleName <- NULL

# Transpose the matrix
pcrMat <- t(pcrMat)

# One gene is misspelled, rename it
rownames(pcrMat)[rownames(pcrMat) == "TFGB2"] <- "TGFB2"

## Which TSPs
i <- 1:nrow(tcgaKTSP$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=tcgaKTSP$TSPs)

## Assemble
dfTspPCR <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), out)
}, x=pcrMat, g=pcrGroup)

names(dfTspPCR) <- rownames(tcgaKTSP$TSPs)

# Change the names of elements inside each element in dfTspTrain (For Plotting)  
for(i in seq_along(dfTspPCR)) names(dfTspPCR[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")

## Reduce
datPCR <- Reduce("rbind", dfTspPCR)

## Scatter plots
png(filename = "./Figs/ScatterPlots/KTSP_PCR_ScatterPlot.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datPCR), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 12),
        axis.title=element_text(face="bold", size = 12),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=17.5),
        legend.title = element_text(face="bold", size=17.5),
        strip.text.x = element_text(face="bold", size=11))
sctplt
dev.off()

##################################################################################################
### Heatmaps
##################################################################################################
#### Make a heatmap of the TSPs expression in the 3 datasets

## 1- Arrays
# Get the prediction stats
ktspStatsArray <- SWAP.KTSP.Statistics(
  inputMat = ArrayMat,
  classifier = ArrayKTSP,
  CombineFunc = sum)

Array_Stats <- ktspStatsArray$comparisons
Array_Stats <- Array_Stats*1

# Order the stats by the sum of votes
NewOrder <- order(rowSums(Array_Stats))
Array_Stats <- Array_Stats[NewOrder, ]

#GroupCol <- ArrayGroup[NewOrder]
#levels(GroupCol) <- c("red", "blue")

# Order the true class labels
ArrayGroup <- ArrayGroup[NewOrder]

# Get the predicted class labels
PredClass <- ifelse(rowSums(Array_Stats) >= 3, "NEG", "POS")
PredClass <- factor(PredClass, levels = c("POS", "NEG")) 
levels(PredClass) <- c("red", "blue")
# Check order 
all(rownames(Array_Stats) == names(PredClass)) # TRUE:: No Need to order them (Already ordered by sum of votes) 

# Plot
png(filename = "./Figs/Heatmaps/ArrayKTSP_Heatmap.png", width = 3000, height = 2000, res = 200)
superheat(Array_Stats, col.dendrogram = F, yr = rowSums(Array_Stats), yr.plot.type  = "bar", left.label.text.col = c("red", "blue"), membership.rows = ArrayGroup, bottom.label.text.size = 3.5, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the array dataset")
dev.off()

########
## 2- TCGA
# Get the prediction stats
ktspStatsTCGA <- SWAP.KTSP.Statistics(
  inputMat = tcgaMat,
  classifier = tcgaKTSP,
  CombineFunc = sum)

TCGA_Stats <- ktspStatsTCGA$comparisons
TCGA_Stats <- TCGA_Stats*1

# Order the stats by the sum of votes
NewOrder <- order(rowSums(TCGA_Stats))
TCGA_Stats <- TCGA_Stats[NewOrder, ]

# Order the true class labels
tcgaGroup <- tcgaGroup[NewOrder]

# Get the predicted class labels
PredClass <- ifelse(rowSums(TCGA_Stats) >= 3, "NEG", "POS")
PredClass <- factor(PredClass, levels = c("POS", "NEG")) 
#X <- X[NewOrder]
levels(PredClass) <- c("red", "blue")
# Check order 
all(rownames(TCGA_Stats) == names(PredClass)) # TRUE:: No Need to order them (Already ordered by sum of votes) 


png(filename = "./Figs/Heatmaps/TcgaKTSP_Heatmap.png", width = 3000, height = 2000, res = 200)
superheat(TCGA_Stats, col.dendrogram = F, yr = rowSums(TCGA_Stats), yr.plot.type  = "bar", left.label.text.col = c("red", "blue"), membership.rows = tcgaGroup, bottom.label.text.size = 3.5, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the TCGA dataset")
dev.off()

########
## 3- PCR
# Load the RT-PCR data
pcrData <- read.delim("./Data/RT_PCR/deltaCt.txt")

# Get the node status
pcrGroup <- pcrData$NodeStatus
table(pcrGroup)
# Rename the groups from 0,1 to neg,pos
pcrGroup[pcrGroup == 0] <- "NEG"
pcrGroup[pcrGroup == 1] <- "POS"
pcrGroup <- factor(pcrGroup, levels = c("POS", "NEG"))

# Get the matrix
pcrMat <- pcrData[, -14]
rownames(pcrMat) <- pcrMat$SampleName
pcrMat$SampleName <- NULL

# Transpose the matrix
pcrMat <- t(pcrMat)

# One gene is misspelled, rename it
rownames(pcrMat)[rownames(pcrMat) == "TFGB2"] <- "TGFB2"

# Get the prediction stats
ktspStatsPCR <- SWAP.KTSP.Statistics(
  inputMat = pcrMat,
  classifier = tcgaKTSP,
  CombineFunc = sum)

PCR_Stats <- ktspStatsPCR$comparisons
PCR_Stats <- PCR_Stats*1

# Order the stats by the sum of votes
NewOrder <- order(rowSums(PCR_Stats))
PCR_Stats <- PCR_Stats[NewOrder, ]

# Order the true class labels
pcrGroup <- pcrGroup[NewOrder]

# Get the predicted class labels
PredClass <- ifelse(rowSums(PCR_Stats) >= 3, "NEG", "POS")
PredClass <- factor(PredClass, levels = c("POS", "NEG")) 
#X <- X[NewOrder]
levels(PredClass) <- c("red", "blue")
# Check order 
all(rownames(PCR_Stats) == names(PredClass)) # TRUE:: No Need to order them (Already ordered by sum of votes) 


png(filename = "./Figs/Heatmaps/PcrKTSP_Heatmap.png", width = 3000, height = 2000, res = 200)
superheat(PCR_Stats, col.dendrogram = F, yr = rowSums(PCR_Stats), yr.plot.type  = "bar", left.label.text.col = c("red", "blue"), membership.rows = pcrGroup, bottom.label.text.size = 3.5, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the PCR (Testing) dataset")
dev.off()

#################################################################################################
## MDS Plots

# 1- Arrays
Array_Stats <- ktspStatsArray$comparisons
Array_Stats <- Array_Stats*1
Array_Stats <- as.data.frame(Array_Stats)
Array_Stats$Group <- ArrayGroup

Array_Stats_WithoutLabels <- Array_Stats
Array_Stats_WithoutLabels$Group <- NULL

Array_Stats_Dist <- dist(as.matrix(Array_Stats_WithoutLabels))

MDS_Array_Stats_WithoutLabels <- cmdscale(Array_Stats_Dist)

x <- MDS_Array_Stats_WithoutLabels[,1]
y <- MDS_Array_Stats_WithoutLabels[,2]

png(filename = "./Figs/MDS_Plots/MDS_Array.png", width = 2500, height = 2000, res = 300)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Array MDS", col = as.numeric(ArrayGroup), pch = 16)
#abline(0,1)
text(x, y, labels = ktspStatsArray$statistics, cex=1 ,pos = 1, col = as.numeric(ArrayGroup),)
dev.off()
#######
# 2- TCGA
Tcga_Stats <- ktspStatsTCGA$comparisons
Tcga_Stats <- Tcga_Stats*1
Tcga_Stats <- as.data.frame(Tcga_Stats)
Tcga_Stats$Group <- tcgaGroup

Tcga_Stats_WithoutLabels <- Tcga_Stats
Tcga_Stats_WithoutLabels$Group <- NULL

Tcga_Stats_Dist <- dist(as.matrix(Tcga_Stats_WithoutLabels))

MDS_Tcga_Stats_WithoutLabels <- cmdscale(Tcga_Stats_Dist)

x <- MDS_Tcga_Stats_WithoutLabels[,1]
y <- MDS_Tcga_Stats_WithoutLabels[,2]

png(filename = "./Figs/MDS_Plots/MDS_Tcga.png", width = 2500, height = 2000, res = 300)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Tcga MDS", col = as.numeric(tcgaGroup), pch = 16)
#abline(0,1)
text(x, y, labels = ktspStatsTCGA$statistics, cex=1 ,pos = 1, col = as.numeric(tcgaGroup),)
dev.off()

#######
# 2- PCR
PCR_Stats <- ktspStatsPCR$comparisons
PCR_Stats <- PCR_Stats*1
PCR_Stats <- as.data.frame(PCR_Stats)
PCR_Stats$Group <- pcrGroup

PCR_Stats_WithoutLabels <- PCR_Stats
PCR_Stats_WithoutLabels$Group <- NULL

PCR_Stats_Dist <- dist(as.matrix(PCR_Stats_WithoutLabels))

MDS_PCR_Stats_WithoutLabels <- cmdscale(PCR_Stats_Dist)

x <- MDS_PCR_Stats_WithoutLabels[,1]
y <- MDS_PCR_Stats_WithoutLabels[,2]

png(filename = "./Figs/MDS_Plots/MDS_PCR.png", width = 2500, height = 2000, res = 300)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="PCR MDS", col = as.numeric(pcrGroup), pch = 16)
#abline(0,1)
text(x, y, labels = ktspStatsPCR$statistics, cex=1 ,pos = 1, col = as.numeric(pcrGroup),)
dev.off()

#################################################################################################

##################################################################################################
### Boxplots 
##################################################################################################
## Make paired boxplot for the array data

## Which TSPs
i <- 1:nrow(ArrayKTSP$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ArrayKTSP$TSPs)

## Assemble
dfTspArray <- lapply(tsp, function(i, x, g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=ArrayMat, g=ArrayGroup)

## Reduce
datArray <- Reduce("rbind", dfTspArray)

## Rename columns
colnames(datArray)[colnames(datArray) %in% c("variable", "value")] <- c("Gene", "Expression")

#####
## Make paired boxplot
png("./Figs/BoxPlots/KTSP_Array_Boxplot.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datArray), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

##################################################################################################
## Make paired boxplot for the array data

## Which TSPs
i <- 1:nrow(tcgaKTSP$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=tcgaKTSP$TSPs)

## Assemble
dfTspTcga <- lapply(tsp, function(i, x, g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=tcgaMat, g=tcgaGroup)

## Reduce
datTcga <- Reduce("rbind", dfTspTcga)

## Rename columns
colnames(datTcga)[colnames(datTcga) %in% c("variable", "value")] <- c("Gene", "Expression")

#####
## Make paired boxplot
png("./Figs/BoxPlots/KTSP_TCGA_Boxplot.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datTcga), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

#################################################################################
## Make paired boxplot for the PCR data

## Which TSPs
i <- 1:nrow(ArrayKTSP$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ArrayKTSP$TSPs)

## Assemble
dfTspPCR <- lapply(tsp, function(i, x, g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=pcrMat, g=pcrGroup)

## Reduce
datPCR <- Reduce("rbind", dfTspPCR)

## Rename columns
colnames(datPCR)[colnames(datPCR) %in% c("variable", "value")] <- c("Gene", "Expression")

#####
## Make paired boxplot
png("./Figs/BoxPlots/KTSP_PCR_Boxplot.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datPCR), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()


################################################################################################
#########################################################################
######################################################

## GGPLOT to compare the two classifiers

# Get the stats
ktspStatsArray <- SWAP.KTSP.Statistics(
  inputMat = ArrayMat,
  classifier = ArrayKTSP,
  CombineFunc = sum)

ktspStatsTCGA <- SWAP.KTSP.Statistics(
  inputMat = tcgaMat,
  classifier = tcgaKTSP,
  CombineFunc = sum)

ktspStatsPCR <- SWAP.KTSP.Statistics(
  inputMat = pcrMat,
  classifier = ArrayKTSP,
  CombineFunc = sum)

### Prepare the legend
forLegend_KTSP <- apply(rbind(
  ci(roc(ArrayGroup, ktspStatsArray$statistics)),
  ci(roc(tcgaGroup, ktspStatsTCGA$statistics)),
  ci(roc(pcrGroup, ktspStatsPCR$statistics))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


#################################################################
### ROC curves Using ggplot2

### Training
datArray_KTSP <- melt(data.frame(
  ## Training Group
  ArrayGroup= ArrayGroup,
  ## Mechanistic KTSP SUM training
  ktspStatsArray= ktspStatsArray$statistics))
### Change Colnames
colnames(datArray_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")


### Testing
datTcga_KTSP <- melt(data.frame(
  ## Testing group
  Testing= tcgaGroup,
  ## Mechanistic KTSP SUM training
  ktspStatsTCGA=ktspStatsTCGA$statistics))
### Change Colnames
colnames(datTcga_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")

### Testing
datPCR_KTSP <- melt(data.frame(
  ## Testing group
  Testing= pcrGroup,
  ## Mechanistic KTSP SUM training
  ktspStatsPCR=ktspStatsPCR$statistics))
### Change Colnames
colnames(datPCR_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")

### Combine
dat_KTSP <- rbind(datArray_KTSP, datTcga_KTSP, datPCR_KTSP)
dat_KTSP$Status <- as.numeric(dat_KTSP$Status)-1

### Replace levels
levels(dat_KTSP$KTSP_type) <- gsub("Stats", "", levels(dat_KTSP$KTSP_type))
levels(dat_KTSP$KTSP_type) <- paste(levels(dat_KTSP$KTSP_type), forLegend_KTSP)

#################################################################
### Plot Curve
png("./Figs/CompareAUCggplot.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "K-TSP Performance in the training and testing datasets"
#legendTitle <- paste("Training (", "Array)",
#                     " - Training (", "TCGA)",  
#                     " - Testing (", "PCR)", sep="")
### Plot
basicplot_KTSP <- ggplot(dat_KTSP, aes(d=Status, m=KTSP_sum, color=KTSP_type,
                                       linetype = KTSP_type)) +
  geom_roc(cutoffs.at = seq(1,20,1)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  scale_color_manual(values = c("red", "red", "darkgreen")) +
  #scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(legendTitle, values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_KTSP
### Close device
dev.off()

save(basicplot_KTSP, file = "./Objs/KTSP/BasicPlot_KTSP.rda")


####################################################################################################
####################################################################################################
