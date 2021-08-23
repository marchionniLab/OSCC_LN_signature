
rm(list = ls())
setwd("/Volumes/Macintosh/Research/Projects/HNSCC")

library(ggplot2)
library(switchBox)
library(superheat)
library(reshape2)
library(pROC)
library(plotROC)
library(dplyr)
library(ggforce)

## Load the 2 classifiers and data
load("./Objs/FinalClassifiers.rda")


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


# Samples <- rownames(Array_Stats)
# Array_Stats <- apply(Array_Stats, 2, as.integer)
# rownames(Array_Stats) <- Samples

# Plot
png(filename = "./Figs/Heatmaps/ArrayKTSP_Heatmap.png", width = 3000, height = 2000, res = 200)
superheat(Array_Stats, col.dendrogram = F, heat.pal = c("royalblue4", "darkolivegreen3"), heat.pal.values = c(0, 1), yr = rowSums(Array_Stats), yr.plot.type  = "bar", left.label.text.col = c("red", "blue"), membership.rows = ArrayGroup, bottom.label.text.size = 3.5, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the array dataset")
dev.off()

pdf("./Figs/Heatmaps/SuppFigure2.pdf", width = 10, height = 9, onefile = F)
superheat(Array_Stats, col.dendrogram = F, heat.pal = c("royalblue4", "darkolivegreen3"), heat.pal.values = c(0, 1), yr = rowSums(Array_Stats), yr.plot.type  = "bar", left.label.text.col = c("red", "blue"), membership.rows = ArrayGroup, bottom.label.text.size = 2, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap of the TSPs votes in the array dataset")
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
superheat(TCGA_Stats, col.dendrogram = F, heat.pal.values = c(0, 1), yr = rowSums(TCGA_Stats), yr.plot.type  = "bar", left.label.text.col = c("red", "blue"), membership.rows = tcgaGroup, bottom.label.text.size = 3.5, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the TCGA dataset")
dev.off()

pdf("./Figs/Heatmaps/SuppFigure3.pdf", width = 10, height = 9, onefile = F)
superheat(TCGA_Stats, col.dendrogram = F, heat.pal = c("royalblue4", "darkolivegreen3"), heat.pal.values = c(0, 1), yr = rowSums(TCGA_Stats), yr.plot.type  = "bar", left.label.text.col = c("red", "blue"), membership.rows = tcgaGroup, bottom.label.text.size = 2, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the TCGA dataset")
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


pdf("./Figs/Heatmaps/Figure1.pdf", width = 10, height = 9, onefile = F)
superheat(PCR_Stats, col.dendrogram = F,  heat.pal = c("royalblue4", "darkolivegreen3"), heat.pal.values = c(0, 1), yr = rowSums(PCR_Stats), yr.plot.type  = "bar", left.label.text.col = c("red", "blue"), membership.rows = pcrGroup, bottom.label.text.size = 2, yr.num.ticks = 4, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the PCR (Testing) dataset")
dev.off()

#################################################################################################
################################################################################################

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
forLegend_Train <- apply(rbind(
  ci(roc(ArrayGroup, ktspStatsArray$statistics, levels = c("NEG", "POS"), direction = ">")),
  ci(roc(tcgaGroup, ktspStatsTCGA$statistics, levels = c("NEG", "POS"), direction = ">"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})

forLegend_Test <- apply(rbind(
  ci(roc(pcrGroup, ktspStatsPCR$statistics, levels = c("NEG", "POS"), direction = ">"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})
#################################################################
### ROC curves Using ggplot2

### Training1
datArray_KTSP <- melt(data.frame(
  ## Training Group
  ArrayGroup= factor(ArrayGroup, levels = c("NEG", "POS")),
  ## Mechanistic KTSP SUM training
  ktspStatsArray= ktspStatsArray$statistics))
### Change Colnames
colnames(datArray_KTSP) <- c("Status", "Dataset", "KTSP_sum")


### Training2
datTcga_KTSP <- melt(data.frame(
  ## Testing group
  Testing= factor(tcgaGroup, levels = c("NEG", "POS")),
  ## Mechanistic KTSP SUM training
  ktspStatsTCGA=ktspStatsTCGA$statistics))
### Change Colnames
colnames(datTcga_KTSP) <- c("Status", "Dataset", "KTSP_sum")

### Combine
dat_KTSP_Train <- rbind(datArray_KTSP, datTcga_KTSP)
dat_KTSP_Train$Status <- as.numeric(dat_KTSP_Train$Status)-1


### Testing
datPCR_KTSP <- melt(data.frame(
  ## Testing group
  Testing= factor(pcrGroup, levels = c("NEG", "POS")),
  ## Mechanistic KTSP SUM training
  ktspStatsPCR=ktspStatsPCR$statistics))
### Change Colnames
colnames(datPCR_KTSP) <- c("Status", "Dataset", "KTSP_sum")

### Combine
datPCR_KTSP$Status <- as.numeric(datPCR_KTSP$Status)-1

### Replace levels
levels(dat_KTSP_Train$Dataset) <- c("Training (Array)", "Training (TCGA)")
levels(dat_KTSP_Train$Dataset) <- paste(levels(dat_KTSP_Train$Dataset), forLegend_Train)

levels(datPCR_KTSP$Dataset) <- c("Testing (RT-PCR)")
levels(datPCR_KTSP$Dataset) <- paste(levels(datPCR_KTSP$Dataset), forLegend_Test)

#################################################################
### Plot Curve
png("./Figs/CompareAUCggplot_Training.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "K-TSP Performance in the training datasets"
#legendTitle <- "Dataset"
### Plot
basicplot_KTSP <- ggplot(dat_KTSP_Train, aes(d=Status, m=KTSP_sum, color=Dataset,
                                       linetype = Dataset)) +
  geom_roc(increasing = F) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  scale_color_manual(values = c("darkred", "darkblue")) +
  #scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_KTSP
### Close device
dev.off()

png("./Figs/CompareAUCggplot_Testing.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "K-TSP Performance in the RT-PCR dataset"
#legendTitle <- "Dataset"
### Plot
basicplot_PCR <- ggplot(datPCR_KTSP, aes(d=Status, m=KTSP_sum, color=Dataset,
                                             linetype = Dataset)) +
  geom_roc(increasing = F) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  scale_color_manual(values = c("darkgreen")) +
  #scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_PCR
### Close device
dev.off()


pdf("./Figs/Heatmaps/Figure1.pdf", width = 10, height = 9, onefile = F)
basicplot_KTSP
dev.off()


pdf("./Figs/Heatmaps/Figure3.pdf", width = 10, height = 9, onefile = F)
basicplot_PCR
dev.off()

save(basicplot_KTSP, file = "./Objs/BasicPlot_KTSP.rda")

####################################################################################################
##################
### Boxplots

## Array

## Which TSPs
i <- 1:nrow(ArrayKTSP$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ArrayKTSP$TSPs)



# Assemble
dfTspArray <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out$diff = out[, 1] - out[, 2]
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=ArrayMat, g=ArrayGroup)

## Reduce
datArray <- Reduce("rbind", dfTspArray)

datArray_diff <- datArray[datArray$variable == 'diff', ]

## Rename columns
colnames(datArray)[colnames(datArray) %in% c("variable", "value")] <- c("Gene", "Expression")
colnames(datArray_diff)[colnames(datArray_diff) %in% c("variable", "value")] <- c("Gene", "Expression")

## Make paired boxplot
png("./Figs/BoxPlots/ArrayBoxPlot_diff.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datArray_diff), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

########################
## TCGA

## Assemble
dfTspTCGA <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out$diff = out[, 1] - out[, 2]
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=tcgaMat, g=tcgaGroup)

## Reduce
datTCGA <- Reduce("rbind", dfTspTCGA)
datTCGA_diff <- datTCGA[datTCGA$variable == 'diff', ]


## Rename columns
colnames(datTCGA)[colnames(datTCGA) %in% c("variable", "value")] <- c("Gene", "Expression")
colnames(datTCGA_diff)[colnames(datTCGA_diff) %in% c("variable", "value")] <- c("Gene", "Expression")

## Make paired boxplot
png("./Figs/BoxPlots/TCGABoxPlot_diff.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datTCGA_diff), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

########################
## PCR

## Load the RT-PCR data
pcrData <- read.delim("./Data/RT_PCR/deltaCt.txt")

# Get the node status
pcrGroup <- pcrData$NodeStatus
table(pcrGroup)
# Rename the groups from 0,1 to neg,pos
pcrGroup[pcrGroup == 0] <- "NEG"
pcrGroup[pcrGroup == 1] <- "POS"

table(pcrGroup)
pcrGroup = factor(pcrGroup, levels = c("POS", "NEG"))

# Get the matrix
pcrMat <- pcrData[, -14]
rownames(pcrMat) <- pcrMat$SampleName
pcrMat$SampleName <- NULL

# Transpose the matrix
pcrMat <- t(pcrMat)

# One gene is misspelled, rename it
rownames(pcrMat)[rownames(pcrMat) == "TFGB2"] <- "TGFB2"

pcrMat <- log2(pcrMat+1)
#pcrMat <- t(scale(t(pcrMat), scale = T, center = T))

## Assemble
dfTspPCR <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out$diff = out[, 1] - out[, 2]
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=pcrMat, g=pcrGroup)

## Reduce
datPCR <- Reduce("rbind", dfTspPCR)
datPCR_diff <- datPCR[datPCR$variable == 'diff', ]

## Rename columns
colnames(datPCR_diff)[colnames(datPCR_diff) %in% c("variable", "value")] <- c("Gene", "Expression")
colnames(datPCR_diff)[colnames(datPCR_diff) %in% c("variable", "value")] <- c("Gene", "Expression")


## Make paired boxplot
png("./Figs/BoxPlots/PCRBoxPlot_diff.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datPCR_diff), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  #boxplot_framework(upper_limit = 1,
  #                  lower_limit = -1, fill_var = 'Group', logY = T) +
  #geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
  facet_wrap(~pair, scales = "free_y", nrow = 2) +
  #scale_y_log10() +
  coord_cartesian(ylim = c(-0.5, 0.5))
  #scale_y_continuous(limits = function(x){c(0.1, max(0.1, x))}) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()




