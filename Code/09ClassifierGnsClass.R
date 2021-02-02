
rm(list = ls())
setwd("/Volumes/Macintosh/Research/Projects/HNSCC")


load("./Objs/allTSPs_Annot.rda")
load("./Objs/myTSPs.rda")
load("./Objs/FinalClassifiers.rda")


ClassifierGoodGns <- ArrayKTSP$TSPs[,1]
ClassifierBadGns <- ArrayKTSP$TSPs[,2]


AllGoodGns <- allTSPs[, c(1,3)]
AllGoodGns <- as.data.frame(AllGoodGns)


GoodGnsType <- AllGoodGns[AllGoodGns$GoodGene %in% ClassifierGoodGns, ]
table(GoodGnsType$GoodGene, GoodGnsType$Type)



AllBadGns <- allTSPs[, c(2,3)]
AllBadGns <- as.data.frame(AllBadGns)


BadGnsType <- AllBadGns[AllBadGns$BadGene %in% ClassifierBadGns, ]
table(BadGnsType$BadGene, BadGnsType$Type)

summary(ClassifierBadGns %in% myTSPs[,2])
summary(ClassifierBadGns %in% myTSPs[,1])

summary(ClassifierGoodGns %in% myTSPs[,1])

# PDLIM4 should be good  >> found in classifer bad &  SKI > found in both
# GSTP1 should be good >> found in classifer bad  &   CTSB > found in both

################
# Which bad gene lists contain SKI and CTSP
SKI_CTSP_Class <- AllBadGns[AllBadGns$BadGene %in% ClassifierGoodGns, ]
table(SKI_CTSP_Class$BadGene, SKI_CTSP_Class$Type)
# SKI < TAGs (oncogene)
# CTSP < emt (pro-invasion)

################
# Which good gene lists contain PDLIM4 and GSTP1
PDLIM4_GSTP1_Class <- AllGoodGns[AllGoodGns$GoodGene %in% ClassifierBadGns, ]
table(PDLIM4_GSTP1_Class$GoodGene, PDLIM4_GSTP1_Class$Type)
# PDLIM4 < immune-surveillance (immune-surveillance on) and TAGs (tumor suppressor)
# GSTP1 < emt (anti-invasion) and growth (anti-growth)
# TGFB2 < apoptosis (pro-apoptotic)



###############
## Which mechanism ??
SKI_PDLIM4 <- allTSPs[which(allTSPs[,2] == "SKI" & allTSPs[,1] == "PDLIM4"), ]
SKI_PDLIM4  #TAGs : SKI bad ! (oncogene) & PDLIM4 good ! (tumor suppressor)

STK11_INO80 <- allTSPs[which(allTSPs[,1] == "STK11" & allTSPs[,2] == "INO80"), ]
STK11_INO80 # Growth : STK11 good & INO80 Bad

SETMAR_RAB11FIP4 <- allTSPs[which(allTSPs[,1] == "SETMAR" & allTSPs[,2] == "RAB11FIP4"), ]
SETMAR_RAB11FIP4 # Growth: SETMAR: Good & RAB11FIP4 bad

SATB1_TMEM176B <- allTSPs[which(allTSPs[,1] == "SATB1" & allTSPs[,2] == "TMEM176B"), ]
SATB1_TMEM176B # Immune : SATB1: Good & TMEM: Bad

PAPPA_TGFB2 <- allTSPs[which(allTSPs[,1] == "PAPPA" & allTSPs[,2] == "TGFB2"), ]
PAPPA_TGFB2 # EMT: PAPPA: Good & TGFB2: Bad


CTSB_GSTP1 <- allTSPs[which(allTSPs[,2] == "CTSB" & allTSPs[,1] == "GSTP1"), ]
CTSB_GSTP1 # EMT: CTSB : Bad !! & GSTP1 : Good !!


