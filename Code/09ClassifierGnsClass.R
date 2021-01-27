

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

