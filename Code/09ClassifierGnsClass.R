

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

