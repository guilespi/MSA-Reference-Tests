library(lattice)
library(reshape)
library(Hmisc)
library(car)

scores <- read.csv("/Users/guilespi/Documents/Development/interrupted/bioinformatics/bio1/multiple_seq_aligner/results/scores.csv")
#omit lines where bali score failed
scores <- na.omit(scores)
#segmentation of balibase references
ref1 <- subset(scores, scores$referid == "ref1")
ref2 <- subset(scores, scores$referid == "ref2")
ref3 <- subset(scores, scores$referid == "ref3")
ref4 <- subset(scores, scores$referid == "ref4")
ref5 <- subset(scores, scores$referid == "ref5")
ref7 <- subset(scores, scores$referid == "ref7")
ref8 <- subset(scores, scores$referid == "ref8")


#comparison of scores for the different references
bwplot(algorithm~sp, data=ref1, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score", main="results ref1")
bwplot(algorithm~sp, data=ref2, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score", main="results ref2")
bwplot(algorithm~sp, data=ref3, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score", main="results ref3")
bwplot(algorithm~sp, data=ref4, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score", main="results ref4")
bwplot(algorithm~sp, data=ref5, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score", main="results ref5")
bwplot(algorithm~sp, data=ref7, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score", main="results ref7")
bwplot(algorithm~sp, data=ref8, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score", main="results ref8")

#comparison of execution times
bwplot(algorithm~time, data=ref1, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución", main="times ref1")
bwplot(algorithm~time, data=ref2, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución", main="times ref2")
bwplot(algorithm~time, data=ref3, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución", main="times ref3")
bwplot(algorithm~time, data=ref4, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución", main="times ref4")
bwplot(algorithm~time, data=ref5, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución", main="times ref5")
bwplot(algorithm~time, data=ref8, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución", main="times ref8")
