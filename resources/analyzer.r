library(lattice)
library(reshape)
library(Hmisc)
library(car)

scores <- read.csv("/Users/guilespi/Documents/Development/interrupted/bioinformatics/bio1/multiple_seq_aligner/results/scores.csv")
#omit lines where bali score failed
scores <- na.omit(scores)
#segmentation of 
ref1 <- subset(scores, scores$referid == "ref1")
ref2 <- subset(scores, scores$referid == "ref2")
ref3 <- subset(scores, scores$referid == "ref3")
ref4 <- subset(scores, scores$referid == "ref4")
ref5 <- subset(scores, scores$referid == "ref5")

par(mfrow=c(1,3))
bwplot(tc ~ algorithm, data = ref2)
#comparison of scores for the different references
bwplot(algorithm~sp, data=ref1, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score")
bwplot(algorithm~sp, data=ref2, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score")
bwplot(algorithm~sp, data=ref3, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score")
bwplot(algorithm~sp, data=ref4, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score")
bwplot(algorithm~sp, data=ref5, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="score")
#comparison of execution times
bwplot(algorithm~time, data=ref1, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución")
bwplot(algorithm~time, data=ref2, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución")
bwplot(algorithm~time, data=ref3, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución")
bwplot(algorithm~time, data=ref4, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución")
bwplot(algorithm~time, data=ref5, panel=panel.bpplot, datadensity=TRUE,ylab="algoritmo", xlab="tiempo de ejecución")


scatterplot(sp~tc, data=ref1)

qqPlot(ref1$time, main="Normal QQ plot de temperatura")
ref1.poa = subset(ref1, ref1$algorithm == "poa")
plot(density(ref1.poa$sp, na.rm=T))
ref1.clustalw2  = subset(ref1, ref1$algorithm == "clustalw2")
lines(density(ref1.clustalw2$sp, na.rm=T))
