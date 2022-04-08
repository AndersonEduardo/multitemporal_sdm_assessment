######################################################################
####SCRIPT FOR THE STUDY 'Assessing multitemporal calibration for ####
####              species distribution models'                    ####
######################################################################


###Some necessary Parameters###
spsTypes = c('spHW', 'spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
outputData = list() #tabela de dados de saida
vetor.nomes = vector()
projectFolder = getwd() #pasta do projeto


### AUC and TSS of the models

##multitemporal, spHW

spHWmulti = read.csv(paste(projectFolder,'/models/multitemporal/spHW/StatisticalResults-spHW.csv', sep=''), header=TRUE)

##multitemporal, spCD
spCDmulti = read.csv(paste(projectFolder,'/models/multitemporal/spCD/StatisticalResults-spCD.csv', sep=''), header=TRUE)

##monotemporal, spHW
spHWmono = read.csv(paste(projectFolder,'/models/monotemporal/spHW/StatisticalResults-spHW.csv', sep=''), header=TRUE)

##monotemporal, spCD
spCDmono = read.csv(paste(projectFolder,'/models/monotemporal/spCD/StatisticalResults-spCD.csv', sep=''), header=TRUE)


## boxplots models X AUC and TSS, full dataset

jpeg('results/boxplotModelos&Acuracia_dadosTotais.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,1,1), cex=1.3)
boxplot(rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$AUC ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType, ylim=c(0,1), ylab='AUC')
boxplot(rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$TSS ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType, ylim=c(0,1), ylab='TSS')
dev.off()

##tests (both statistically significative)

kruskal.test( rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$AUC ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType )
kruskal.test( rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$TSS ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType )



## boxplots models X AUC and TSS, species

jpeg('results/boxplotModelos&Acuracia_sps.jpeg', height=600)
par(mfrow=c(2,2), las=2, mar=c(8,5,2,1), cex=1.1)
boxplot(rbind(spHWmulti,spHWmono)$AUC ~  rbind(spHWmulti,spHWmono)$modelType, ylim=c(0,1), ylab=c('AUC'), main='HW species')
boxplot(rbind(spCDmulti,spCDmono)$AUC ~  rbind(spCDmulti,spCDmono)$modelType, ylim=c(0,1), ylab=c('AUC'), main='CD species')
boxplot(rbind(spHWmulti,spHWmono)$TSS ~  rbind(spHWmulti,spHWmono)$modelType, ylim=c(0,1), ylab=c('TSS'), main='HW species')
boxplot(rbind(spCDmulti,spCDmono)$TSS ~  rbind(spCDmulti,spCDmono)$modelType, ylim=c(0,1), ylab=c('TSS'), main='CD species')
dev.off()

##tests (all statistically significative)

kruskal.test(rbind(spHWmulti,spHWmono)$AUC ~  rbind(spHWmulti,spHWmono)$modelType)
kruskal.test(rbind(spCDmulti,spCDmono)$AUC ~  rbind(spCDmulti,spCDmono)$modelType)
kruskal.test(rbind(spHWmulti,spHWmono)$TSS ~  rbind(spHWmulti,spHWmono)$modelType)
kruskal.test(rbind(spCDmulti,spCDmono)$TSS ~  rbind(spCDmulti,spCDmono)$modelType)


## boxplots sp X AUC and TSS, full dataset

jpeg('results/boxplotSps&Acuracia.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1),cex=1.1)
boxplot(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sp, ylim=c(0,1), ylab=c('AUC'), main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sp, ylim=c(0,1), ylab=c('AUC'), main='Monotemporal')
boxplot(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sp, ylim=c(0,1), ylab=c('TSS'), main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sp, ylim=c(0,1), ylab=c('TSS'), main='Monotemporal')
dev.off()


##tests (both AUC and TSS statistically significative only for SDMmulti)

kruskal.test(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sp) 
kruskal.test(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sp)
kruskal.test(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sp)
kruskal.test(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sp)


## boxplots sampleSize X AUC aand TSS, full dataset

jpeg('results/boxplotSampleSize&Acuracia_dadosTotais.jpeg')
par(mfrow=c(2,2))
boxplot(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
dev.off()

## boxplots sampleSize X AUC, species

jpeg('results/boxplotSampleSize&Acuracia_spHW.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1))
boxplot(rbind(spHWmulti)$AUC ~  rbind(spHWmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spHWmono)$AUC ~  rbind(spHWmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spHWmulti)$TSS ~  rbind(spHWmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spHWmono)$TSS ~  rbind(spHWmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
title('spHW', outer=TRUE, line=-1)
dev.off()

## boxplots modelo X AUC, species

jpeg('results/boxplotSampleSize&Acuracia_spCD.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1))
boxplot(rbind(spCDmulti)$AUC ~  rbind(spCDmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spCDmono)$AUC ~  rbind(spCDmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spCDmulti)$TSS ~  rbind(spCDmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spCDmono)$TSS ~  rbind(spCDmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
title('spCD', outer=TRUE, line=-1)
dev.off()


### niche overlap

outputData = read.csv(file=paste(projectFolder,'/results/output.csv',sep=''), header=TRUE)

## Schoener e Hellinger for full dataset
jpeg('results/boxplotDadosTotais.jpeg', width=600)
par(mfrow=c(1,2), mar=c(8,3,3,1), cex=1.4, las=2)
boxplot(outputData$Schoeners_D_simi~ outputData$sdmType, ylim=c(0,1), main="Schoeners' D")
boxplot(outputData$Hellinger_I_simi~ outputData$sdmType, ylim=c(0,1), main='Hellinger')
dev.off()

##tests (both statistically non-significative)
kruskal.test(Schoeners_D_simi ~ sdmType, data = outputData)
kruskal.test(Hellinger_I_simi ~ sdmType, data = outputData)


## bosplots for sps
jpeg('results/boxplotSps.jpeg', height=650)
par(mfrow=c(2,2), mar=c(7,4.5,6,1), cex=1.1, las=2)# cex.axis=2.5, cex.lab=3, cex.main=3)
boxplot(outputData[outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spHW',]$sdmType, ylim=c(0,1), ylab="Schoener's D", main='HW species')
boxplot(outputData[outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spHW',]$sdmType, ylim=c(0,1), ylab="Hellinger", main='HW species')
boxplot(outputData[outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spCD',]$sdmType, ylim=c(0,1), ylab="Schoener's D",  main='CD species')
boxplot(outputData[outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spCD',]$sdmType, ylim=c(0,1), ylab="Hellinger", main='CD species')
dev.off()

##tests (all statistically non-significative)
kruskal.test(outputData[outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spHW',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spHW',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spCD',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spCD',]$sdmType)

##Density for full dataset
jpeg('results/densidadeDadosTotais.jpeg', width=600, height = 400)
par(mfrow=c(1,2), lwd=2, cex=1)
plot(density(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi),ylim=c(0,5), lwd=2, col='red', main='', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi),ylim=c(0,5), lwd=2, col='red', main='', xlab='Hellinger', ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi), lwd=2)
##
legend(x='topright', legend=c('Multitemporal calibration', 'Monotemporal calibration'), lty=1, col=c('red','black'), bty='n')
dev.off()

##Density for full dataset
jpeg('results/densidade_sps.jpeg')
par(mfrow=c(2,2), mar=c(5,4,3,1), lwd=2, cex=1)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi),ylim=c(0,7), lwd=2, col='red', main='HW species', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi),ylim=c(0,5), lwd=2, col='red', main='HW species', xlab="Hellinger", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi), lwd=2)
legend(x='topright', legend=c('Multitemporal calibration', 'Monotemporal calibration'), lty=1, col=c('red','black'), bty='n', cex=0.8)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi),ylim=c(0,5), lwd=2, col='red', main='CD species', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi),ylim=c(0,5), col='red', lwd=2, main='CD species', xlab="Hellinger", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi), lwd=2)
dev.off()

##Schoener and Hellinger for time
jpeg('results/Shoener&HellingerXtempo.jpeg',width=600, height=600)
par(mfrow=c(2,2), mar=c(4,4,4,1), cex=1.2)
plot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal',]$kyrBP),type='p',ylab="Schoeners' D", xlab="Time (kyr BP)", ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), main='Multitemporal')
#
plot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Multitemporal')
#
plot(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal',]$kyrBP),type='p',ylab="Schoeners' D",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Monotemporal')
#
plot(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Monotemporal')
dev.off()

## Shoener e Hellinger for tempo - species
jpeg('results/Shoener&HellingerXtempo_sps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='HW species', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)", main='HW species',ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$kyrBP),type='p',ylab="Schoeners' D",xlab="Time (kyr BP)", main='CD species', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)", main='CD species',ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
dev.off()

## Schoener X sample size X time - spHW
jpeg('results/SchoenerXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()

## Schoener X sample size X time - spCD
jpeg('results/SchoenerXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()


## Hellinger X sample size X time - spHW
jpeg('results/HellingerXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()

## Hellinger X sample size X time - spCD
jpeg('results/HellingerXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()



## Sample size - full dataset
jpeg('results/boxplot_sampleSize_dadosTotais.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Sample size - spHW
jpeg('results/boxplot_sampleSize_spHW.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Sample size - spCD
jpeg('results/boxplot_sampleSize_spCD.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Shoener's D e Hellinger X number of time layers - SDMmulti
jpeg('results/boxplot_NumberOfTimeLayers.jpeg', height=1000, width=600)
par(mfrow=c(3,2), cex=1.3)
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='Full dataset')
##
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='Full dataset')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='HW species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='HW species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='CD species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='CD species')
dev.off()


##graphics for clamping

library(raster)

projectFolder = getwd()
sdmTypes = c("multitemporal", "monotemporal")
spsTypes = c("spHW", "spCD")
sampleSizes = c(10) #c(10, 50, 100)
numRep = 2 #5
clampList = list()
territory = list()


for(h in 1:length(sdmTypes)){
  for(i in 1:length(spsTypes)){
    for(m in 1:length(sampleSizes)){
      for(n in 1:numRep){
        for(l in c('000', '006', '012')){ #1:24){

          ##clamping maps
          sdmClampPath = paste(projectFolder,'/models/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[m],'.replicate',n,'/proj_',l,'kyr/','proj_',l,'kyr_ClampingMask.grd',sep='') #path to suitability mapsfrom SDM
          clampLayer_i = raster(sdmClampPath)
          scenName = paste(sdmTypes[h],'_proj_',l,'kyr_',spsTypes[i],'.sample',sampleSizes[m],'.replicate',n,sep='')
          clampList[[scenName]] = clampLayer_i
          clamping = (sum(getValues(clampLayer_i)>0, na.rm=TRUE)/ncell(getValues(clampLayer_i))) * 100

          ##map of species distribution
          sdmDistPath = paste(projectFolder,'/models/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[m],'.replicate',n,'/proj_',l,'kyr/','proj_',l,'kyr_',spsTypes[i],'.sample',sampleSizes[m],'.replicate',n,'_TSSbin.grd',sep='') #path to suitability mapsfrom SDM
          distLayer_i = raster(sdmDistPath)

          ##computing extent of clamping being performed in modelled distribution
          distUnderClamp = (clampLayer_i + distLayer_i)==2
          distUnderClamp = ( freq(distUnderClamp, value=1)/sum(freq(distUnderClamp)[1:2,2]) ) * 100

          ##output table
          outputDF = outputData[ which(outputData$sdmType==sdmTypes[h] & outputData$sp==spsTypes[i] & outputData$sampleSize==sampleSizes[m] & outputData$replicate==n & outputData$kyrBP==as.numeric(l)), ]
          
          if(nrow(outputDF) > 0){
            territory[[scenName]] = data.frame( outputDF,
                                                clamping = clamping,
                                                distUnderClamp = distUnderClamp )

          }
        }
      }
    }
  }
}


##transforming the list into stack of gridfiles
clampStack = stack(clampList)


##multitemporal

## all clamping maps - multitemporal, spHW, sample size 10 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample10.*replicate1', x=names(clampStack), value=TRUE) #separating the names
scenNames =  grep(pattern='^multitemporal.*spHW.*sample100.*replicate1', x=scenNames, value=TRUE, invert=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Multitemporal_SpHW_sample10_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions = colorRampPalette(c("lightgrey","red")),
                     main = 'spHW, sample size = 10',
                     names.attr = c(paste('spHW ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spCD, sample 10 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample10.*replicate1', x=names(clampStack), value=TRUE) #separating the names
scenNames =  grep(pattern='^multitemporal.*spCD.*sample100.*replicate1', x=scenNames, value=TRUE, invert=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Multitemporal_SpCD_sample10_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 10',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spHW, sample 50 pts (all time layers)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample50.*replicate1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Multitemporal_SpHW_sample50_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 50',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spCD, sample 50 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample50.*replicate1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Multitemporal_SpCD_sample50_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 50',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spHW, sample 100 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample100.*replicate1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Multitemporal_SpHW_sample100_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 100',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spCD, sample 100 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample100.*replicate1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Multitemporal_SpCD_sample100_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 100',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()


###monotemporal

## all clamping maps - monotemporal, spHW, sample 10 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample10.*replicate1', x=names(clampStack), value=TRUE) #separating the names
scenNames =  grep(pattern='^monotemporal.*spHW.*sample100.*replicate1', x=scenNames, value=TRUE, invert=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Monotemporal_SpHW_sample10_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 10',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spCD, sample 10 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample10.*replicate1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Monotemporal_SpCD_sample10_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 10',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spHW, sample 50 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample50.*replicate1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Monotemporal_SpHW_sample50_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 50',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spCD, sample 50 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample50.*replicate1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Monotemporal_SpCD_sample50_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 50',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spHW, sample 100 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample100.*replicate1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Monotemporal_SpHW_sample100_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 100',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spCD, sample 100 pts (all temporal layers)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample100.*replicate1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('results/clamp_Monotemporal_SpCD_sample100_replicate1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 100',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()


### trend graphs ###


##dataClamp = data.frame(clamping=as.numeric(territory), scenario=names(territory))
dataClamp = do.call('rbind', territory)
dataClamp$scenario = names(territory)
rownames(dataClamp) = NULL

######
## dataClamp$sdm = NA
## dataClamp[ grep('monotemporal', dataClamp$scenario), ]$sdm = 'monotemporal'
## dataClamp[ grep('multitemporal', dataClamp$scenario), ]$sdm = 'multitemporal'
## dataClamp$kyr = c(0:23)
## dataClamp$sp = NA
## dataClamp[grep('spHW', dataClamp$scenario),]$sp = 'spHW'
## dataClamp[grep('spCD', dataClamp$scenario),]$sp = 'spCD'
## dataClamp$sample = NA
## dataClamp[grep('10.replica', dataClamp$scenario),]$sample = 10
## dataClamp[grep('50.replica', dataClamp$scenario),]$sample = 50
## dataClamp[grep('100.replica', dataClamp$scenario),]$sample = 100
######


## clamping for America do Sul X Schoener's D (D = overlap between real X modelled niche-based distribution) ##

jpeg('results/clampAmSulXschoenerXsampleXsps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping), ylim=c(0,1), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping), ylim=c(0,1), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spCD multitemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping), ylim=c(0,1), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping), ylim=c(0,1), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spHW monotemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spHW multitemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
dev.off()



## clamping across modelled distribution X Schoener's D (D = overlap between real X modelled niche-based distribution) ##

jpeg('results/clampSpsDistXschoenerXsampleXsps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3,1), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spCD multitemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>=0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>=0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,0), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, ylab="Schoener's D", col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW monotemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-10,5), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW multitemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-10,0), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
dev.off()



## clamping for America do Sul X clamping across modelled distributions (D = overlap between real X modelled niche-based distribution) ##

jpeg('results/clampSpsDistXclampAmSul.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3.5,1), ylim=c(-2,3), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spCD multitemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,-0.5), ylim=c(-3,0.5), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW monotemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3.5,1), ylim=c(-2,3), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW multitemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,-0.5), ylim=c(-3,0.5), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
dev.off()



## clamping trends across geological time


jpeg('results/clampXtempo.jpeg', width=900)
par(mfrow=c(1,2))
plot(dataClamp[which(dataClamp$sdm=='multitemporal'),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal'),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='South América area (in %)', main='Multitemporal')
plot(dataClamp[which(dataClamp$sdm=='monotemporal'),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal'),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='South América area (in %)', main='Monotemporal')
dev.off()




## clamping trends across geological time, accounting for species, sample size e model calibration (mono and multitemporal)

## graphics for spCD
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.7, cex.axis=2, cex.main=2)
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['mono']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['multi']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='South América area (in %)', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
dev.off()


## graphics for spHW
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.7, cex.axis=2, cex.main=2)
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['mono']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['multi']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='South América area (in %)', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
dev.off()




## boxplot clamping (full dataset and sps)
jpeg('results/BoxplotClamp.jpeg', width=1200)
par(mfrow=c(1,3), cex=1.3)
boxplot(log(dataClamp$clamping) ~ dataClamp$sdm, ylab='log(South America area (in %))', main='Full dataset')
##
boxplot(log(dataClamp[which(dataClamp$sp == 'spHW'),]$clamping) ~ dataClamp[which(dataClamp$sp == 'spHW'),]$sdm, ylab='log(South America area (in %))', main='HW species')
##
boxplot(log(dataClamp[which(dataClamp$sp == 'spCD'),]$clamping) ~ dataClamp[which(dataClamp$sp == 'spCD'),]$sdm, ylab='log(South America area (in %))', main='CD species')
dev.off()



##correlation between Schoener and clamping across sps' modelled distributions

##monotemporal - spCD
corSpCDmono10 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpCDmono50 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpCDmono100 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##monotemporal - spHW
corSpHWmono10 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpHWmono50 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpHWmono100 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##multitemporal - spCD
corSpCDmulti10 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpCDmulti50 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpCDmulti100 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                           dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##multitemporal - spHW
corSpHWmulti10 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpHWmulti50 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpHWmulti100 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                           dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')



##table
corTable = data.frame( scenario = c('corFullDataset','corMulti','corMono','corSpCDmono10','corSpCDmono50','corSpCDmono100','corSpHWmono10','corSpHWmono50','corSpHWmono100','corSpCDmulti10','corSpCDmulti50','corSpCDmulti100','corSpHWmulti10','corSpHWmulti50','corSpHWmulti100'),
                       correlation = c(as.numeric(corFullDataset$estimate),as.numeric(corMulti$estimate),as.numeric(corMono$estimate),as.numeric(corSpCDmono10$estimate),as.numeric(corSpCDmono50$estimate),as.numeric(corSpCDmono100$estimate),as.numeric(corSpHWmono10$estimate),as.numeric(corSpHWmono50$estimate),as.numeric(corSpHWmono100$estimate),as.numeric(corSpCDmulti10$estimate),as.numeric(corSpCDmulti50$estimate),as.numeric(corSpCDmulti100$estimate),as.numeric(corSpHWmulti10$estimate),as.numeric(corSpHWmulti50$estimate),as.numeric(corSpHWmulti100$estimate)),
                       p.value = c(as.numeric(corFullDataset$p.value),as.numeric(corMulti$p.value),as.numeric(corMono$p.value),as.numeric(corSpCDmono10$p.value),as.numeric(corSpCDmono50$p.value),as.numeric(corSpCDmono100$p.value),as.numeric(corSpHWmono10$p.value),as.numeric(corSpHWmono50$p.value),as.numeric(corSpHWmono100$p.value),as.numeric(corSpCDmulti10$p.value),as.numeric(corSpCDmulti50$p.value),as.numeric(corSpCDmulti100$p.value),as.numeric(corSpHWmulti10$p.value),as.numeric(corSpHWmulti50$p.value),as.numeric(corSpHWmulti100$p.value)) )

corTable[,'p.value'] = round(corTable[,'p.value'], 3)

write.csv(corTable, paste(projectFolder,'correlationTable.csv'), row.names=FALSE)



### maps to compare climatic conditions between 0 and 22 kyrBP

temp0kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_01.asc')
temp0kyr = mask(temp0kyr, AmSulShape)
preci0kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_12.asc')
preci0kyr = mask(preci0kyr, AmSulShape)

temp22kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/022/bioclim_01.asc')
temp22kyr = mask(temp22kyr, AmSulShape)
preci22kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/022/bioclim_12.asc')
preci22kyr = mask(preci22kyr, AmSulShape)

jpeg('results/temp0&22kyr.jpeg', width=800)
par(mfrow=c(1,2), mar=c(5,5,5,6))
plot(temp0kyr, main='0 kyr BP'); grid()
plot(temp22kyr, main='22 kyr BP'); grid()
dev.off()

jpeg('results/preci0&22kyr.jpeg', width=800)
par(mfrow=c(1,2), mar=c(5,5,5,6))
plot(preci0kyr, main='0 kyr BP'); grid()
plot(preci22kyr, main='22 kyr BP'); grid()
dev.off()

clamp0 = mask(clampStack$multitemporal_proj_0kyr_spHW.sample50.replica1, AmSulShape)
clamp22 = mask(clampStack$multitemporal_proj_22kyr_spHW.sample50.replica1, AmSulShape)

par(mfrow=c(1,2))
plot(clamp0, main=c(paste('spHW ',0,'kyr BP',sep='')), col=c('lightgrey','red'), legend=FALSE)
plot(clamp22, main=c(paste('spHW ',22,'kyr BP',sep='')), col=c('lightgrey','red'), legend=FALSE)

dev.off()



##mapped distributions for present, interglacial and maximum glacial

library(raster)
library(maptools)
library(RColorBrewer)
AmSulShape = readShapePoly("/home/anderson/shapefiles/Am_Sul/borders.shp")

### MULTITEMPORAL ###

## spHW ##

## real distributions
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.2
HW22Real = raster(paste(projectFolder,'NichoReal/spHW/022.asc',sep='')) > 0.2

## multitemporal
HWModel_0kyrSample10 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample10.replica1/proj_0kyr/proj_0kyr_spHW.sample10.replica1_TSSbin.grd',sep=''))
HWModel_0kyrSample50 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample50.replica1/proj_0kyr/proj_0kyr_spHW.sample50.replica1_TSSbin.grd',sep=''))
HWModel_0kyrSample100 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample100.replica1/proj_0kyr/proj_0kyr_spHW.sample100.replica1_TSSbin.grd',sep=''))
##
HWModel_22kyrSample10 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample10.replica1/proj_22kyr/proj_22kyr_spHW.sample10.replica1_TSSbin.grd',sep=''))
HWModel_22kyrSample50 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample50.replica1/proj_22kyr/proj_22kyr_spHW.sample50.replica1_TSSbin.grd',sep=''))
HWModel_22kyrSample100 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample100.replica1/proj_22kyr/proj_22kyr_spHW.sample100.replica1_TSSbin.grd',sep=''))

## spCD ##

## real distributions
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.2
CD22Real = raster(paste(projectFolder,'NichoReal/spCD/022.asc',sep='')) > 0.2

## SDM multitemporal
CDModel_0kyrSample10 = raster(paste(projectFolder,'maxent/multitemporal/spCD/spCD.sample10.replica1/proj_0kyr/proj_0kyr_spCD.sample10.replica1_TSSbin.grd', sep=''))
CDModel_0kyrSample50 = raster(paste(projectFolder,'maxent/multitemporal/spCD/spCD.sample50.replica1/proj_0kyr/proj_0kyr_spCD.sample50.replica1_TSSbin.grd', sep=''))
CDModel_0kyrSample100 = raster(paste(projectFolder,'/maxent/multitemporal/spCD/spCD.sample100.replica1/proj_0kyr/proj_0kyr_spCD.sample100.replica1_TSSbin.grd', sep=''))
##
CDModel_22kyrSample10 = raster(paste(projectFolder,'/maxent/multitemporal/spCD/spCD.sample10.replica1/proj_22kyr/proj_22kyr_spCD.sample10.replica1_TSSbin.grd', sep='')) 
CDModel_22kyrSample50 = raster(paste(projectFolder,'/maxent/multitemporal/spCD/spCD.sample50.replica1/proj_22kyr/proj_22kyr_spCD.sample50.replica1_TSSbin.grd', sep=''))
CDModel_22kyrSample100 = raster(paste(projectFolder,'/maxent/multitemporal/spCD/spCD.sample100.replica1/proj_22kyr/proj_22kyr_spCD.sample100.replica1_TSSbin.grd', sep=''))

##overlaps spHW

jpeg(filename='results/sobreposicoesHWmulti.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(HWcurrentReal*1+HWModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legend
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(HW22Real*1+HWModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('HW species',outer=TRUE,cex=4)
dev.off()

##overlap spCD

jpeg(filename='results/sobreposicoesCDmulti.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(CDcurrentReal*1+CDModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legend
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(CD22Real*1+CDModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('CD species',outer=TRUE,cex=4)
dev.off()


### MONOTEMPORAL ###


## spHW ##

## real distributions
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.2
HW22Real = raster(paste(projectFolder,'NichoReal/spHW/022.asc',sep='')) > 0.2

## monotemporal
HWModel_0kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample10.replica1/proj_0kyr/proj_0kyr_spHW.sample10.replica1_TSSbin.grd')
HWModel_0kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample50.replica1/proj_0kyr/proj_0kyr_spHW.sample50.replica1_TSSbin.grd')
HWModel_0kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample100.replica1/proj_0kyr/proj_0kyr_spHW.sample100.replica1_TSSbin.grd')
##
HWModel_22kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample10.replica1/proj_22kyr/proj_22kyr_spHW.sample10.replica1_TSSbin.grd') 
HWModel_22kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample50.replica1/proj_22kyr/proj_22kyr_spHW.sample50.replica1_TSSbin.grd')
HWModel_22kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample100.replica1/proj_22kyr/proj_22kyr_spHW.sample100.replica1_TSSbin.grd')

## spCD ##

## real distributions
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.2
CD22Real = raster(paste(projectFolder,'NichoReal/spCD/022.asc',sep='')) > 0.2

## SDM monotemporal
CDModel_0kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample10.replica1/proj_0kyr/proj_0kyr_spCD.sample10.replica1_TSSbin.grd')
CDModel_0kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample50.replica1/proj_0kyr/proj_0kyr_spCD.sample50.replica1_TSSbin.grd')
CDModel_0kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample100.replica1/proj_0kyr/proj_0kyr_spCD.sample100.replica1_TSSbin.grd')
##
CDModel_22kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample10.replica1/proj_22kyr/proj_22kyr_spCD.sample10.replica1_TSSbin.grd') 
CDModel_22kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample50.replica1/proj_22kyr/proj_22kyr_spCD.sample50.replica1_TSSbin.grd')
CDModel_22kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample100.replica1/proj_22kyr/proj_22kyr_spCD.sample100.replica1_TSSbin.grd')

##overlaps spHW
jpeg(filename='results/sobreposicoesHWmono.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(HWcurrentReal*1+HWModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legend
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(HW22Real*1+HWModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('HW species.',outer=TRUE,cex=4)
dev.off()

##overlap spCD
jpeg(filename='results/sobreposicoesCDmono.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(CDcurrentReal*1+CDModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legend
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(CD22Real*1+CDModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('CD species',outer=TRUE,cex=4)
dev.off()