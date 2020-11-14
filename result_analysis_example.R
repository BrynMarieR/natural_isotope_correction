skewed_results <- read.csv(file="./skewed.csv")
classical_results <- read.csv(file="./classical.csv")
input <- read.csv(file="./samplefile.csv")

par(mfrow=c(2,2))
expt1.mat <- as.matrix(rbind(input[1,3:6], classical_results[1,3:6], skewed_results[1,3:6]))
barplot(expt1.mat,beside=T,col=c("black","red","blue"),main="Experiment 1",
        xlab="Mass peak",ylab="Fractional abundance")
legend("topright",legend=c("Obs","Classical","Skewed"),fill=c("black","red","blue"))

expt2.mat <- as.matrix(rbind(input[2,3:6], classical_results[2,3:6], skewed_results[2,3:6]))
barplot(expt2.mat,beside=T,col=c("black","red","blue"),main="Experiment 2",
        xlab="Mass peak",ylab="Fractional abundance")
legend("topright",legend=c("Obs","Classical","Skewed"),fill=c("black","red","blue"))

perc.change.expt1 <- apply(expt1.mat, 1, function(x) { x / expt1.mat[1,] })
barplot(t(perc.change.expt1)[c(2,3),],beside=T,col=c("red","blue"),ylim=c(0,2),
        main="Fold change from observed peak in experiment 1",
        xlab="Mass peak",ylab="Fold change")
abline(h=1,col="black")
legend("topright",legend=c("Classical","Skewed"),fill=c("red","blue"))

perc.change.expt2 <- apply(expt2.mat, 1, function(x) { x / expt2.mat[1,] })
barplot(t(perc.change.expt2)[c(2,3),],beside=T,col=c('red','blue'),ylim=c(0,2),
        main="Fold change from observed peak in experiment 2",
        xlab="Mass peak",ylab="Fold change")
abline(h=1,col="black")
legend("topright",legend=c("Classical","Skewed"),fill=c("red","blue"))
