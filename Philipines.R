setwd("C:/Users/jntrcs/Desktop/Statistical\ Computation/BurrDistribution")
philip<-read.csv("householdincome.csv")
#income=philip$Total.Household.Income/sd(philip$Total.Household.Income)
plot(density(income))
source("burrFunctions.R")


curve(dburr(x, c, k), xlim=c(0,10), ylim=c(0,1.3), main="Comparison of fitted MLE and observed values", 
      xlab="Household income ")
lines(density(income), col="Blue")
#fake<-rburr(41544, c,k)
#lines(density(fake, from=0, to=25, cut=150), col="red")

income=philip$Total.Household.Income/100000
plot(density(income))
mle<-burr.mle(income, 1)
c.mle<-mle[1]
k.mle<-mle[2]
percentile_est<-percentile_matching(income)
c.pr<-percentile_est[1]
k.pr<-percentile_est[2]
c.pr
k.pr
#bayes.draws<-bayes_estimates_sim(income, .02, nDraws = 7000, nBurn=1000, check=T)
load("BayesResults.RData")
c.bayes<-mean(bayes.draws[,1])
k.bayes<-mean(bayes.draws[,2])
save(bayes.draws, file="BayesResults.RData")
bayes.int.c<-quantile(bayes.draws[,1], c(.025, .975))
bayes.int.k<-quantile(bayes.draws[,2], c(.025,.975))
plot(1, main=expression(1%~~%2))

pdf(file="IncomeModel.pdf")
par(lwd=2)
curve(dburr(x, c.mle, k.mle), xlim=c(0,8), ylim=c(0,.6), main="Comparison of Density Curves using Esimators and Observed Values", 
      xlab=expression(paste("Household Income (",1%~~% 2000," USD)")), ylab="Density", col="Blue")
lines(density(income), col="black")
curve(dburr(x, c.bayes, k.bayes), add=T, col="red", lty=2)
curve(dburr(x, c.pr, k.pr), add=T, col="gold")
legend("topright", lty=c(1,1,2,1), col=c("black", "blue", "red", "gold"), legend = c("Observed Incomes",
                                                                                     "MLE Estimates",
                                                                                     "Bayes Estimates",
                                                                                     "Percentile Matching"))
dev.off()
mle.boostraps<-matrix(0, nrow=1000, ncol=2)
pr.boostraps<-matrix(0, nrow=1000, ncol=2)
for (i in 1:nrow(mle.boostraps)){
  boot<-sample(income, length(income), replace=T)
  mle.boostraps[i,]<-burr.mle(boot,3)
  pr.boostraps[i,]<-percentile_matching(boot)
}
save(mle.boostraps, pr.boostraps, file="Bootstraps.RData")
c.int.mle<-quantile(mle.boostraps[,1],c(.025,.975))
k.int.mle<-quantile(mle.boostraps[,2],c(.025,.975))
c.int.pr<-quantile(pr.boostraps[,1],c(.025,.975))
k.int.pr<-quantile(pr.boostraps[,2],c(.025,.975))
bayes.int.c
bayes.int.k
k.int.mle
