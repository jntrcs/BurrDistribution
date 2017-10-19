##Model mispecification
alpha<-1.5
beta<-3
mis_spec<-rgamma(100, alpha, beta)
mle<-burr.mle(mis_spec, 2)
pr<-percentile_matching(mis_spec)
bayes<-bayes_estimates_sim(mis_spec, d=.35, check=T)
bayes.est<-apply(bayes,2,mean)

pdf(file="Mispecification-Model.pdf")
curve(dgamma(x, alpha,beta), xlim=c(0,3), ylim=c(0,1.8), main="Misspecified Model", ylab="Density", lty=2)
lines(density(mis_spec))
curve(dburr(x, mle[1],mle[2]), add=T, col="blue")
curve(dburr(x, pr[1], pr[2]), add=T, col="gold")
curve(dburr(x, bayes.est[1], bayes.est[2]), col="red", add=T)
legend("topright", lty=c(2,1,1,1,1), lwd=c(1,1,1,1,1), col=c("black","black", "blue", "red", "gold"), 
                                                                                    legend = c("Theoretical Density",
                                                                                                      "Observed Random Draws",
                                                                                                     "MLE Estimates",
                                                                                                     "Bayes Estimates",
                                                                                                     "Percentile Matching"))
dev.off()

alpha<-5
beta<-2
mis_spec<-rgamma(100, alpha, beta)
mle<-burr.mle(mis_spec, 2)
pr<-percentile_matching(mis_spec)
bayes<-bayes_estimates_sim(mis_spec, d=.15, check=T)
bayes.est<-apply(bayes,2,mean)

pdf(file="Misspec2.pdf")
curve(dgamma(x, alpha,beta), xlim=c(0,6), ylim=c(0,.7), main="Misspecified Model (X~Gamma(5,2))", ylab="Density", lty=2)
lines(density(mis_spec))
curve(dburr(x, mle[1],mle[2]), add=T, col="blue")
#curve(dburr(x, pr[1], pr[2]), add=T, col="gold")
curve(dburr(x, bayes.est[1], bayes.est[2]), col="red", add=T)
legend("topright", lty=c(2,1,1,1), lwd=c(1,1,1,1), col=c("black","black", "blue", "red"), 
       legend = c("Theoretical Density",
                  "Observed Random Draws",
                  "MLE Estimates",
                  "Bayes Estimates"
                  ))
dev.off()