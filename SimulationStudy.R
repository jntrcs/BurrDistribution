setwd("C:/Users/jntrcs/Desktop/Statistical\ Computation/DistributionProject")
source("burrFunctions.R")

parameter_cases<-matrix(c(2,8,7,3, 5,5), byrow = T, nrow=3) #k bigger than c, c bigger than k, both equal 
colnames(parameter_cases)<-c("c", "k")

n.data<-c(20,100,1000)

mh_step_size<-matrix(c(1.6, .6, .2,
                       1.9, 1, .3,
                       2.4, 1.2,.3), byrow = T, nrow = 3) #To get good sampling, we need an appropriate step size to be 
#passed to our Metropolis Hastings algorithm. These were found empirically 

sim.results<-list()
sim.results$bias<-list()
sim.results$mse<-list()
sim.results$bias$c<-list()
sim.results$bias$k<-list()
sim.results$mse$c<-list()
sim.results$mse$k<-list()
sim.results$bias$c$scenario1<-matrix(0, nrow=3, ncol=3)
sim.results$bias$c$scenario2<-matrix(0, nrow=3, ncol=3)
sim.results$bias$c$scenario3<-matrix(0, nrow=3, ncol=3)
sim.results$bias$k$scenario1<-matrix(0, nrow=3, ncol=3)
sim.results$bias$k$scenario2<-matrix(0, nrow=3, ncol=3)
sim.results$bias$k$scenario3<-matrix(0, nrow=3, ncol=3)
sim.results$mse$c$scenario1<-matrix(0, nrow=3, ncol=3)
sim.results$mse$c$scenario2<-matrix(0, nrow=3, ncol=3)
sim.results$mse$c$scenario3<-matrix(0, nrow=3, ncol=3)
sim.results$mse$k$scenario1<-matrix(0, nrow=3, ncol=3)
sim.results$mse$k$scenario2<-matrix(0, nrow=3, ncol=3)
sim.results$mse$k$scenario3<-matrix(0, nrow=3, ncol=3)



for (i in 1:length(n.data)){
  for (j in 1:nrow(parameter_cases)){
    nreps<-1000
    mle.c<-rep(0, nreps)
    mle.k<-rep(0, nreps)
    pr.c<-rep(0, nreps)
    pr.k<-rep(0,nreps)
    bayes.c<-rep(0, nreps)
    bayes.k<-rep(0, nreps)
    for (reps in 1:nreps){
      dat<-rburr(n.data[i], parameter_cases[j, 1], parameter_cases[j,2])
      mle<- burr.mle(dat, 2)
      bayes<-bayes_estimates_sim(dat, mh_step_size[j,i])
      perc<-percentile_matching(dat)
      mle.c[reps]<-mle[1]
      mle.k[reps]<-mle[2]
      pr.c[reps]<-perc[1]
      pr.k[reps]<-perc[2]
      bayes.c[reps]<-bayes[1]
      bayes.k[reps]<-bayes[2]
    }
    sim.results$mse$c[[j]][1,i]<-mean((mle.c-parameter_cases[j, 1])^2)
    sim.results$mse$k[[j]][1,i]<-mean((mle.k-parameter_cases[j, 2])^2)
    sim.results$mse$c[[j]][2,i]<-mean((pr.c-parameter_cases[j, 1])^2)
    sim.results$mse$k[[j]][2,i]<-mean((pr.k-parameter_cases[j, 2])^2)
    sim.results$mse$c[[j]][3,i]<-mean((bayes.c-parameter_cases[j, 1])^2)
    sim.results$mse$k[[j]][3,i]<-mean((bayes.k-parameter_cases[j, 2])^2)
    sim.results$bias$c[[j]][1,i]<-mean(mle.c)-parameter_cases[j,1]
    sim.results$bias$k[[j]][1,i]<-mean(mle.k)-parameter_cases[j,2]
    sim.results$bias$c[[j]][2,i]<-mean(pr.c)-parameter_cases[j,1]
    sim.results$bias$k[[j]][2,i]<-mean(pr.k)-parameter_cases[j,2]
    sim.results$bias$c[[j]][3,i]<-mean(bayes.c)-parameter_cases[j,1]
    sim.results$bias$k[[j]][3,i]<-mean(bayes.k)-parameter_cases[j,2]
    
  }
}

sim.results
save(sim.results, file="Simulation-Results.RData")
