

dburr<- function(x,c,k){
  c*k*x^(c-1)/(1+x^c)^(k+1)
}


meanBurr<-function(c,k){
  k*beta(k-1/c, 1+1/c)
}

varBurr<-function(c,k){
  k*beta((c*k-2)/c, (c+2)/c)-meanBurr(c, k)^2
}



burr_nth_moment<-function(n, c,k){
  k*beta((c*k-n)/c,(c+n)/c)
}

pburr<-function(x, c, k){
  1-(1+x^c)^-k
}



finv_burr<-function(y,c,k){
  ((1/(1-y))^(1/k)-1)^(1/c)
}

rburr<-function(n, c, k){
  finv_burr(runif(n),c,k)
}

burr_gradient<-function(data, c, k){
  n<-length(data)
  dc<-n/c + sum(log(data))-(k+1)*sum(data^c*log(data)/(data^c+1))
  dk<-n/k - sum(log(data^c+1))
  c(dc, dk)
}

burr_hessian<-function(data, c,k){
  n<-length(data)
  xc<-data^c
  logx<-log(data)
  dc2<--n/c^2 -(k-1)*sum(xc*logx^2/(xc+1)^2)
  dcdk<--sum(xc*logx/(xc+1))
  dk2<--n/k^2
  matrix(c(dc2, dcdk, dcdk, dk2), byrow=2, nrow=2, ncol=2)
}



burr_c_deriv<-function(data, c, k){
  n<-length(data)
  dc<-n/c + sum(log(data))+(-k-1)*sum(data^c*log(data)/(data^c+1))
  dc
}

burr_c_2nd_deriv<-function(data, c, k){
  n<-length(data)
  xc<-data^c
  logx<-log(data)
  dc2<--n/c^2 -(k-1)*sum(xc*logx^2/(xc+1)^2)
  dc2
}


burr.mle<-function(data, starting_c){
  #handle invalid input
  if (starting_c<0){
    stop("Starting value must be greater than zero.")
  }
  if (any(data<0)){
    stop("Data values must be greater than 0")
  }
  
  
  tol<-.0001
  c<-starting_c
  k<-length(data)/sum(log(data^c+1))
  #print(abs(burr_c_deriv(data,c, k)))
  while(abs(burr_c_deriv(data,c, k))>tol){
    c<- c-burr_c_deriv(data, c, k)/ burr_c_2nd_deriv(data, c,k)
    k<-length(data)/sum(log(data^c+1))
  }
  return(c(c_mle=c, k_mle=k))
}

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}



percentile_matching<-function(data){
  med<-median(data)
  s<-quantile(data,.75)
  
  c_as_function<-function(c, med, s){
    (4^(log(med^c+1)/log(2))-1)^(1/c)-s
  }
  vari=1
  tryCatch(c<-uniroot(c_as_function, interval=c(.5,30), med=med,s=s)$root,
           error = function(c) {print("No root available")
             return(c(NA,NA))
             vari=5}
           
           
  )
  
  if(is.numeric(c)){
    k=log(4)/log(s^c+1)
    return(c(c_est=c, k_est=k))
  }
  return(c(NA,NA))
}

#check allows you to visually examine your mixing, but it will be set to false for the simulation
bayes_estimates_sim<-function(data, d, nDraws=2000, nBurn=500, check=FALSE){
  ##c and k will follow a gamma prior with alpha=2.5 and beta =.5
  sumlogdata<-sum(log(data))
  n<-length(data)
  log_posterior<-function(params){
    if (any(params<0))return(-Inf)
      c<-params[1]
      k<-params[2]
      n*log(c)+n*log(k)+(c-1)*sumlogdata-(k+1)*sum(log(1+data^c))+log(c)-.5*c+log(k)-.5*k
  }
  draws<-matrix(0, ncol=2, nrow=nDraws+nBurn)
  draws[1,]<-c(5,5)
  old<-log_posterior(draws[1,])
  accept<-0
  for (i in 2:nrow(draws)){
    cands<-draws[i-1,]+rnorm(2, 0, d)#
    new<-log_posterior(cands)
    if (log(runif(1))<new-old){
      draws[i,]<-cands
      accept<-accept+1
      old<-new
    }else{draws[i,]<-draws[i-1,]}
  }
  draws<-draws[-(1:nBurn),]
  if (check){
    print(accept/(nDraws+nBurn))
    plot(draws[,1])
    plot(draws[,2])
    plot(draws[,1],draws[,2])
    return(draws)
  }

  if (accept/(nDraws+nBurn)>.8 |accept/(nDraws+nBurn)<.05)return(NA) #something went wrong and results should not be trusted
  m<-apply(draws, 2, mean)
  return(c(c=m[1], k=m[2]))
}


