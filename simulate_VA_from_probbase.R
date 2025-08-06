## Simulate VA example
library(mnormt)

# n: number of simulations
# pi: disease frequency
# pb: probbase
#simulate_va=function(n,pi,pb,blocks,covs)

n=5

pi=c(0.3,0.5,0.2)

pb=matrix(runif(length(pi)*10),10,length(pi))

blocks=list(c(1:3),c(4:6),c(7:10))

dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc 

covs=list(dcov(3,0.1),dcov(3,0.1),dcov(4,0.1))





# Number of causes of death
n_cod=length(pi)

# Number of questions
n_q=nrow(pb)

# Number of blocks
n_block=length(blocks)


# Firstly simulate under conditional independence
va=matrix(NA,n,n_q)
for (i in 1:n) {
  cause_of_death=sample(1:n_cod,1,prob=pi)
  
  qprobs=pb[,cause_of_death]
  
  answers=rep(NA,n_q)
  
  for (b in 1:n_block) {
    wq=blocks[[b]] # indices of questions in block b
    bcov=covs[[b]] # covariance matrix for latent variable
    nb=length(wq) # number of questions in block
    
    latent_q=rmnorm(1,mean=rep(0,nb),varcov=bcov)
    
    bprobs=qprobs[wq]
    thresholds=qnorm(1-bprobs)
    
    answers[wq]=latent_q>thresholds
  }
  va[i,]=answers
}


