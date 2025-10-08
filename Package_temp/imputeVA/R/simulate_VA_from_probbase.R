##' Simulate a set of VA answers consistent with a given probbase, prior on CoDs, conditional dependency structure, and dependency matrices. 
##' 
##' @name simulatedVA
##' @description This function generates simulated completed verbal autopsy 
##'  questionnaires. It does this assuming a given probbase is 'correct'; that 
##'  is, gives the correct probability of answering a particular question given
##'  a particular cause of death.
##'  
##' It also assumes a dependency structure. We recommend using the preset 
##'  dependency structure, or learning a dependency structure from an existing
##'  set of VA answers using the function learn_VA_dependency. 
##' 
##' Technically, given a block k of answers b_1,b_2,...,b_i, which are taken to 
##'  be conditionally independent of other answers given CoD, we simulate 
##'  answers according to the hierarchical model:
##'  
##'  c~ Multinomial(pi)   # Cause of death
##'  z_1,z_2,... z_i ~ MultivariateNormal(0,V_{k,c})  # 'Latent' answers
##'  b_j=1{z_j > qnorm(p_{c,j})}
##' 
##' where pi is a prior on cause of death (CoD) frequencies, V_{k,c} is a 
##'  variance matrix associated with block k and CoD c, and p_{c,j}, coming 
##'  from a probbase, describes the probability of answering question j as 'yes'
##'  given CoD c. 
##'  
##' The input to our simulation is a probbase, a prior, and a list of sub-lists
##'  V where V[[k]][[c]] is the covariance matrix for block k and CoD c. 
##'   
##' @param n number of simulated VA answers
##' @param probbase A probbase; defaults to null in which case a standard probbase is used. Must be of the same form as data(probbase).
##' @param letterprob Correspondences between probbase letters and probabilities. Defaults to null, in which case a standard lookup is used.
##' @param prior A prior on cause-of-death types; defaults to null in which case a standardised prior is used.
##' @param V A list of lists of variance matrices as above; defaults to null in which case a standardised list is used.
##' @return InterVA output
##' @export
##' @examples
##' 
simulatedVA=function(n,probbase=NULL, letterprob=NULL,prior=NULL,V=NULL) {
  
  # Establish defaults for missing inputs
  if (is.null(probbase)) {
    
  }
  
}


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


