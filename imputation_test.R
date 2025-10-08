##**********************************************************************
##  Calibration of marginal probabilities in verbal autopsy         ####
##   Test of imputation proposal
##  Nathan Higgins, James Liley, Eilidh Cowan                                                         
##  2025                                                              
##**********************************************************************


##**********************************************************************
##  Functions, packages, scripts, data                              ####
##**********************************************************************

# Packages
library(openVA)           # Standard VA suite
library(nbc4va)           # Naive Bayes VA
library(mnormt)           # Multivariate normal distribution
library(SPARRAfairness)   # ROC curves
library(latex2exp)        # Plotting

# Main functions: to be incorporated in an R package in time
source("Package_temp/imputeVA/R/simulations.R")

# Datasets
data(RandomPhysician) # From openVA package: WHO2012 VAs with physician-coded CoDs
load("Package_temp/imputeVA/data/standard_simulation.RData") # Standard probbase, prior etc
load("Package_temp/imputeVA/data/cod_subcategories.RData") # Correspondence between CoDs in RandomPhysician and probbase



##**********************************************************************
##  Switches and settings                                           ####
##**********************************************************************

# Number of simulations
nsim=10000 # run with 1500 or 10000

# Degrees of pertubation
ptb=0.2 # Perturb probbase entries by this much (SD)

# Manner of pertubation: uniform or Gaussian
uniform_pertubation = TRUE

# Deviation to use for approximation of partial derivative
epsilon = 0.01

# Perturb this many elements
nptb=500

# Perturb only elements with high probbase entries or not
perturb_all=FALSE

# Random seed
seed=3726152
set.seed(seed)

# Draw plots every so often while running?
running_figures=TRUE

# Save plots as PDFs
save_pdf=TRUE

##**********************************************************************
##  Load data if already simulated                                  ####
##**********************************************************************

lx=list.files("../Images/",pattern=paste0("imputation_test_temp_",nsim,"*"),full.names=TRUE)
n_done=length(lx)

save_file=paste0("../Images/imputation_test_temp_",nsim,"_",1+n_done,".Rdata")
set.seed(seed + n_done)

save.image(save_file)

##**********************************************************************
##  Learn parameters for simulation                                 ####
##**********************************************************************


code=RandomPhysician$code1 # Physician coded causes of death
va=RandomPhysician[,1:246] # VA questionnaire data

# Learn covariance structure from data
rparam=learn_parameters_from_va(va,code)
probbase=standard_simulation$probbase

# Set cause-of-death specific covariance matrices.
V_rp=list()
for (i in 1:ncol(standard_simulation$probbase)) {
  xi=which(cod_subcategories[,1]==colnames(standard_simulation$probbase)[i])
  subcat=cod_subcategories[xi,2]
  if (subcat!="NA") V_rp[[i]]=rparam$V[[subcat]] else V_rp[[i]]=rparam$V0
}
names(V_rp)=colnames(probbase)

# Conditional independence blocks
cblocks=rparam$blocks



##**********************************************************************
##  Simulate VA questionnaires under given model                    ####
##**********************************************************************

# Simulate a VA database and compute posteriors
vax=simulateVA_WHO2012(nsim,seed=32212,V=V_rp)
rva=vax$VA

# Params
blocks=rparam$blocks
va=rva 
probbase_numeric=vax$parameters$probbase_numeric
prior=vax$parameters$prior
wf=which(is.finite(probbase_numeric))



##**********************************************************************
##  Evaluate baseline imputation accuracy (log-loss)                ####
##**********************************************************************


# With real probbase
i_real=imputation_accuracy(va,probbase_numeric,prior,blocks)


##**********************************************************************
##  Evaluate imputation accuracy for very perturbed probbases       ####
##**********************************************************************


# With random probbase
probbase_rand=probbase_numeric
probbase_rand[wf]=runif(length(wf))
i_rand=imputation_accuracy(va,probbase_rand,prior,blocks)

# With scrambled probbase
probbase_scrambled=probbase_numeric
probbase_scrambled[wf]=sample(probbase_numeric[wf])
i_scrambled=imputation_accuracy(va,probbase_scrambled,prior,blocks)

# With globally perturbed probbase
probbase_perturbed=probbase_numeric
probbase_perturbed[wf]=pmax(pmin(probbase_numeric[wf]+rnorm(wf,sd=0.05),1),0)
i_perturbed=imputation_accuracy(va,probbase_perturbed,prior,blocks)

# With slightly perturbed probbase
probbase_almost=probbase_numeric
w2=sample(wf,20)
probbase_almost[w2]=pmax(pmin(probbase_numeric[w2]+rnorm(w2,sd=0.05),1),0)
i_almost=imputation_accuracy(va,probbase_almost,prior,blocks)


##**********************************************************************
##  Test slight probbase pertubation                                ####
##**********************************************************************


is_real= imputation_accuracy(va,probbase_numeric,prior,blocks)

nsim_perturb=100 # Number of simulations
is_test=rep(0,nsim_perturb)
for (i in 1:length(is_test)) {
  set.seed(seed + i)
  probbase_almost=probbase_numeric
  w2=sample(wf,20) # 20 elements perturbed
  probbase_almost[w2]=pmax(pmin(probbase_numeric[w2]+rnorm(w2,sd=0.05),1),0)
  is_test[i]=imputation_accuracy(va,probbase_almost,prior,blocks)
}


##**********************************************************************
##  See if perturbed elements can be found                          ####
##**********************************************************************


# Choose elements to perturb
if (perturb_all) {
  wp=1:length(probbase_numeric)
} else {
  w1=which(probbase_numeric >= ptb & probbase_numeric <= 1-ptb)
  pr=outer(rep(1,dim(probbase)[1]),prior); pr[,c(1:16,77:81)]=0
  w2=which(pr>= 1e-3)
  wp=intersect(w1,w2)
}
test_elements=sample(wp,nptb)
nontest_elements=setdiff(wp,test_elements)

# Shuffle indices (makes no difference but means effect can be seen sooner)
test_elements=sample(test_elements)
nontest_elements=sample(nontest_elements)

# Perturb these elements by +/- ptb (about 0.1)
probbase_test=probbase_numeric
orig=probbase_test[test_elements]
if (uniform_pertubation) {
  probbase_test[test_elements] = pmax(pmin(orig+sample(c(ptb,-ptb),length(test_elements),rep=T),1),0) 
} else {
  probbase_test[test_elements] = pmax(pmin(orig+rnorm(length(orig),sd=ptb),1),0) 
}

# Overall imputation accuracy with pertubation
i_test=imputation_accuracy(va,probbase_test,prior,blocks)

# Compare partial derivatives
# Let P[i,j,x] be the perturbed probbase with the i,j th element changed by adding x.
# Let I(P) be the imputation accuracy of a probbase P. 
# The matrix action_mat_plus is defined as action_mat[i,j] = I(P[i,j,0.1]-P[i,j,0])
action_mat_plus=0*probbase_numeric 
action_mat_minus=0*probbase_numeric 
action_mat_true=probbase_test-probbase_numeric
print(range(action_mat_true)) # as a check

for (i in 1:length(test_elements)) {
  oval=probbase_test[test_elements[i]]
  pb_minus=probbase_test; pb_plus=probbase_test
  pb_plus[test_elements[i]]=pmax(oval+epsilon,0)
  pb_minus[test_elements[i]]=pmax(oval-epsilon,0)
  val_plus=imputation_accuracy(va,pb_plus,prior,blocks)
  val_minus=imputation_accuracy(va,pb_minus,prior,blocks)
  action_mat_plus[test_elements[i]]=(val_plus-i_test)/epsilon
  action_mat_minus[test_elements[i]]=(val_minus-i_test)/epsilon
  if ((i%%10)==0) {
    if (running_figures) {
      par(mfrow=c(1,2))
      image(action_mat_plus)
      image(action_mat_true)
    }
    save.image(file=save_file)
    print(i)
  }
}
for (i in 1:length(nontest_elements)) {
  oval=probbase_test[nontest_elements[i]]
  pb_minus=probbase_test; pb_plus=probbase_test
  pb_plus[nontest_elements[i]]=pmax(oval+epsilon,0)
  pb_minus[nontest_elements[i]]=pmax(oval-epsilon,0)
  val_plus=imputation_accuracy(va,pb_plus,prior,blocks)
  val_minus=imputation_accuracy(va,pb_minus,prior,blocks)
  action_mat_plus[nontest_elements[i]]=(val_plus-i_test)/epsilon
  action_mat_minus[nontest_elements[i]]=(val_minus-i_test)/epsilon
  if ((i%%10)==0) {
    if (running_figures) {
      par(mfrow=c(1,2))
      image(action_mat_plus)
      image(action_mat_true)
    }
    save.image(file=save_file)
    print(i)
  }
}

save.image(file=save_file)





##**********************************************************************
##  Test perturbed individual elements                              ####
##**********************************************************************


i_test=imputation_accuracy(va,probbase_numeric,prior,blocks)

n_ex=20 # Number of samples
ns=20 # Resolution
del=0.2 # zoom in this much

dat_example=list()
for (ii in 1:n_ex) {
  set.seed(seed + ii)
  
  s=sample(test_elements,1)
  sq=seq(0,1,length=ns)
  ival=rep(0,length(sq))
  for (i in 1:length(sq)) {
    pt2=probbase_numeric
    pt2[s]=sq[i]
    isub=imputation_accuracy(va,pt2,prior,blocks)
    ival[i]=isub-i_test
    print(i)
  }
  
  # Refine
  s2=sq[which.min(ival)]
  sq2=seq(pmax(s2-del,0),pmin(s2+del,1),length=ns)
  ival2=rep(0,length(sq))
  for (i in 1:length(sq2)) {
    pt2=probbase_numeric
    pt2[s]=sq2[i]
    isub=imputation_accuracy(va,pt2,prior,blocks)
    ival2[i]=isub-i_test
    print(i)
  }
  
  sq=c(sq,sq2)
  ival=c(ival,ival2)
  ox=order(sq)
  sq=sq[ox]; ival=ival[ox]
  
  dat_example[[ii]]=list(sample=s,test=sq,imp=ival)
  
}

save.image(file=save_file)

