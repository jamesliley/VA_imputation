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
library(openVA)
library(nbc4va)

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
ptb=0.1 # Perturb probbase entries by this much

# Perturb this many elements
nptb=200

# Random seed
seed=3726152
set.seed(seed)

##**********************************************************************
##  Learn parameters for simulation                                 ####
##**********************************************************************

code=RandomPhysician$code1 # Physician coded causes of death
va=RandomPhysician[,1:246] # VA questionnaire data

# Learn covariance structure from data
rparam=learn_parameters_from_va(va,code)

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

nsim_perturb=100
i_test=rep(0,nsim_perturb)
for (i in 1:length(i_test)) {
  probbase_almost=probbase_numeric
  w2=sample(wf,20)
  probbase_almost[w2]=pmax(pmin(probbase_numeric[w2]+rnorm(w2,sd=0.05),1),0)
  i_test[i]=imputation_accuracy(va,probbase_almost,prior,blocks)
}
print(length(which(i_test>i_real)))


##**********************************************************************
##  See if perturbed elements can be found                          ####
##**********************************************************************

# Choose elements to perturb
wp=which(probbase_numeric>ptb & probbase_numeric<1-ptb)
test_elements=sample(wp,nptb)
nontest_elements=setdiff(wp,test_elements)

# Shuffle indices (makes no difference but means effect can be seen sooner)
test_elements=sample(test_elements)
nontest_elements=sample(nontest_elements)

# Perturb these elements by +/- ptb (about 0.1)
probbase_test=probbase_numeric
orig=probbase_test[test_elements]
probbase_test[test_elements] = orig+sample(c(ptb,-ptb),length(test_elements),rep=T)

# Overall imputation accuracy with pertubation
i_test=imputation_accuracy(va,probbase_test,prior,blocks)

# Compare partial derivatives
# Let P[i,j,x] be the perturbed probbase with the i,j th element changed by adding x.
# Let I(P) be the imputation accuracy of a probbase P. 
# The matrix action_mat_plus is defined as action_mat[i,j] = I(P[i,j,0.1]-P[i,j,0])
action_mat_plus=0*probbase_numeric 
action_mat_minus=0*probbase_numeric 
action_mat_true=probbase_test-probbase_numeric
for (i in 1:length(test_elements)) {
  oval=probbase_test[test_elements[i]]
  pb_minus=probbase_test; pb_plus=probbase_test
  pb_plus[test_elements[i]]=pmax(oval+ptb,0)
  pb_minus[test_elements[i]]=pmax(oval-ptb,0)
  val_plus=imputation_accuracy(va,pb_plus,prior,blocks)
  val_minus=imputation_accuracy(va,pb_minus,prior,blocks)
  action_mat_plus[test_elements[i]]=val_plus-i_test
  action_mat_minus[test_elements[i]]=val_minus-i_test
  if ((i%%10)==0) {
    par(mfrow=c(1,2))
    image(action_mat_plus)
    image(action_mat_true)
    save.image(file="imputation_test_temp.RData")
    print(i)
  }
}
for (i in 1:length(nontest_elements)) {
  oval=probbase_test[nontest_elements[i]]
  pb_minus=probbase_test; pb_plus=probbase_test
  pb_plus[nontest_elements[i]]=pmax(oval+ptb,0)
  pb_minus[nontest_elements[i]]=pmax(oval-ptb,0)
  val_plus=imputation_accuracy(va,pb_plus,prior,blocks)
  val_minus=imputation_accuracy(va,pb_minus,prior,blocks)
  action_mat_plus[nontest_elements[i]]=val_plus-i_test
  action_mat_minus[nontest_elements[i]]=val_minus-i_test
  if ((i%%10)==0) {
    par(mfrow=c(1,2))
    image(action_mat_plus)
    image(action_mat_true)
    save.image(file="imputation_test_temp.RData")
    print(i)
  }
}

##**********************************************************************
##  Display results                                                 ####
##**********************************************************************

# Predictor
pred=(action_mat_plus-action_mat_minus)

# Output
y=(action_mat_true)

# ROC curve: differentiation of too high/too low
pdf(paste0("Figures/roc_",nsim,".pdf"),width=4,height=4)
yd=y[test_elements]>0; pd=pred[test_elements]
rocd=getroc(yd,pd)
plot(rocd,main="ROC for +/- pertubation",addauc=TRUE)
print(wilcox.test(pd[which(yd==1)],pd[which(yd==0)]))
dev.off()

# Densities of values
pdf(paste0("Figures/density_",nsim,".pdf"),width=6,height=4)
if (nsim>5000) {
  xr=10;yr=0.03; dr=c(-100,100)
} else {
  xr=2; yr=0.17; dr=c(-12,12)
}
par(mar=c(4,4,2,4))
d0=density(pred[nontest_elements],from=dr[1],to=dr[2])
dp=density(pd[which(yd>0)],from=dr[1],to=dr[2])
dn=density(pd[which(yd==0)],from=dr[1],to=dr[2])
plot(d0,xlim=dr,main="",xlab=expression(gamma[jl]))
lines(dp,col="red")
lines(dn,col="blue")
dpr=dp$y/dn$y; dnr=1/dpr
lines(dp$x,dpr*max(d0$y)/max(dpr),col="red",lty=2)
lines(dp$x,dnr*max(d0$y)/max(dpr),col="blue",lty=2)
pxx=pretty(seq(0,max(dpr),length=5))
axis(4,at=seq(min(d0$y),max(d0$y),length=length(pxx)),labels=pxx)
mtext("Density ratio",4,line=2)

legend(xr,yr,c("Pos.","Neg.","None","Pos/Neg","Neg/Pos"),
       lty=c(1,1,1,2,2),col=c("red","blue","black","red","blue"),
       bg="white",bty="n")
dev.off()

# ROC curve: identification of too high
roch=getroc(yp[wp]>0,pred[wp])
plot(rocd,main="ID: high",addauc=TRUE)

# ROC curve: identification of too low
roch=getroc(yp[wp]>0,pred[wp])
plot(rocd,main="ID: high",addauc=TRUE)



# Heatmaps: not very easy to read
par(mfrow=c(1,2))
disp_mat=NA*probbase_numeric
disp_mat[wp]=0
disp_mat[test_elements]=probbase_numeric[test_elements]-probbase_test[test_elements]
crange=colorRampPalette(c("blue","gray","red"))(100)
image(disp_mat[,17:76],xlab="Questions",ylab="CoDs",col=crange,main="Perturbed values")
#legend("bottomright",c("Too high","Too low","Candidate"),col=c("red","blue","gray"),pch=16)

dval=action_mat_plus-action_mat_minus; nd=length(dval)
high_thresh=quantile(dval,1-nptb/(2*nd))
low_thresh=quantile(dval,nptb/(2*nd))
dec_mat=NA*probbase_numeric
dec_mat[wp]=0
dec_mat[which(dval>high_thresh)]=-1
dec_mat[which(dval<low_thresh)]=1

image(dec_mat[,17:76],xlab="Questions",ylab="CoDs",col=crange,main="Identified values")
#legend("bottomright",c("Too high","Too low","Candidate"),col=c("red","blue","gray"),pch=16)

