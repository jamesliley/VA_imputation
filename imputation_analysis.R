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

# Save figures to PDF?
save_pdf=TRUE

##**********************************************************************
##  Load data as lists                                              ####
##**********************************************************************

lx=list.files("../Images/",pattern=paste0("imputation_test_temp_",nsim,"*"),full.names=TRUE)
n_done=length(lx)

all_data=list()
for (i in 1:n_done) {
  load(lx[i],  temp_env <- new.env())
  all_data[[i]] = as.list(temp_env)
}


##**********************************************************************
##  Load some common data                                           ####
##**********************************************************************

# Imputation accuracy with no perturbation
is_real=all_data[[1]]$is_real
probbase=all_data[[1]]$probbase
probbase_numeric=all_data[[1]]$probbase_numeric


##**********************************************************************
##  Evaluate imputation accuracy for very perturbed probbases       ####
##**********************************************************************


i_rand=unlist(lapply(all_data,function(x) x$i_rand))   # With random probbase
i_scrambled=unlist(lapply(all_data,function(x) x$i_scrambled))   # With scrambled probbase
i_perturbed=unlist(lapply(all_data,function(x) x$i_perturbed))    # With globally perturbed probbase
i_almost=unlist(lapply(all_data,function(x) x$i_perturbed))    # With slightly perturbed probbase

# Test runs of slight perturbations
is_test=unlist(lapply(all_data,function(x) x$is_test))


# Print results
cat("Perturbation of twenty elements by N(0,0.05^2) led to higher score in this proportion of cases: \n")
print(paste0(length(which(is_test>is_real))," / ",length(is_test)))

xtab=c(nsim,
       signif(c(mean(i_rand),
       mean(i_scrambled),
       mean(i_perturbed),
       mean(is_test),
       quantile(is_test,c(0.25,0.75)),
       is_real)/(nsim*dim(probbase)[1]),digits=4))

names(xtab)=c("N_VA","Random","Scrambled","All_Perturbed",
              "Slight_perturbed","Slight_25q","Slight_75q","Correct")

print(xtab)


##**********************************************************************
##  Display results                                                 ####
##**********************************************************************

action_mat_plus=unlist(lapply(all_data,function(x) x$action_mat_plus))
action_mat_minus=unlist(lapply(all_data,function(x) x$action_mat_minus))
action_mat_true=unlist(lapply(all_data,function(x) x$action_mat_true))

test_elements=c(); nontest_elements=c()
for (i in 1:length(all_data)) {
  test_elements=c(test_elements,(i-1)*246*81 + all_data[[i]]$test_elements)
  nontest_elements=c(nontest_elements,(i-1)*246*81 + all_data[[i]]$nontest_elements)
}


# Predictor
pred=(action_mat_plus-action_mat_minus)

# Output
y=(action_mat_true)

## *ROC curves ####
# ROC curve: differentiation of too high/too low
if (save_pdf) pdf(paste0("Figures/roc_",nsim,"_posneg.pdf"),width=4,height=4)
yd=y[test_elements]>0; pd=pred[test_elements]
rocd=suppressWarnings(getroc(yd,pd))
plot(rocd,main=paste0("ROC for +/- perturb., n. VAs: ",nsim),addauc=TRUE)
if (save_pdf) dev.off()
cat(paste0("Positive vs Negative","\n"))
wpn=wilcox.test(pd[which(yd==1)],pd[which(yd==0)])
print(wpn)
print(paste0("P-val: ",wpn$p.value))

# ROC curve: differentiation of too high/unperturbed
if (save_pdf) pdf(paste0("Figures/roc_",nsim,"_posnone.pdf"),width=4,height=4)
yd2=y[c(test_elements,nontest_elements)]>0; pd2=pred[c(test_elements,nontest_elements)]
rocd2=suppressWarnings(getroc(yd2,pd2))
plot(rocd2,main=paste0("ROC for +/0 perturb., n. VAs: ",nsim),addauc=TRUE)
if (save_pdf) dev.off()
cat(paste0("Positive vs unperturbed","\n"))
wpu=wilcox.test(pd2[which(yd2==1)],pd2[which(yd2==0)])
print(wpu)
print(paste0("P-val: ",wpu$p.value))

# ROC curve: differentiation of too low/unperturbed
if (save_pdf) pdf(paste0("Figures/roc_",nsim,"_negnone.pdf"),width=4,height=4)
yd3=y[c(test_elements,nontest_elements)]<0; pd3=-pred[c(test_elements,nontest_elements)]
rocd3=suppressWarnings(getroc(yd3,pd3))
plot(rocd3,main=paste0("ROC for -/0 perturb., n. VAs: ",nsim),addauc=TRUE)
if (save_pdf) dev.off()
cat(paste0("Negative vs unperturbed","\n"))
wnu=wilcox.test(pd3[which(yd3==1)],pd3[which(yd3==0)])
print(wnu)
print(paste0("P-val: ",wnu$p.value))






## *Densities of values ####
if (save_pdf) pdf(paste0("Figures/density_",nsim,".pdf"),width=6,height=4)

# Data
sc=(2*nsim*dim(all_data[[1]]$probbase)[1]*all_data[[1]]$epsilon) # scale
pred=(action_mat_plus-action_mat_minus)/sc # Predictor
y=(action_mat_true) # Target
py=pred[c(test_elements,nontest_elements)]
yr=y[c(test_elements,nontest_elements)]

# Transform (arctan)
sc1=100
py=atan(py*sc1)

# Quantities to plot
y0=py[which(abs(yr)<0.05)]
yp=py[which(yr> 0.05)]
yn=py[which(yr< -0.05)]


# Densities and frequencies
b_adj=2
d0=density(y0,from=-pi/2,to=pi/2,adjust=b_adj)
dn=density(yn,from=-pi/2,to=pi/2,adjust=b_adj)
dp=density(yp,from=-pi/2,to=pi/2,adjust=b_adj)

# Set up plot
par(mar=c(4,4,2,4))

# Plot densities
plot(d0,xlim=(pi/2)*c(-1,1),
     xlab=expression(paste("f(",gamma[jl],")")),
     main=paste0("Density of predictors; n. VAs: ",nsim))
lines(dp,col="red")
lines(dn,col="blue")

# Ratios
dpr=dp$y/dn$y; dnr=1/dpr; 
d1=max(c(d0$y,dp$y,dn$y))
d2=quantile(c(dpr,dnr),0.9)
dlim=d1/d2
lines(dp$x,dpr*dlim,col="red",lty=2)
lines(dp$x,dnr*dlim,col="blue",lty=2)
pxx=pretty(seq(0,d2,length=5))
axis(4,at=seq(0,d1,length=length(pxx)),labels=pxx)
mtext("Density ratio",4,line=2)

# Add legend
legend("topright",c("Pos.","Neg.","None","Pos/Neg","Neg/Pos"),
       lty=c(1,1,1,2,2),col=c("red","blue","black","red","blue"),
       bg="white",bty="n")

# Close
if (save_pdf) dev.off()


##**********************************************************************
##  Prediction with extreme values                                  ####
##**********************************************************************

pred=(action_mat_plus-action_mat_minus)/sc # Predictor
y=(action_mat_true) # Target
py=pred[c(test_elements,nontest_elements)]
yr=y[c(test_elements,nontest_elements)]
yr[which(yr>1e-5)]=1; yr[which(yr < -1e-5)]=-1 # Binary
thresh=1/5

print(paste0("Perturbation of values (out of total: ",length(yr),") for which estimated partial derivative > ",thresh,": "))
print(table(yr[which(py>thresh)]))


print(paste0("Perturbation of values (out of total: ",length(yr),") for which estimated partial derivative < -",thresh,": "))
print(table(yr[which(py< -thresh)]))

print(paste0("Proportion of values for which absolute estimated partial derivative < -",thresh,": "))
print(paste0(length(which(abs(py)>thresh))," / ",length(py)))


##**********************************************************************
##  Test perturbed individual elements                              ####
##**********************************************************************


dat_example=list()
for (i in 1:length(all_data)) dat_example=c(dat_example,all_data[[i]]$dat_example)
n_ex=length(dat_example)

for (ii in 1:n_ex) {
  
  if (save_pdf) pdf(paste0("Figures/Examples/grad_",nsim,"_",ii,".pdf"),width=4,height=4)
  plot(dat_example[[ii]]$test,dat_example[[ii]]$imp,
       type="l",xlab="Probbase element",
       ylab=expression(paste(Delta," imputation accuracy")))
  abline(v=probbase_numeric[dat_example[[ii]]$sample],col="red")
  legend("topright","True value",lty=1,col="red")
  
  if (save_pdf) dev.off()
  
}


