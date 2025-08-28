##**********************************************************************
##  Calibration of symptom frequencies in VA using imputation       ####
##   proposal. 
##   
##  Working folder for VA_imputations
##  
##  Nathan Higgins, James Liley, Eilidh Cowan                                                        
##  July 2025                                                              
##**********************************************************************
##
## 

## Working area before functions added to functions file

library("openVA")
library("nbc4va")


##**********************************************************************
#General R package creation code
##**********************************************************************
##
#install.packages("devtools")
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

#Guide: https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/comment-page-1/#comments

#Creating an r package
setwd(path.expand("~"))  # Sets working directory to your home folder
create("calibrateVA")

#In calibrateVA, go in R folder, create new r file and add new function and documentation
#as in website(or create new folder each time)

#Process your documentation
setwd("./calibrateVA")
document()

#Install
setwd("..")
install("calibrateVA")






##**********************************************************************
#First attempt at simulating VA data- 10 arbitary question and 3 arbitary causes
##**********************************************************************
##


#Simulation of VA data

## Simulate VA example
library(mnormt)

# n: number of simulations
# pi: disease frequency
# pb: probbase
#simulate_va=function(n,pi,pb,blocks,covs)

n=5
#Take 3 diseases- we choose pi
pi=c(0.3,0.5,0.2)

# Take probbase for 3 causes and 10 questions- we choose
pb=matrix(runif(length(pi)*10),10,length(pi))

#List blocks of questions which are conditionally independent
blocks=list(c(1:3),c(4:6),c(7:10))


dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc
#Covariance matrix for each block of questions
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
va









##**********************************************************************
#First attempt at creating table of blocks of questions
##**********************************************************************
##

#Defintely needs input from an expert in best way to group questions-I have done best I can with my limited knowledge in this area

#This is from form of randomVA1 questions
data("probbase")
question_text<- probbase[-1,2:3]

data("RandomVA1")
names<-list(colnames(RandomVA1)[1:10], colnames(RandomVA1)[11:13],colnames(RandomVA1)[14:17],colnames(RandomVA1)[18:20],colnames(RandomVA1)[21:22],
            colnames(RandomVA1)[23],colnames(RandomVA1)[24],colnames(RandomVA1)[25],colnames(RandomVA1)[26],colnames(RandomVA1)[27],colnames(RandomVA1)[28],
            colnames(RandomVA1)[29],colnames(RandomVA1)[30],colnames(RandomVA1)[31],colnames(RandomVA1)[32],colnames(RandomVA1)[33],
            colnames(RandomVA1)[34],colnames(RandomVA1)[35],colnames(RandomVA1)[36],colnames(RandomVA1)[37],colnames(RandomVA1)[38],colnames(RandomVA1)[39:40],colnames(RandomVA1)[41:42],colnames(RandomVA1)[43:46],
            colnames(RandomVA1)[47:52],colnames(RandomVA1)[53:56],colnames(RandomVA1)[57:62],
            colnames(RandomVA1)[63:64], colnames(RandomVA1)[65],colnames(RandomVA1)[66:70],colnames(RandomVA1)[71:72],colnames(RandomVA1)[73:76],
            colnames(RandomVA1)[77:79], colnames(RandomVA1)[80:82],colnames(RandomVA1)[83],colnames(RandomVA1)[84:86],colnames(RandomVA1)[87:91],colnames(RandomVA1)[92:94],
            colnames(RandomVA1)[95:96],colnames(RandomVA1)[97:100],colnames(RandomVA1)[101:104],colnames(RandomVA1)[105:106],colnames(RandomVA1)[107:108],colnames(RandomVA1)[109:116],
            colnames(RandomVA1)[117],colnames(RandomVA1)[118],colnames(RandomVA1)[119],colnames(RandomVA1)[120:121],colnames(RandomVA1)[122:125],colnames(RandomVA1)[126],colnames(RandomVA1)[127:128],
            colnames(RandomVA1)[129:130],colnames(RandomVA1)[131:134],colnames(RandomVA1)[135:138],colnames(RandomVA1)[139], colnames(RandomVA1)[140:146],colnames(RandomVA1)[147], colnames(RandomVA1)[148:149],
            colnames(RandomVA1)[150:152],colnames(RandomVA1)[153:157],colnames(RandomVA1)[158:159],colnames(RandomVA1)[160:161],colnames(RandomVA1)[162:164],colnames(RandomVA1)[165:167],colnames(RandomVA1)[168:170],
            colnames(RandomVA1)[171:178],colnames(RandomVA1)[179:180],colnames(RandomVA1)[181:184],colnames(RandomVA1)[185:186],colnames(RandomVA1)[187:188],colnames(RandomVA1)[189:194],colnames(RandomVA1)[195:196],
            colnames(RandomVA1)[197],colnames(RandomVA1)[198:200],colnames(RandomVA1)[201:202],colnames(RandomVA1)[203:207],colnames(RandomVA1)[208:211],colnames(RandomVA1)[c(212, 215:217)],colnames(RandomVA1)[213:214],
            colnames(RandomVA1)[c(218, 222)],colnames(RandomVA1)[219:221],colnames(RandomVA1)[223],colnames(RandomVA1)[224],colnames(RandomVA1)[225],colnames(RandomVA1)[226],colnames(RandomVA1)[227],colnames(RandomVA1)[228:233],
            colnames(RandomVA1)[234:235],colnames(RandomVA1)[236],colnames(RandomVA1)[237:246])

# Your existing list of grouped column names
grouped_names <- names

# Compute lengths of each group
group_lengths <- lengths(grouped_names)

# Create list of indices corresponding to each group
index_list <- mapply(
  function(start, len) seq(from = start, length.out = len),
  start = cumsum(c(1, head(group_lengths, -1))),
  len = group_lengths,
  SIMPLIFY = FALSE
)

Index<-lapply(index_list, as.numeric)


block_text<- data.frame(Names = I(names),
                        Index = I(Index))



#Testing to see if works

#How to input blocks from list of blocks
str(block_text$Index[[2]])
data(RandomVA1)
out=recoverAnswers(RandomVA1,block=block_text$Index[[2]],method="OpenVA",
                   data.type = "WHO2012",model = "InterVA",
                   version = "4.03", HIV = "h", Malaria = "h")

head(out$Original)
head(out$Imputed)
names[[3]]



## Impute blocks of questions 18-20 and 21-24, use probbase condprobnum and change element (1,1) of probbase by 0.1
data("condprobnum")
data("RandomVA1")
out<- gradient_probbase(RandomVA1[1:100,],block=block_text$Index[2:3],method="OpenVA",data.type = "WHO2012",
                        model = "InSilicoVA",CondProbNum=condprobnum, probbase_element = list(1,1) , epsilon = 0.1)





##**********************************************************************
#Simulating VA without first 9 questions(age and sex)
##**********************************************************************
##



#Eventually use probbase of form condprobnum- 60 causes and 245 questions
#Simulate n VAs, each of which answer 245 questions

#First try with real probbase and real subset of questions
#Use condprobnum and subset of blocks from other file

data(condprobnum)

## Simulate VA example
library(mnormt)

# n: number of simulations
# pi: disease frequency
# pb: probbase
#simulate_va=function(n,pi,pb,blocks,covs)

n=5
#Take 3 diseases- we choose pi
#Fixing pi
#pi=c(0.2,0.5,0.2,0.1)

#Randomly assign pi
random_vector <- runif(60)
# Normalize to sum to 1
pi <- random_vector / sum(random_vector)


# Take probbase for 3 causes and first 10 'proper' questions(ignore first 9 for now as they are different type)- we choose
#pb= condprobnum[10:245,]
#Or random probbase
pb=matrix(runif(length(pi)*236),nrow=236,ncol=length(pi))

#List blocks of questions which are conditionally independent
blocks=block_text$Index

dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc
#Covariance matrix for each block of questions
covs<-vector("list", length = length(blocks))
for (k in 1:length(blocks)) {
  covs[[k]]<- dcov(length(blocks[[k]]),0.1)
}



# Number of causes of death
n_cod=as.numeric(length(pi))

# Number of questions
n_q=as.numeric(nrow(pb))

# Number of blocks
n_block=as.numeric(length(blocks))


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
    
    bprobs=qprobs[wq-10]
    thresholds=qnorm(1-bprobs)
    
    answers[wq-10]=latent_q>thresholds
  }
  va[i,]=answers
}
va

colnames(va)<- unlist(block_text$Names[-1])
va <- ifelse(va, "Y", "")

#Now add first 9 question simulation

#Need to simulate id column, age, and sex



#Need to generalise for all causes, then assign col/row names, then add in first 9 questions in randomVA1 form

#Then need to fix question which men/women only, then find correct covariance structure for binary question etc

#Upload all code on simulation/blocks/r packages to working folder for james to check- ask help on blocks/choices in simulation








##**********************************************************************
#Simulating VA including age and sex- via probbase
##**********************************************************************
##


#Version that simulates all questions- can give silly answers for first 9-need to simulate seperately, ie cant have not male or female
data(condprobnum)

## Simulate VA example
library(mnormt)

# n: number of simulations
# pi: disease frequency
# pb: probbase
#simulate_va=function(n,pi,pb,blocks,covs)

n=5
#Take 3 diseases- we choose pi
#Fixing pi
#pi=c(0.2,0.5,0.2,0.1)

#Randomly assign pi
random_vector <- runif(60)
# Normalize to sum to 1
pi <- random_vector / sum(random_vector)


# Take probbase for 3 causes and first 10 'proper' questions(ignore first 9 for now as they are different type)- we choose
#pb= condprobnum[1:245,]
#Or random probbase
pb=matrix(runif(length(pi)*245),nrow=245,ncol=length(pi))

#List blocks of questions which are conditionally independent
blocks=block_text$Index

dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc
#Covariance matrix for each block of questions
covs<-vector("list", length = length(blocks))
for (k in 1:length(blocks)) {
  covs[[k]]<- dcov(length(blocks[[k]]),0.1)
}



# Number of causes of death
n_cod=as.numeric(length(pi))

# Number of questions
n_q=as.numeric(nrow(pb))

# Number of blocks
n_block=as.numeric(length(blocks))


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
    
    bprobs=qprobs[wq-10]
    thresholds=qnorm(1-bprobs)
    
    answers[wq-10]=latent_q>thresholds
  }
  va[i,]=answers
}
va

colnames(va)<- unlist(block_text$Names)[-1]
va <- ifelse(va, "Y", "")

#Need to generalise for all causes, then assign col/row names, then add in first 9 questions in randomVA1 form

#Then need to fix question which men/women only, then find correct covariance structure for binary question etc

#Upload all code on simulation/blocks/r packages to working folder for james to check- ask help on blocks/choices in simulation



##**********************************************************************
#Simulating VA including age and sex- via seperate simulation
##**********************************************************************
##

#Eventually use probbase of form condprobnum- 60 causes and 245 questions
#Simulate n VAs, each of which answer 245 questions

#First try with real probbase and real subset of questions
#Use condprobnum and subset of blocks from other file

data(condprobnum)

## Simulate VA example
library(mnormt)

# n: number of simulations
# pi: disease frequency
# pb: probbase
#simulate_va=function(n,pi,pb,blocks,covs)

n=20
#Take 3 diseases- we choose pi
#Fixing pi
#pi=c(0.2,0.5,0.2,0.1)

#Randomly assign pi
random_vector <- runif(60)
# Normalize to sum to 1
pi <- random_vector / sum(random_vector)


# Take probbase for 3 causes and first 10 'proper' questions(ignore first 9 for now as they are different type)- we choose
#pb= condprobnum[10:245,]
#Or random probbase
pb=matrix(runif(length(pi)*236),nrow=236,ncol=length(pi))

#List blocks of questions which are conditionally independent
blocks=block_text$Index

dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc
#Covariance matrix for each block of questions
covs<-vector("list", length = length(blocks))
for (k in 1:length(blocks)) {
  covs[[k]]<- dcov(length(blocks[[k]]),0.1)
}



# Number of causes of death
n_cod=as.numeric(length(pi))

# Number of questions
n_q=as.numeric(nrow(pb))

# Number of blocks
n_block=as.numeric(length(blocks))


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
    
    #-10 as we ignoring first 9 questions at the moment and simulate that seperately
    bprobs=qprobs[wq-10]
    thresholds=qnorm(1-bprobs)
    
    answers[wq-10]=latent_q>thresholds
  }
  va[i,]=answers
}
va

#Apply correct coloumn names
colnames(va)<- unlist(block_text$Names[-1])

#Change TRUE/FALSE to Y/""
va <- ifelse(va, "Y", "")

#Now add first 9 question simulation

#Need to simulate id column, age, and sex

#Simulate ID
ID<- rep(0,n)
for(i in 1:n){
  ID[i]<-paste0("d",i)
}

#Simulate age

elder<-rep("",n)
midage<- rep("",n)
adult<- rep("",n)
child<- rep("",n)
under5<- rep("",n)
infant<- rep("",n)
neonate<- rep("",n)
for (i in 1:n) {
  x<- sample(c(0,1,2,3,4,5,6),1,prob = c(1/3,1/3,1/3,0,0,0,0))
  if(x==0){
    elder[i]<-"Y"
  }
  if(x==1){
    midage[i]<-"Y"
  }
  if(x==2){
    adult[i]<-"Y"
  }
  if(x==3){
    child[i]<-"Y"
  }
  if(x==4){
    under5[i]<-"Y"
  }
  if(x==5){
    infant[i]<-"Y"
  }
  if(x==6){
    neonate[i]<-"Y"
  }
}

#Simulate sex
male<-rep("",n)
female<- rep("",n)
for (i in 1:n) {
  x<- sample(c(0,1),1,prob = c(0.5,0.5))
  if(x==0){
    male[i]<-"Y"
  }
  if(x==1){
    female[i]<-"Y"
  }
}

#Combining to get final VA
va<- cbind(ID,elder,midage,adult,child,under5,infant,neonate,male,female,va)






##**********************************************************************
#Simulating VA function-second attempt
##**********************************************************************
##

#Simulation of VA data

##' Simulation of VA data (WHO2012)
##' @name simulateVA
##' @description Simulates VA data. Format is WHO2012. Example; data(RandomVA1). Compatible for InterVA4 and InSilicoVA models
##' @param n Number of simulations
##' @param blocks A list of blocks of questions. Each block is a set of indices for questions, corresponding to rows in VA. MUST EXCLUDE the first block containing 9 questions for sex/age. See data(RandomVA1) for explicit coloumns.
##' @param probBase A customised probability matrix with 245 row symptomns on 60 column causes in the same order as InterVA4 specification. For example input see condprobnum.
##' @param pi Vector of prior probabilities for the 60 causes of death- must sum to 1
##' @param covs List of Covariance matrices for each block of questions in argument blocks
##' @param ageProb Vector of probability for age in order c(elder,midage,adult,child,under5,infant,neonate) 
##' @param sexProb Vector of probability for sex in order c(male,female) 
##' @param woman_block A subset of the list blocks for which the questions are for female deaths only
##' @param neonate_block A subset of the list blocks for which the questions are for neonates deaths only
##' @return n simulations of WHO2012 VA data
##' @export
##' @examples
##'#Number of simulations
##'n=100
##'
##'#Randomly assign pi
##'random_vector <- runif(60)
##'# Normalize to sum to 1
##'pi <- random_vector / sum(random_vector)
##'
##'
##'# Take known probbase (60 causes and first 245 questions)
##'data(condprobnum)
##'pb= condprobnum
##'
##'#Or random probbase
##'#pb=matrix(runif(length(pi)*245),nrow=245,ncol=length(pi))
##'#Different to previous attempts as we want input probBase to be full probbase(exactly same as condprobnum), not a subset of it(or letter probBase)
##'
##'#List blocks of questions which are conditionally independent-leaving out first block(age/sex questions)- needed for this version to work
##'blocks=block_text$Index[-1]
##'
##'
##'#Covariance function
##'dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc
##'
##'#Covariance matrix for each block of questions
##'covs<-vector("list", length = length(blocks))
##'for (k in 1:length(blocks)) {
##'covs[[k]]<- dcov(length(blocks[[k]]),0.1)
##'}
##'
##'
##'#Age probs
##'ageProb<- c(1/3,1/3,1/3,0,0,0,0)
##'
##'#Sex probs
##'sexProb<-c(0.5,0.5)
##'
##'#Blocks of questions just for women deaths
##'woman_block<-blocks[c(1,48,50:61)]
##'woman_block_index<-c(1,48,50:61)
##'woman_block_names<- block_text$Names[woman_block_index +1]
##'
##'Blocks of questions just for neonate deaths
##'neonate_block<-blocks[c(2,62:76)]
##'neonate_block_index<-c(2,62:76)
##'neonate_block_names<- block_text$Names[neonate_block_index +1]
##'
##'simVA<-simulateVA(n,pi,pb,blocks,covs,ageProb,sexProb,woman_block,neonate_block)
##'




simulateVA<- function(n,pi,probBase,blocks,covs, ageProb, sexProb, woman_block,neonate_block){
  
  ## package for multivariate normal simulation
  library(mnormt)
  
  # n: number of simulations
  # pi: disease frequency
  # pb: probbase
  #covs: Covariance matrix for each block of questions
  
  #Probbase must be numeric- need to extract probbase without Q1-9 for this method
  #Hence we say argument probbase specify full probbase, ie condprobnum then pb is just this without first 9 questions
  
  pb<- probBase[10:245,]
  
  # Number of causes of death
  n_cod=as.numeric(length(pi))
  
  # Number of questions
  n_q=as.numeric(nrow(pb))
  
  # Number of blocks(input excludes initial block of age/sex questions-simulated separately)
  n_block=as.numeric(length(blocks))
  
  #Simulate all questions apart from first 9
  
  # Firstly simulate under conditional independence
  va=matrix(NA,n,n_q)
  for (i in 1:n) {
    cause_of_death=sample(1:n_cod,1,prob=pi)
    
    qprobs=pb[,cause_of_death]
    
    answers=rep(NA,n_q)
    #Adjusting for the fact treating first block(9 questions) seperately
    for (b in 1:n_block) {
      wq=blocks[[b]] # indices of questions in block b
      bcov=covs[[b]] # covariance matrix for latent variable
      nb=length(wq) # number of questions in block
      
      latent_q=rmnorm(1,mean=rep(0,nb),varcov=bcov)
      
      #-10 as we ignoring first 9 questions at the moment and simulate that seperately
      bprobs=qprobs[wq-10]
      thresholds=qnorm(1-bprobs)
      
      answers[wq-10]=latent_q>thresholds
    }
    va[i,]=answers
  }
  
  #Apply correct coloumn names
  data("RandomVA1", package = "InSilicoVA")
  colnames(va)<- colnames(RandomVA1)[-c(1:10)]
  
  #Change TRUE/FALSE to Y/""
  va <- ifelse(va, "Y", "")
  
  #Now add first 9 question simulation
  
  #Need to simulate id column, age, and sex
  
  #Simulate ID
  ID<- rep(0,n)
  for(i in 1:n){
    ID[i]<-paste0("d",i)
  }
  
  #Simulate age
  
  elder<-rep("",n)
  midage<- rep("",n)
  adult<- rep("",n)
  child<- rep("",n)
  under5<- rep("",n)
  infant<- rep("",n)
  neonate<- rep("",n)
  for (i in 1:n) {
    x<- sample(c(0,1,2,3,4,5,6),1,prob = ageProb)
    if(x==0){
      elder[i]<-"Y"
    }
    if(x==1){
      midage[i]<-"Y"
    }
    if(x==2){
      adult[i]<-"Y"
    }
    if(x==3){
      child[i]<-"Y"
    }
    if(x==4){
      under5[i]<-"Y"
    }
    if(x==5){
      infant[i]<-"Y"
    }
    if(x==6){
      neonate[i]<-"Y"
    }
  }
  
  #Simulate sex
  male<-rep("",n)
  female<- rep("",n)
  for (i in 1:n) {
    x<- sample(c(0,1),1,prob = sexProb)
    if(x==0){
      male[i]<-"Y"
    }
    if(x==1){
      female[i]<-"Y"
    }
  }
  
  #Combining to get final VA
  va<- cbind(ID,elder,midage,adult,child,under5,infant,neonate,male,female,va)
  va<-as.data.frame(va)
  #Removing answers which don't make sense- due to age/sex
  for(i in 1:n){
    if(va$female[i] == ""){
      va[i,unlist(woman_block)]<- ""
    }
    if(va$neonate[i] == ""){
      va[i,unlist(neonate_block)]<- ""
    }
  }
  return(va)
}





##**********************************************************************
#Methods for generating covariance matrices for blocks
##**********************************************************************
##

#Covariance for blocks

#Method 1: all blocks have same covariance structure
dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc
#Covariance matrix for each block of questions
covs<-vector("list", length = length(blocks))

#Method 2: Introduce negative covariance for "binary" blocks
#Negative correlation for binary blocks of questions: dcov(2,-0.95)
#Need to loop around certain blocks

dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc

#Blocks of binary questions and names for reference
binary_blocks<-c(1,2,3,4,22,29,32,33,37,38,50,51,62,71,73,78)
binary_blocks_names<- block_text$Names[binary_blocks +1]

#Covariance matrix for each block of questions
covs<-vector("list", length = length(blocks))
for (k in 1:length(blocks)) {
  if(k %in% binary_blocks){
  covs[[k]]<- dcov(length(blocks[[k]]),-0.95)
  }
  if(!(k %in% binary_blocks)){
    covs[[k]]<- dcov(length(blocks[[k]]),0.1)
  }
}


#Method 3: Random covariance matrices for each block
#Need number of simulations n to be bigger than d, dimension of cov matrix(ie length of each block)
#library(mnormt)
#X=rmnorm(20,mean=rep(0,10),varcov=diag(10))
#M=cor(X)
#Not singular as non zero eigenvalue
#min(eigen(M)$values)
#Need to loop this around blocks


dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc

#Blocks of binary questions and names for reference
binary_blocks<-c(1,2,3,4,22,29,32,33,37,38,50,51,62,71,73,78)
binary_blocks_names<- block_text$Names[binary_blocks +1]

# Constants
max_len <- max(sapply(blocks, length)) + 20
covs <- vector("list", length = length(blocks))

library(mnormt)
for (k in seq_along(blocks)) {
  n_vars <- length(blocks[[k]])
  
  # Simulate data
  X <- rmnorm(max_len, mean = rep(0, n_vars), varcov = diag(n_vars))
  X <- matrix(X, ncol = n_vars)
  # Choose sign for correlations
  if (k %in% binary_blocks) {
    M <- -abs(cor(X))  # Negative correlations
  } else {
    M <- abs(cor(X))   # Positive correlations
  }
  
  diag(M) <- 1  # Ensure correlation matrix diagonal = 1
  
  covs[[k]] <- M
}





##**********************************************************************
#Attempting to use simulated data on gradient function-Issues
##**********************************************************************
##

#Simulation of VA data

##' Simulation of VA data (WHO2012)
##' @name simulateVA
##' @description Simulates VA data. Format is WHO2012. Example; data(RandomVA1). Compatible for InterVA4 and InSilicoVA models
##' @param n Number of simulations
##' @param blocks A list of blocks of questions. Each block is a set of indices for questions, corresponding to rows in VA. MUST EXCLUDE the first block containing 9 questions for sex/age. See data(RandomVA1) for explicit coloumns.
##' @param probBase A customised probability matrix with 245 row symptomns on 60 column causes in the same order as InterVA4 specification. For example input see condprobnum.
##' @param pi Vector of prior probabilities for the 60 causes of death- must sum to 1
##' @param covs List of Covariance matrices for each block of questions in argument blocks
##' @param ageProb Vector of probability for age in order c(elder,midage,adult,child,under5,infant,neonate) 
##' @param sexProb Vector of probability for sex in order c(male,female) 
##' @param woman_block A subset of the list blocks for which the questions are for female deaths only
##' @param neonate_block A subset of the list blocks for which the questions are for neonates deaths only
##' @return n simulations of WHO2012 VA data
##' @export
##' @examples
##'#Number of simulations
##'n=100
##'
##'#Randomly assign pi
##'random_vector <- runif(60)
##'# Normalize to sum to 1
##'pi <- random_vector / sum(random_vector)
##'
##'
##'# Take known probbase (60 causes and first 245 questions)
##'data(condprobnum)
##'pb= condprobnum
##'
##'#Or random probbase
##'#pb=matrix(runif(length(pi)*245),nrow=245,ncol=length(pi))
##'#Different to previous attempts as we want input probBase to be full probbase(exactly same as condprobnum), not a subset of it(or letter probBase)
##'
##'#List blocks of questions which are conditionally independent-leaving out first block(age/sex questions)- needed for this version to work
##'blocks=block_text$Index[-1]
##'
##'
##'#Covariance function
##'dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc
##'
##'#Covariance matrix for each block of questions
##'covs<-vector("list", length = length(blocks))
##'for (k in 1:length(blocks)) {
##'covs[[k]]<- dcov(length(blocks[[k]]),0.1)
##'}
##'
##'
##'#Age probs
##'ageProb<- c(1/3,1/3,1/3,0,0,0,0)
##'
##'#Sex probs
##'sexProb<-c(0.5,0.5)
##'
##'#Blocks of questions just for women deaths
##'woman_block<-blocks[c(1,48,50:61)]
##'woman_block_index<-c(1,48,50:61)
##'woman_block_names<- block_text$Names[woman_block_index +1]
##'
##'Blocks of questions just for neonate deaths
##'neonate_block<-blocks[c(2,62:76)]
##'neonate_block_index<-c(2,62:76)
##'neonate_block_names<- block_text$Names[neonate_block_index +1]
##'
##'simVA<-simulateVA(n,pi,pb,blocks,covs,ageProb,sexProb,woman_block,neonate_block)
##'




simulateVA<- function(n,pi,probBase,blocks,covs, ageProb, sexProb, woman_block,neonate_block){
  
  ## package for multivariate normal simulation
  library(mnormt)
  
  # n: number of simulations
  # pi: disease frequency
  # pb: probbase
  #covs: Covariance matrix for each block of questions
  
  #Probbase must be numeric- need to extract probbase without Q1-9 for this method
  #Hence we say argument probbase specify full probbase, ie condprobnum then pb is just this without first 9 questions
  
  pb<- probBase[10:245,]
  
  # Number of causes of death
  n_cod=as.numeric(length(pi))
  
  # Number of questions
  n_q=as.numeric(nrow(pb))
  
  # Number of blocks(input excludes initial block of age/sex questions-simulated separately)
  n_block=as.numeric(length(blocks))
  
  #Simulate all questions apart from first 9
  
  # Firstly simulate under conditional independence
  va=matrix(NA,n,n_q)
  for (i in 1:n) {
    cause_of_death=sample(1:n_cod,1,prob=pi)
    
    qprobs=pb[,cause_of_death]
    
    answers=rep(NA,n_q)
    #Adjusting for the fact treating first block(9 questions) seperately
    for (b in 1:n_block) {
      wq=blocks[[b]] # indices of questions in block b
      bcov=covs[[b]] # covariance matrix for latent variable
      nb=length(wq) # number of questions in block
      
      latent_q=rmnorm(1,mean=rep(0,nb),varcov=bcov)
      
      #-10 as we ignoring first 9 questions at the moment and simulate that seperately
      bprobs=qprobs[wq-10]
      thresholds=qnorm(1-bprobs)
      
      answers[wq-10]=latent_q>thresholds
    }
    va[i,]=answers
  }
  
  #Apply correct coloumn names
  data("RandomVA1", package = "InSilicoVA")
  colnames(va)<- colnames(RandomVA1)[-c(1:10)]
  
  #Change TRUE/FALSE to Y/""
  va <- ifelse(va, "Y", "")
  
  #Now add first 9 question simulation
  
  #Need to simulate id column, age, and sex
  
  #Simulate ID
  ID<- rep(0,n)
  for(i in 1:n){
    ID[i]<-paste0("d",i)
  }
  
  #Simulate age
  
  elder<-rep("",n)
  midage<- rep("",n)
  adult<- rep("",n)
  child<- rep("",n)
  under5<- rep("",n)
  infant<- rep("",n)
  neonate<- rep("",n)
  for (i in 1:n) {
    x<- sample(c(0,1,2,3,4,5,6),1,prob = ageProb)
    if(x==0){
      elder[i]<-"Y"
    }
    if(x==1){
      midage[i]<-"Y"
    }
    if(x==2){
      adult[i]<-"Y"
    }
    if(x==3){
      child[i]<-"Y"
    }
    if(x==4){
      under5[i]<-"Y"
    }
    if(x==5){
      infant[i]<-"Y"
    }
    if(x==6){
      neonate[i]<-"Y"
    }
  }
  
  #Simulate sex
  male<-rep("",n)
  female<- rep("",n)
  for (i in 1:n) {
    x<- sample(c(0,1),1,prob = sexProb)
    if(x==0){
      male[i]<-"Y"
    }
    if(x==1){
      female[i]<-"Y"
    }
  }
  
  #Combining to get final VA
  va<- cbind(ID,elder,midage,adult,child,under5,infant,neonate,male,female,va)
  va<-as.data.frame(va)
  #Removing answers which don't make sense- due to age/sex
  for(i in 1:n){
    if(va$female[i] == ""){
      va[i,unlist(woman_block)]<- ""
    }
    if(va$neonate[i] == ""){
      va[i,unlist(neonate_block)]<- ""
    }
  }
  return(va)
}
 



 #Inputs
  
  #Creating table of blocks of questions
  #This is from form of randomVA1 questions
  data("probbase")
  question_text<- probbase[-1,2:3]
  
  data("RandomVA1")
  names<-list(colnames(RandomVA1)[1:10], colnames(RandomVA1)[11:13],colnames(RandomVA1)[14:17],colnames(RandomVA1)[18:20],colnames(RandomVA1)[21:22],
              colnames(RandomVA1)[23],colnames(RandomVA1)[24],colnames(RandomVA1)[25],colnames(RandomVA1)[26],colnames(RandomVA1)[27],colnames(RandomVA1)[28],
              colnames(RandomVA1)[29],colnames(RandomVA1)[30],colnames(RandomVA1)[31],colnames(RandomVA1)[32],colnames(RandomVA1)[33],
              colnames(RandomVA1)[34],colnames(RandomVA1)[35],colnames(RandomVA1)[36],colnames(RandomVA1)[37],colnames(RandomVA1)[38],colnames(RandomVA1)[39:40],colnames(RandomVA1)[41:42],colnames(RandomVA1)[43:46],
              colnames(RandomVA1)[47:52],colnames(RandomVA1)[53:56],colnames(RandomVA1)[57:62],
              colnames(RandomVA1)[63:64], colnames(RandomVA1)[65],colnames(RandomVA1)[66:70],colnames(RandomVA1)[71:72],colnames(RandomVA1)[73:76],
              colnames(RandomVA1)[77:79], colnames(RandomVA1)[80:82],colnames(RandomVA1)[83],colnames(RandomVA1)[84:86],colnames(RandomVA1)[87:91],colnames(RandomVA1)[92:94],
              colnames(RandomVA1)[95:96],colnames(RandomVA1)[97:100],colnames(RandomVA1)[101:104],colnames(RandomVA1)[105:106],colnames(RandomVA1)[107:108],colnames(RandomVA1)[109:116],
              colnames(RandomVA1)[117],colnames(RandomVA1)[118],colnames(RandomVA1)[119],colnames(RandomVA1)[120:121],colnames(RandomVA1)[122:125],colnames(RandomVA1)[126],colnames(RandomVA1)[127:128],
              colnames(RandomVA1)[129:130],colnames(RandomVA1)[131:134],colnames(RandomVA1)[135:138],colnames(RandomVA1)[139], colnames(RandomVA1)[140:146],colnames(RandomVA1)[147], colnames(RandomVA1)[148:149],
              colnames(RandomVA1)[150:152],colnames(RandomVA1)[153:157],colnames(RandomVA1)[158:159],colnames(RandomVA1)[160:161],colnames(RandomVA1)[162:164],colnames(RandomVA1)[165:167],colnames(RandomVA1)[168:170],
              colnames(RandomVA1)[171:178],colnames(RandomVA1)[179:180],colnames(RandomVA1)[181:184],colnames(RandomVA1)[185:186],colnames(RandomVA1)[187:188],colnames(RandomVA1)[189:194],colnames(RandomVA1)[195:196],
              colnames(RandomVA1)[197],colnames(RandomVA1)[198:200],colnames(RandomVA1)[201:202],colnames(RandomVA1)[203:207],colnames(RandomVA1)[208:211],colnames(RandomVA1)[c(212, 215:217)],colnames(RandomVA1)[213:214],
              colnames(RandomVA1)[c(218, 222)],colnames(RandomVA1)[219:221],colnames(RandomVA1)[223],colnames(RandomVA1)[224],colnames(RandomVA1)[225],colnames(RandomVA1)[226],colnames(RandomVA1)[227],colnames(RandomVA1)[228:233],
              colnames(RandomVA1)[234:235],colnames(RandomVA1)[236],colnames(RandomVA1)[237:246])
  
  # Your existing list of grouped column names
  grouped_names <- names
  
  # Compute lengths of each group
  group_lengths <- lengths(grouped_names)
  
  # Create list of indices corresponding to each group
  index_list <- mapply(
    function(start, len) seq(from = start, length.out = len),
    start = cumsum(c(1, head(group_lengths, -1))),
    len = group_lengths,
    SIMPLIFY = FALSE
  )
  
  Index<-lapply(index_list, as.numeric)
  
  
  block_text<- data.frame(Names = I(names),
                          Index = I(Index))
  
  
  #Number of simulations
  n=100
  
  #Randomly assign pi
  random_vector <- runif(60)
  # Normalize to sum to 1
  pi <- random_vector / sum(random_vector)
  
  
  # Take known probbase (60 causes and first 245 questions)
  data(condprobnum)
  pb= condprobnum
  
  #Or random probbase
  #pb=matrix(runif(length(pi)*245),nrow=245,ncol=length(pi))
  #Different to previous attempts as we want input probBase to be full probbase(exactly same as condprobnum), not a subset of it(or letter probBase)
  
  #List blocks of questions which are conditionally independent-leaving out first block(age/sex questions)- needed for this version to work
  blocks=block_text$Index[-1]
  
  
  #Covariance function
  dcov=function(n,xc=0.1) diag(n) + (1-diag(n))*xc
  
  #Covariance matrix for each block of questions
  covs<-vector("list", length = length(blocks))
  for (k in 1:length(blocks)) {
    covs[[k]]<- dcov(length(blocks[[k]]),0.1)
  }
  
  
  #Age probs
  ageProb<- c(1/3,1/3,1/3,0,0,0,0)
  
  #Sex probs
  sexProb<-c(0.5,0.5)
  
  #Blocks of questions just for women deaths
  woman_block<-blocks[c(1,48,50:61)]
  woman_block_index<-c(1,48,50:61)
  woman_block_names<- block_text$Names[woman_block_index +1]
  
  #Blocks of questions just for neonate deaths
  neonate_block<-blocks[c(2,62:76)]
  neonate_block_index<-c(2,62:76)
  neonate_block_names<- block_text$Names[neonate_block_index +1]
  
  epsilon = 0.4

  

#Both work RandomVA1 and va_data(simulation) work individually(not in loop below)
  probBase<- condprobnum
  va_data<-simulateVA(n,pi,pb,blocks,covs,ageProb,sexProb,woman_block,neonate_block)
 out1<-gradientProbbase(va_data,blocks = blocks,epsilon = epsilon, probbase_element = list(10,25),
                                                   CondProbNum = probBase,model = "InSilicoVA")
probBase<- condprobnum
out2<- gradientProbbase(RandomVA1[1:100,],blocks = blocks,epsilon = epsilon, probbase_element = list(1,60),
                                                  CondProbNum = probBase,model = "InSilicoVA")

  
  #InSilico Simulation
  #Then for each element of probbase, simulate data and then find gradient:
probBase<- condprobnum
gradient_probbase_matrix<- matrix(0, nrow =length(rownames(probBase)), ncol = length(colnames(probBase)))

  
#ISSUE IN RUNNING - when using va_data- from simulateVA- Error in out$Imputed[k, j] : subscript out of bounds
#b starts at 10 to avoid sex/age questions
  for (a in 1:4) {
    for (b in 10:13) {
      va_data<-simulateVA(n,pi,pb,blocks,covs,ageProb,sexProb,woman_block,neonate_block)
      gradient_probbase_matrix[b,a]<- gradientProbbase(va_data,blocks = blocks,epsilon = epsilon, probbase_element = list(b,a),
                                                       CondProbNum = probBase,model = "InSilicoVA")
    }
  }
  
#Works for RandomVAA1- using subset of probBase to test as othewise too long to compute
  for (a in 1:4) {
    for (b in 10:13) {
      va_data<-simulateVA(n,pi,pb,blocks,covs,ageProb,sexProb,woman_block,neonate_block)
      gradient_probbase_matrix[b,a]<- gradientProbbase(RandomVA1[1:100,],blocks = blocks,epsilon = epsilon, probbase_element = list(b,a),
                                                       CondProbNum = probBase,model = "InSilicoVA")
    }
  }
  
  
  
  
  
#Try interVA


probBase<- probbase
va_data<-simulateVA(n,pi,pb,blocks,covs,ageProb,sexProb,woman_block,neonate_block)
out1<-gradientProbbase(va_data,blocks = blocks,epsilon = epsilon, probbase_element = list(10,25),
                       probBase = probBase,model = "InterVA")
#Doesn't work for simulated data

probBase<- probbase
out2<- gradientProbbase(RandomVA1[1:100,],blocks = blocks,epsilon = epsilon, probbase_element = list(10,25),
                        probBase = probBase,model = "InterVA")



#Evaluate probbase function for reference- as I think issue may be here

#Evaluating a probbase

#Function to take several different blocks, impute in turn, find accuracy, then average to get probbase overall
#accuracy
##' Evaluating a probbase for InSilico and InterVA models
##' @name evaluateProbbase
##' @description Takes a dataset of verbal autopsies, a set of multiple blocks of questions and a VA algorithm.
##' Sets the answers to those questions to 'missing', and attempts to impute their values for each block in turn.
##' Then calculates overall accuracy of the probbase by averaging imputed values across all answers set to missing. 
##' @param VA VA dataset, in the format of (e.g.) RandomVA1
##' @param blocks A list of blocks of questions. Each block is a set of indices for questions, corresponding to rows in VA
##' @param probBase Default null; uses data(probbase). Any custom probbase must be of exact same format as data(probbase). Only use this argument for custom probbase for InterVA model; see CondProbNum function for custom probbase for InSilicoVA model
##' @param letterProb A data frame with first coloumn letter and second coloumn the corresponding probabilities. Should only be specified if using InterVA model and want to change default. If left null then will use the following default: data.frame(grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),value = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
##' @param probbase_element An element of the probbase you want to change by epsilon indexed as a list
##' @param epsilon A small value you want to add to the selected probbase element. If causes probability to go out of range (0,1) then reverse operation(addition or subtraction) will be made so within the range. Absolute value of epsilon should be <=0.1. Epsilon can be positive or negative.
##' @param method Method for attaining cause-of-death distribution from VA answers; default 'OpenVA'
##' @param data.type Format of the VA data; default 'WHO2012'
##' @param model VA algorithm used; at the moment supports 'InterVA' and "InSilicoVA"
##' @param CondProbNum A customised probability matrix with 245 row symptomns on 60 column causes in the same order as InterVA4 specification. For example input see condprobnum. Used for InSilico method only
##' @param version Version of InterVA used; default '4.03'
##' @param HIV Set as 'h','m','l' for InterVA
##' @param Malaria Set as 'h','m','l' for InterVA
##' @return Accuracy of probbase
##' @export
##' @examples
##' ## Impute blocks of questions 18-20 and 21-24 and evaluate probbase for InterVA(used fixed probbase) and InSilicoVA(using condprobnum)
##' 
##' data(RandomVA1)
##' evaluateProbbase(RandomVA1[1:100,],blocks=list(18:20, 21:24),method="OpenVA",data.type = "WHO2012",
##' model = "InterVA",version = "4.03", HIV = "h", Malaria = "h")
##' 
##' data("condprobnum")
##' evaluateProbbase(RandomVA1[1:100,],blocks=list(18:20, 21:24),method="OpenVA",data.type = "WHO2012",
##' model = "InSilicoVA",CondProbNum=condprobnum)
##' 
##' 
##' ## Impute blocks of questions 18-20 and 21-24 and evaluate probbase for InterVA(used custom probbase and letter probabilities)
##' data("probbase")
##' probbase2<- probbase
##' probbase2[10,22]<-"A"
##' evaluateProbbase(RandomVA1[1:100,],blocks=list(18:20, 21:24),method="OpenVA",data.type = "WHO2012",
##' model = "InterVA",version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))
##' 
##' 
##' data("probbase")
##' probbase2<- probbase
##' probbase2[10,22]<-"A"
##' evaluateProbbase(RandomVA1[1:10,],blocks=list(18:20, 21:24),method="OpenVA",data.type = "WHO2012",
##' model = "InterVA",version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)), probbase_element =list(10,25) , epsilon = 0.4)





                 
evaluateProbbase<-function(VA,blocks, method = "OpenVA",data.type = "WHO2012", model, CondProbNum=NULL, probBase=NULL, letterProb=NULL, probbase_element=NULL, epsilon=NULL,
                                version = NULL, HIV = NULL, Malaria = NULL){
  #For each block in list of block, run recover answers, extract imputed values where we removed a Yes
  #Once complete take mean of imputed values
  library(openVA)
  library(nbc4va)
  count<-0
  sum<- 0
  for(i in 1:length(blocks)){
    block<- blocks[[i]]
    if(model=="InterVA"){
      out=recoverAnswers3(VA,block,method=method,
                          data.type = data.type,model = model,
                          version = version, HIV = HIV, Malaria = Malaria, probBase=probBase, letterProb=letterProb,
                          probbase_element = probbase_element, epsilon = epsilon)
    }
    if(model=="InSilicoVA"){
      out=recoverAnswers2(VA,block,method=method,
                          data.type = data.type,model = model,
                          CondProbNum=CondProbNum)
    }
    if(length(block) == 1){
      for (k in 1:length(out$Imputed[,1])) {
        row<- out$Original[k]
        for (j in 1:length(out$Imputed[1,])) {
          if(row[j]=="Y"){
            count<- count+1
            sum<- sum + out$Imputed[k]
          }
        }
      }
    }
    else{
      for (k in 1:length(out$Original[,1])) {
        row<- out$Original[k,]
        for (j in 1:length(out$Original)) {
          if(row[j]=="Y"){
            count<- count+1
            sum<- sum + out$Imputed[k,j]
          }
        }
      }
      }
  return(Accuracy=sum/count)
  }
}


