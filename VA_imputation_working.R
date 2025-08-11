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
#Simulating VA including age and sex
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
#Custom InterVA- allow custom probbase/letter probabilities
##**********************************************************************
##


#Allow custom probbase/letter probabilities for InterVA2- working

##' Provide InterVA4 analysis on the data input(with custom probbase/letter probability options)
##' @name InterVA2
##' @description This function implements the algorithm in the InterVA4 software. 
##' It produces individual cause of death and population cause-specific mortality fractions. Allows for customised probbase of the exact form of data(probbase) and also allows custom letter probabilities for the probbase.
##' @param probBase Default null; uses data(probbase). Any custom probbase must be of exact same format as data(probbase).
##' @param letterProb A data frame with first coloumn letter and second coloumn the corresponding probabilities. If left null then will use the following default: data.frame(grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),value = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
##' @param Input A matrix input, or data read from csv files in the same format as required by InterVA4. Sample input is included as data(SampleInput).
##' @param HIV An indicator of the level of prevalence of HIV. The input should be one of the following: "h"(high),"l"(low), or "v"(very low).
##' @param Malaria An indicator of the level of prevalence of Malaria. The input should be one of the following: "h"(high),"l"(low), or "v"(very low).
##' @param directory The directory to store the output from InterVA4. It should either be an existing valid directory, or a new folder to be created. If no path is given, the current working directory will be used.
##' @param filename The filename the user wish to save the output. No extension needed. The output is in .csv format by default.
##' @param output "classic": The same deliminated output format as InterVA4; or "extended": deliminated output followed by full distribution of cause of death proability.
##' @param append A logical value indicating whether or not the new output should be appended to the existing file.
##' @param groupcode A logical value indicating whether or not the group code will be included in the output causes.
##' @param replicates A logical value indicating whether or not the calculation should replicate original InterVA4 software (version 4.02) exactly. If replicates = F, causes with small probability are not dropped out of calculation in intermediate steps, and a possible bug in original InterVA4 implementation is fixed. If replicates=T, then the output values will be exactly as they would be from calling the InterVA4 program (version 4.02). If replicates=F, the output values will be the same as calling the InterVA4 program (version 4.03). Since version 1.7.3, setting replicates to be FALSE also includes changes to data checking rules and pre-set conditional probabilities to be the same as the official version 4.03 software. Since version 1.6, two control variables are added to control the two bugs respectively. Setting this to TRUE will overwrite both to TRUE.
##' @param replicate.bug1 This logical indicator controls whether or not the bug in InterVA4.2 involving the symptom "skin_les" will be replicated or not. It is suggested to set to FALSE.
##' @param replicate.bug2 This logical indicator controls whether the causes with small probability are dropped out of calculation in intermediate steps or not. It is suggested to set to FALSE.
##' @param write A logical value indicating whether or not the output (including errors and warnings) will be saved to file.
##' @param ... not used
##' @return Accuracy of probbase
##' @export
##' @examples
##' 
##' data(SampleInput)
##' to get easy-to-read version of causes of death make sure the column
##' orders match interVA4 standard input this can be monitored by checking
##' the warnings of column names
##' 
##' 
##' #Use a custom probbase
##' data("probbase")
##' probbase2<- probbase
##' probbase2[,19]<-"A"
##' 
##' 
##' #Custom letter probabilities
##' sample.output1 <- InterVA2(SampleInput, HIV = "h", Malaria = "l", write=FALSE, 
##' replicates = FALSE, probBase=probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))
##' 
##' 
##' to get causes of death with group code for further usage
##' sample.output2 <- InterVA2(SampleInput, HIV = "h", Malaria = "l", write=FALSE,
##' replicates = FALSE, groupcode = TRUE, probBase=probbase2, letterProb = NULL)
##' 
##' 
##' Using default probbase and letter probabilities
##' sample.output3 <- InterVA2(SampleInput, HIV = "h", Malaria = "l", write=FALSE,
##' replicates = FALSE, groupcode = TRUE, probBase=NULL, letterProb = NULL)



InterVA2<-function (Input, HIV, Malaria, directory = NULL, filename = "VA_result", 
          output = "classic", append = FALSE, groupcode = FALSE, replicates = FALSE,probBase=NULL, letterProb=NULL,
          replicate.bug1 = FALSE, replicate.bug2 = FALSE, write = TRUE, 
          ...)
{
  va <- function(ID, MALPREV, HIVPREV, PREGSTAT, PREGLIK, PRMAT, 
                 INDET, CAUSE1, LIK1, CAUSE2, LIK2, CAUSE3, LIK3, wholeprob, 
                 ...) {
    ID <- ID
    MALPREV <- as.character(MALPREV)
    HIVPREV <- as.character(HIVPREV)
    PREGSTAT <- paste(PREGSTAT, paste(rep(" ", 5 - nchar(PREGSTAT)), 
                                      collapse = ""), collapse = "")
    PREGLIK <- PREGLIK
    PRMAT <- PRMAT
    INDET <- as.character(INDET)
    wholeprob <- wholeprob
    va.out <- list(ID = ID, MALPREV = MALPREV, HIVPREV = HIVPREV, 
                   PREGSTAT = PREGSTAT, PREGLIK = PREGLIK, PRMAT = PRMAT, 
                   INDET = INDET, CAUSE1 = CAUSE1, LIK1 = LIK1, CAUSE2 = CAUSE2, 
                   LIK2 = LIK2, CAUSE3 = CAUSE3, LIK3 = LIK3, wholeprob = wholeprob)
    va.out
  }
  save.va <- function(x, filename, write) {
    if (!write) {
      return()
    }
    x <- x[-14]
    x <- as.matrix(x)
    filename <- paste(filename, ".csv", sep = "")
    write.table(t(x), file = filename, sep = ",", append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  save.va.prob <- function(x, filename, write) {
    if (!write) {
      return()
    }
    prob <- unlist(x[14])
    x <- x[-14]
    x <- unlist(c(as.matrix(x), as.matrix(prob)))
    filename <- paste(filename, ".csv", sep = "")
    write.table(t(x), file = filename, sep = ",", append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
  if (replicates) {
    warning("option 'replicates' is turned on, all bugs in InterVA-4 is replicated\n", 
            immediate. = TRUE)
    replicate.bug1 <- TRUE
    replicate.bug2 <- TRUE
  }
  if (is.null(directory)) 
    directory = getwd()
  dir.create(directory, showWarnings = FALSE)
  globle.dir <- getwd()
  setwd(directory)
  #NEW SECTION
  if(is.null(probBase)){
  data("probbase", envir = environment())
  probbase <- get("probbase", envir = environment())
  }
  if(!is.null(probBase)){
    probbase<- probBase
  }
  if (!replicates) {
    if(!is.null(probBase)){
      probbase<- probBase
    }
    if(is.null(probBase)){
      data("probbase3", envir = environment())
      probbase <- get("probbase3", envir = environment())
    }
  }
  #END NEW SECTION
  probbase <- as.matrix(probbase)
  data("causetext", envir = environment())
  causetext <- get("causetext", envir = environment())
  if (groupcode) {
    causetext <- causetext[, -2]
  }
  else {
    causetext <- causetext[, -3]
  }
  if (write) {
    cat(paste("Error log built for InterVA", Sys.time(), 
              "\n"), file = "errorlog.txt", append = FALSE)
    cat(paste("Warning log built for InterVA", Sys.time(), 
              "\n"), file = "warnings.txt", append = FALSE)
  }
  Input <- as.matrix(Input)
  if (dim(Input)[1] < 1) {
    stop("error: no data input")
  }
  N <- dim(Input)[1]
  S <- dim(Input)[2]
  if (S != dim(probbase)[1]) {
    stop("error: invalid data input format. Number of values incorrect")
  }
  if (tolower(colnames(Input)[S]) != "scosts") {
    stop("error: the last variable should be 'scosts'")
  }
  data("SampleInput", envir = environment())
  SampleInput <- get("SampleInput", envir = environment())
  valabels = colnames(SampleInput)
  count.changelabel = 0
  for (i in 1:S) {
    if (tolower(colnames(Input)[i]) != tolower(valabels)[i]) {
      warning(paste("Input column '", colnames(Input)[i], 
                    "' does not match InterVA standard: '", valabels[i], 
                    "'", sep = ""), call. = FALSE, immediate. = TRUE)
      count.changelabel = count.changelabel + 1
    }
  }
  if (count.changelabel > 0) {
    warning(paste(count.changelabel, "column names changed in input. \n If the change in undesirable, please change in the input to match standard InterVA4 input format."), 
            call. = FALSE, immediate. = TRUE)
    colnames(Input) <- valabels
  }
  #NEW SECTION
  if(is.null(letterProb)){
    letterProbs <- data.frame(Letter = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
                                         "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
                              Prob = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,
                                       0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
  }
  if(!is.null(letterProb)){
    letterProbs <- letterProb
  }
  probbase[probbase == "I"] <- letterProbs[1,2]
  probbase[probbase == "A+"] <- letterProbs[2,2]
  probbase[probbase == "A"] <- letterProbs[3,2]
  probbase[probbase == "A-"] <- letterProbs[4,2]
  probbase[probbase == "B+"] <- letterProbs[5,2]
  probbase[probbase == "B"] <- letterProbs[6,2]
  probbase[probbase == "B-"] <- letterProbs[7,2]
  probbase[probbase == "B -"] <- letterProbs[8,2]
  probbase[probbase == "C+"] <- letterProbs[9,2]
  probbase[probbase == "C"] <- letterProbs[10,2]
  probbase[probbase == "C-"] <- letterProbs[11,2]
  probbase[probbase == "D+"] <- letterProbs[12,2]
  probbase[probbase == "D"] <- letterProbs[13,2]
  probbase[probbase == "D-"] <- letterProbs[14,2]
  probbase[probbase == "E"] <- letterProbs[15,2]
  probbase[probbase == "N"] <- letterProbs[16,2]
  probbase[probbase == ""] <- letterProbs[17,2]
  #END NEW SECTION
  probbase[1, 1:13] <- rep(0, 13)
  Sys_Prior <- as.numeric(probbase[1, ])
  D <- length(Sys_Prior)
  HIV <- tolower(HIV)
  Malaria <- tolower(Malaria)
  if (!(HIV %in% c("h", "l", "v")) || !(Malaria %in% c("h", 
                                                       "l", "v"))) {
    stop("error: the HIV and Malaria indicator should be one of the three: 'h', 'l', and 'v'")
  }
  if (HIV == "h") 
    Sys_Prior[19] <- 0.05
  if (HIV == "l") 
    Sys_Prior[19] <- 0.005
  if (HIV == "v") 
    Sys_Prior[19] <- 1e-05
  if (Malaria == "h") {
    Sys_Prior[21] <- 0.05
    Sys_Prior[39] <- 0.05
  }
  if (Malaria == "l") {
    Sys_Prior[21] <- 0.005
    Sys_Prior[39] <- 1e-05
  }
  if (Malaria == "v") {
    Sys_Prior[21] <- 1e-05
    Sys_Prior[39] <- 1e-05
  }
  ID.list <- rep(NA, N)
  VAresult <- vector("list", N)
  if (write && append == FALSE) {
    header = c("ID", "MALPREV", "HIVPREV", "PREGSTAT", "PREGLIK", 
               "PRMAT", "INDET", "CAUSE1", "LIK1", "CAUSE2", "LIK2", 
               "CAUSE3", "LIK3")
    if (output == "extended") 
      header = c(header, as.character(causetext[, 2]))
    write.table(t(header), file = paste(filename, ".csv", 
                                        sep = ""), row.names = FALSE, col.names = FALSE, 
                sep = ",")
  }
  nd <- max(1, round(N/100))
  np <- max(1, round(N/10))
  for (i in 1:N) {
    if (i%%nd == 0) {
      cat(".")
    }
    if (i%%np == 0) {
      cat(paste(round(i/N * 100), "% completed\n", sep = ""))
    }
    index.current <- as.character(Input[i, 1])
    Input[i, which(is.na(Input[i, ]))] <- "0"
    Input[i, which(toupper(Input[i, ]) != "Y")] <- "0"
    Input[i, which(toupper(Input[i, ]) == "Y")] <- "1"
    input.current <- as.numeric(Input[i, ])
    input.current[1] <- 0
    if (sum(input.current[2:8]) < 1) {
      if (write) {
        cat(paste(index.current, " Error in age indicator: Not Specified ", 
                  "\n"), file = "errorlog.txt", append = TRUE)
      }
      next
    }
    if (sum(input.current[9:10]) < 1) {
      if (write) {
        cat(paste(index.current, " Error in sex indicator: Not Specified ", 
                  "\n"), file = "errorlog.txt", append = TRUE)
      }
      next
    }
    if (sum(input.current[23:223]) < 1) {
      if (write) {
        cat(paste(index.current, " Error in indicators: No symptoms specified ", 
                  "\n"), file = "errorlog.txt", append = TRUE)
      }
      next
    }
    for (k in 1:2) {
      for (j in 1:(S - 1)) {
        if (input.current[j + 1] == 1) {
          Dont.ask <- probbase[j + 1, 4:11]
          Dont.ask.list <- input.current[match(toupper(Dont.ask), 
                                               toupper(colnames(Input)))]
          Dont.ask.list[is.na(Dont.ask.list)] <- 0
          if (sum(Dont.ask.list) > 0) {
            input.current[j + 1] <- 0
            if (write) {
              cat(index.current, "   ", paste(probbase[j + 
                                                         1, 2], "  value inconsistent with ", 
                                              Dont.ask[which(Dont.ask.list > 0)], " - cleared in working file \n"), 
                  file = "warnings.txt", append = TRUE)
            }
          }
        }
        if (input.current[j + 1] == 1) {
          Ask.if <- probbase[j + 1, 12]
          if (!is.na(match(toupper(Ask.if), toupper(colnames(Input))))) {
            if (input.current[match(toupper(Ask.if), 
                                    toupper(colnames(Input)))] == 0) {
              input.current[match(toupper(Ask.if), toupper(colnames(Input)))] <- 1
              if (write) {
                cat(index.current, "   ", paste(probbase[j + 
                                                           1, 2], "  not flagged in category ", 
                                                Ask.if, " - updated in working file \n"), 
                    file = "warnings.txt", append = TRUE)
              }
            }
          }
        }
      }
    }
    if (replicate.bug1 == TRUE && input.current[84] == 1) {
      input.current[85] <- 1
    }
    reproductiveAge <- 0
    preg_state <- " "
    lik.preg <- " "
    if (input.current[10] == 1 && (input.current[4] == 1 || 
                                   input.current[5] == 1)) 
      reproductiveAge <- 1
    prob <- Sys_Prior[14:D]
    temp <- which(input.current[2:length(input.current)] == 
                    1)
    for (jj in 1:length(temp)) {
      temp_sub <- temp[jj]
      for (j in 14:D) {
        prob[j - 13] <- prob[j - 13] * as.numeric(probbase[temp_sub + 
                                                             1, j])
      }
      if (sum(prob[1:3]) > 0) 
        prob[1:3] <- prob[1:3]/sum(prob[1:3])
      if (sum(prob[4:63]) > 0) 
        prob[4:63] <- prob[4:63]/sum(prob[4:63])
      if (replicate.bug2) {
        prob[prob < 1e-06] <- 0
      }
    }
    names(prob) <- causetext[, 2]
    prob_A <- prob[1:3]
    prob_B <- prob[4:63]
    if (sum(prob_A) == 0 || reproductiveAge == 0) {
      preg_state <- "Indet"
      lik.preg <- 0
    }
    if (which.max(prob_A) == 1 && prob_A[1] != 0 && reproductiveAge == 
        1) {
      preg_state <- "nrp"
      lik.preg <- as.numeric(round(prob_A[1]/sum(prob_A) * 
                                     100))
    }
    if (which.max(prob_A) == 2 && prob_A[2] != 0 && reproductiveAge == 
        1) {
      preg_state <- "pr6w"
      lik.preg <- as.numeric(round(prob_A[2]/sum(prob_A) * 
                                     100))
    }
    if (which.max(prob_A) == 3 && prob_A[3] != 0 && reproductiveAge == 
        1) {
      preg_state <- "preg"
      lik.preg <- as.numeric(round(prob_A[3]/sum(prob_A) * 
                                     100))
    }
    lik_mat <- " "
    if (reproductiveAge == 1 && sum(prob_A) != 0) 
      lik_mat <- as.numeric(round((prob_A[2] + prob_A[3])/sum(prob_A) * 
                                    100))
    if (sum(prob_B) != 0) 
      prob_B <- prob_B/sum(prob_B)
    prob.temp <- prob_B
    if (max(prob.temp) <= 0.4) {
      indet <- "Indet"
      cause1 <- lik1 <- cause2 <- lik2 <- cause3 <- lik3 <- " "
    }
    if (max(prob.temp) > 0.4) {
      indet <- " "
      lik1 <- round(max(prob.temp) * 100)
      cause1 <- names(prob.temp)[which.max(prob.temp)]
      prob.temp <- prob.temp[-which.max(prob.temp)]
      lik2 <- round(max(prob.temp) * 100)
      cause2 <- names(prob.temp)[which.max(prob.temp)]
      if (max(prob.temp) < 0.5 * max(prob_B)) 
        lik2 <- cause2 <- " "
      prob.temp <- prob.temp[-which.max(prob.temp)]
      lik3 <- round(max(prob.temp) * 100)
      cause3 <- names(prob.temp)[which.max(prob.temp)]
      if (max(prob.temp) < 0.5 * max(prob_B)) 
        lik3 <- cause3 <- " "
    }
    ID.list[i] <- index.current
    VAresult[[i]] <- va(ID = index.current, MALPREV = Malaria, 
                        HIVPREV = HIV, PREGSTAT = preg_state, PREGLIK = lik.preg, 
                        PRMAT = lik_mat, INDET = indet, CAUSE1 = cause1, 
                        LIK1 = lik1, CAUSE2 = cause2, LIK2 = lik2, CAUSE3 = cause3, 
                        LIK3 = lik3, wholeprob = c(prob_A, prob_B))
    if (output == "classic") 
      save.va(VAresult[[i]], filename = filename, write)
    if (output == "extended") 
      save.va.prob(VAresult[[i]], filename = filename, 
                   write)
  }
  setwd(globle.dir)
  out <- list(ID = ID.list[which(!is.na(ID.list))], VA = VAresult[which(!is.na(ID.list))], 
              Malaria = Malaria, HIV = HIV)
  class(out) <- "interVA"
  return(out)
}





##**********************************************************************
#Custom codeVA- allow custom probbase/letter probabilities for interVA model within codeVA
##**********************************************************************
##


#Allow custom probbase/letter probabilities for InterVA using codeVA function-working


##' Version of codeVA which provides InterVA4 analysis with custom probbase/letter probability options
##' @name codeVA2
##' @description This function implements the codeVA function but allows for customised probbase of the exact form of data(probbase) and also allows custom letter probabilities for the probbase when using InterVA.
##' @param probBase Default null; uses data(probbase). Any custom probbase must be of exact same format as data(probbase). Only use custom probbase for InterVA model; see insilico function for custom probbase for InSilicoVA model
##' @param letterProb A data frame with first coloumn letter and second coloumn the corresponding probabilities. Should only be specified if using InterVA model and want to change default. If left null then will use the following default: data.frame(grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),value = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
##' @param Input A matrix input, or data read from csv files in the same format as required by InterVA4. Sample input is included as data(SampleInput).
##' @param HIV An indicator of the level of prevalence of HIV. The input should be one of the following: "h"(high),"l"(low), or "v"(very low).
##' @param Malaria An indicator of the level of prevalence of Malaria. The input should be one of the following: "h"(high),"l"(low), or "v"(very low).
##' @param directory The directory to store the output from InterVA4. It should either be an existing valid directory, or a new folder to be created. If no path is given, the current working directory will be used.
##' @param filename The filename the user wish to save the output. No extension needed. The output is in .csv format by default.
##' @param output "classic": The same deliminated output format as InterVA4; or "extended": deliminated output followed by full distribution of cause of death proability.
##' @param append A logical value indicating whether or not the new output should be appended to the existing file.
##' @param groupcode A logical value indicating whether or not the group code will be included in the output causes.
##' @param replicates A logical value indicating whether or not the calculation should replicate original InterVA4 software (version 4.02) exactly. If replicates = F, causes with small probability are not dropped out of calculation in intermediate steps, and a possible bug in original InterVA4 implementation is fixed. If replicates=T, then the output values will be exactly as they would be from calling the InterVA4 program (version 4.02). If replicates=F, the output values will be the same as calling the InterVA4 program (version 4.03). Since version 1.7.3, setting replicates to be FALSE also includes changes to data checking rules and pre-set conditional probabilities to be the same as the official version 4.03 software. Since version 1.6, two control variables are added to control the two bugs respectively. Setting this to TRUE will overwrite both to TRUE.
##' @param replicate.bug1 This logical indicator controls whether or not the bug in InterVA4.2 involving the symptom "skin_les" will be replicated or not. It is suggested to set to FALSE.
##' @param replicate.bug2 This logical indicator controls whether the causes with small probability are dropped out of calculation in intermediate steps or not. It is suggested to set to FALSE.
##' @param write A logical value indicating whether or not the output (including errors and warnings) will be saved to file.
##' @param ... not used
##' @return Accuracy of probbase
##' @export
##' @examples
##' 
##' data(RandomVA3)test <- RandomVA3[1:200, ]
##' train <- RandomVA3[201:400, ]
##' 
##' 
##' #Use for other models same as codeVA
##' fit1 <- codeVA2(data = test, data.type = "customize", model = "InSilicoVA",
##' data.train = train, causes.train = "cause",
##' Nsim=1000, auto.length = FALSE)
##' 
##' 
##' #Custom probbase/letter probabilities for InterVA
##' data("probbase")
##' probbase2<- probbase
##' probbase2[,19]<-"A"
##' fit2 <- codeVA2(data = test, data.type = "customize", model = "InterVA",
##' data.train = train, causes.train = "cause", write=FALSE,
##' version = "4.02", HIV = "h", Malaria = "l", probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))
##' 
##' 
##' #Default probbase and letter probabilities for InterVA
##' fit2 <- codeVA2(data = test, data.type = "customize", model = "InterVA",
##' data.train = train, causes.train = "cause", write=FALSE,
##' version = "4.02", HIV = "h", Malaria = "l")


codeVA2<-function (data, data.type = c("WHO2012", "WHO2016", "PHMRC", 
                              "customize")[2], data.train = NULL, causes.train = NULL, probBase=NULL, letterProb=NULL,
          causes.table = NULL, model = c("InSilicoVA", "InterVA", "Tariff", 
                                         "NBC")[1], Nchain = 1, Nsim = 10000, version = c("4.02", 
                                                                                          "4.03", "5")[2], HIV = "h", Malaria = "h", phmrc.type = c("adult", 
                                                                                                                                                    "child", "neonate")[1], convert.type = c("quantile", 
                                                                                                                                                                                             "fixed", "empirical")[1], ...) 
{
  version <- as.character(version)
  if (version == "5.0") 
    version <- "5"
  args <- as.list(match.call())
  if (data.type == "WHO") {
    data.type <- "WHO2012"
    warning("The argument data.type of 'WHO' is no longer in use. 'WHO2012' or 'WHO2016' needs to be specified. Default change to 'WHO2012' for backward compatibility.\n", 
            immediate. = TRUE)
    args$data.type <- "WHO2012"
  }
  if (data.type == "WHO2016" & model == "InterVA" & version != 
      "5") {
    stop("Error: WHO2016 type input does not work with InterVA 4.02 or 4.03. Consider switching to 5")
  }
  if (data.type == "WHO2012" & model == "InterVA" & version == 
      "5") {
    stop("Error: WHO2012 type input does not work with InterVA 5. Consider switching to 4.03")
  }
  if (data.type %in% c("WHO2012", "WHO2016") && (model == "Tariff" || 
                                                 model == "NBC")) {
    if (is.null(data.train) || is.null(causes.train)) {
      stop("Error: need training data for WHO questionnaire input with Tariff or NBC method.")
    }
  }
  if (data.type == "customize") {
    if (is.null(data.train) || is.null(causes.train)) {
      stop("Error: need training data for customized input.")
    }
    tmp <- data[, -1]
    tmp <- tmp[, colnames(tmp) != causes.train]
    tmp <- toupper(as.character(as.matrix(tmp)))
    if (sum(!tmp %in% c("Y", "", "N", ".", "-")) != 0) {
      stop("Error: customized train/test data need to use ``Y`` to denote ``presence'', ``'' to denote ``absence'', and ``.'' to denote ``missing''.")
    }
  }
  if (data.type == "PHMRC") {
    if (is.null(data.train)) {
      stop("Error: need training data for PHMRC data, possible training data could be obtained at http://ghdx.healthdata.org/record/population-health-metrics-research-consortium-gold-standard-verbal-autopsy-data-2005-2011")
    }
    if (is.null(causes.train)) {
      stop("Error: please specify which column is the cause-of-death in PHMRC input")
    }
    binary <- ConvertData.phmrc(input = data.train, input.test = data, 
                                cause = causes.train, phmrc.type = phmrc.type, convert.type = convert.type, 
                                ...)
    data.train <- binary$output
    data <- binary$output.test
    causes.train <- colnames(data.train)[2]
  }
  if (model == "InSilicoVA") {
    if (is.null(Nsim)) {
      stop("Please specify Nsim: number of iterations to draw from InSilicoVA sampler")
    }
    if (is.null(args$warning.write) && !(is.null(args$write))) {
      args$warning.write <- args$write
    }
    if (is.null(args$burnin)) {
      args$burnin <- round(Nsim/2)
    }
    if (is.null(args$thin)) {
      args$thin <- 10 + 10 * (Nsim >= 10000)
    }
    if (is.null(args$Nsim)) {
      args$Nsim <- Nsim
    }
    if (data.type %in% c("WHO2012", "WHO2016")) {
      fit <- do.call("insilico", pairlist(args)[[1]][-1])
    }
    else if (data.type == "PHMRC" || data.type == "customize") {
      args$data <- as.name("data")
      args$train <- as.name("data.train")
      args$cause <- as.name("causes.train")
      args$type <- convert.type
      fit <- do.call("insilico.train", pairlist(args)[[1]][-1])
    }
    else {
      stop("Error: unknown data type specified")
    }
  }
  #NEW SECTION
  else if (model == "InterVA") {
    if (version == "4.02") {
      replicates = TRUE
    }
    else {
      replicates = FALSE 
    }
    if (data.type == "WHO2012") {
      if (is.null(args$write)) {
        args$write <- FALSE
      }
      fit <- calibrateVA::InterVA2(Input = data, HIV = HIV,probBase=probBase, letterProb = letterProb,
                               Malaria = Malaria, replicates = replicates, ...)
      #END NEW SECTION
    }
    else if (data.type == "WHO2016") {
      if (is.null(args$write)) {
        args$write <- FALSE
      }
      for (i in 1:dim(data)[2]) {
        data[, i] <- as.character(data[, i])
        data[, i][data[, i] == ""] <- "n"
      }
      tmp <- tolower(as.character(as.matrix(data[, -1])))
      if (sum(tmp %in% c("y", "n", "-", ".")) < length(tmp)) {
        stop("InterVA5 input data contains values other than 'y', 'n', '.', or '-'. Please check your input, especially for extra space characters in the cells, or standardize how missing data is coded.")
      }
      fit <- InterVA5::InterVA5(Input = data, HIV = HIV, 
                                Malaria = Malaria, ...)
    }
    else if (data.type == "PHMRC" || data.type == "customize") {
      fit <- interVA_train(data = data, train = data.train, 
                           causes.train = causes.train, causes.table = causes.table, 
                           type = convert.type, ...)
    }
    else {
      stop("Error: unknown data type specified")
    }
  }
  else if (model == "Tariff") {
    if (data.type == "WHO2016") {
      data <- ConvertData(data, yesLabel = c("y", "Y"), 
                          noLabel = c("n", "N"), missLabel = c("-"))
      data.train <- ConvertData(data.train, yesLabel = c("y", 
                                                         "Y"), noLabel = c("n", "N"), missLabel = c("-"))
    }
    data.train[, causes.train] <- as.character(data.train[, 
                                                          causes.train])
    if (data.type %in% c("WHO2012", "WHO2016")) {
      fit <- Tariff::tariff(causes.train = causes.train, 
                            symps.train = data.train, symps.test = data, 
                            causes.table = NULL, ...)
    }
    else if (data.type == "PHMRC" || data.type == "customize") {
      fit <- Tariff::tariff(causes.train = causes.train, 
                            symps.train = data.train, symps.test = data, 
                            causes.table = NULL, ...)
    }
    else {
      stop("Error: unknown data type specified")
    }
  }
  else if (model == "NBC") {
    if (!isTRUE(requireNamespace("nbc4va", quietly = TRUE))) {
      stop("You need to install the packages 'nbc4va'. Please run in your R terminal:\n install.packages('nbc4va')")
    }
    if (data.type == "WHO2016") {
      data <- ConvertData(data, yesLabel = c("y", "Y"), 
                          noLabel = c("n", "N"), missLabel = c("-"))
    }
    data.train[, 1] <- as.character(data.train[, 1])
    data[, 1] <- as.character(data[, 1])
    fit <- nbc4va::ova2nbc(data.train, data, causes.train)
    if (colnames(fit$prob)[1] != "CaseID") {
      temp <- data.frame(CaseID = fit$test.ids)
      fit$prob <- cbind(temp, fit$prob)
    }
    for (i in 1:dim(fit$prob)[1]) {
      if (sum(fit$prob[i, -1]) > 0) {
        fit$prob[i, -1] <- fit$prob[i, -1]/sum(fit$prob[i, 
                                                        -1])
      }
    }
  }
  else {
    stop("Error, unknown model specification")
  }
  return(fit)
}
