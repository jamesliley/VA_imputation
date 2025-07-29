##**********************************************************************
##  Calibration of symptom frequencies in VA using imputation       ####
##   proposal. 
##   
##  Experimentation and testing
##  
##  Nathan Higgins, James Liley, Eilidh Cowan                                                        
##  July 2025                                                              
##**********************************************************************
##
## Run with working directory set to main directory for repository.
## This directory should contain the document VA_imputation_functions.R
## 


##**********************************************************************
## Packages and scripts                                             ####
##**********************************************************************

# Packages
library(openVA)     # General verbal autopsy functions including openVA and insilicoVA
library(nbc4va)     # Naive Bayes classiier

# Scripts
source("./VA_imputation_functions.R")   # Functions


##**********************************************************************
## Seeds and switches                                               ####
##**********************************************************************

seed=382839                     # Uncomment for repeatable simulation
# seed=as.numeric(Sys.time())   # Uncomment for randomised simulation
set.seed(seed)

# Output directory: for figures, tables etc
output_dir="./Outputs/"



##**********************************************************************
## Datasets                                                         ####
##**********************************************************************

# Get dataset RandomVA1 (1000 VAs) from package
data(RandomVA1)


##**********************************************************************
## To be organised                                                  ####
##**********************************************************************


## Brief experiment to look at changes in VA output when some data is set to NA

#Now we are going to group columns 2-8(age realted questions)- only used column 2 as treating 7 columns
#missing currently causing issues
#We will create new data set which turns the answers in column to missing and then compute new posteriror CSMF


new_data<- RandomVA1
new_data[, 11:20]<- "."

fit_inter_who_old <- codeVA(data = RandomVA1, data.type = "WHO2012",
                            model = "InterVA", version = "4.03",
                            HIV = "h", Malaria = "h")
fit_inter_who_new <- codeVA(data = new_data, data.type = "WHO2012",
                            model = "InterVA", version = "4.03",
                            HIV = "h", Malaria = "h")
to_prob=function(vax) {
probs=as.data.frame(
  t(as.matrix(as.data.frame(
    lapply(vax$VA,function(x) x$wholeprob)))))
rownames(probs)=vax$ID
colnames(probs)=names(vax$VA[[1]]$wholeprob)
xprobs=probs[,4:dim(probs)[2]]
return(xprobs)
}

x1=to_prob(fit_inter_who_old)
x2=to_prob(fit_inter_who_new)
## plot(rowSums(xprobs))

fit_ins_who_new <- codeVA(new_data, data.type = "WHO2012", model = "InSilicoVA",
                          Nsim = 10000, auto.length = FALSE)

#Posterior CSMFs
csmf_inter_new <- getCSMF(fit_inter_who_new)
csmf_ins_new <- getCSMF(fit_ins_who_new)






## Function to impute missing VA values. ####
## To be migrated to VA_imputation_functions.R once working.


#Note function should only take single VA interview in argument

my_function<- function(VA,block,method="OpenVA",data.type = "WHO2012",
                       model = "InterVA", version = "4.03", HIV = "h",
                       Malaria = "h"){
  new_data<- VA
  missing<- VA[,block]
  new_data[, block]<- ""
  
  #This is probbase for InterVA4.03, need to convert into matrix of probs
  data(probbase, package = "InterVA4")
  
  #from inside interVA function
  
  probbase[probbase == "I"] <- 1
  probbase[probbase == "A+"] <- 0.8
  probbase[probbase == "A"] <- 0.5
  probbase[probbase == "A-"] <- 0.2
  probbase[probbase == "B+"] <- 0.1
  probbase[probbase == "B"] <- 0.05
  probbase[probbase == "B-"] <- 0.02
  probbase[probbase == "B -"] <- 0.02
  probbase[probbase == "C+"] <- 0.01
  probbase[probbase == "C"] <- 0.005
  probbase[probbase == "C-"] <- 0.002
  probbase[probbase == "D+"] <- 0.001
  probbase[probbase == "D"] <- 5e-04
  probbase[probbase == "D-"] <- 1e-04
  probbase[probbase == "E"] <- 1e-05
  probbase[probbase == "N"] <- 0
  probbase[probbase == ""] <- 0
  #probbase[1, 1:13] <- rep(0, 13)
  
  # Specify the columns to sum
  probbase<-as.data.frame(probbase)
  cols_to_sum <- c("D_LOCATE.C.7", "D_LOGIST.C.7", "D_SERVICE.C.7", "D_PERCEPT.C.7", "D_COST.C.7")
  
  # Create a new column with row sums
  probbase$total <- rowSums(sapply(probbase[, cols_to_sum], as.numeric))
  probbase$total <- probbase$total / 5
  # Remove only the summed columns
  probbase <- probbase[, !names(probbase) %in% cols_to_sum]
  
  #Run VA algorithm
  output_new <- codeVA(data = new_data, data.type = data.type,
                       model = model, version = version,
                       HIV = HIV, Malaria = Malaria)
  
  #CSMF
  CSMFs_new <- as.matrix(getCSMF(output_new))
  assigned_causes_new <- data.frame(
    Cause = rownames(CSMFs_new),
    Value = as.numeric(CSMFs_new[, 1]),
    row.names = NULL
  )
  #Calc P(Ai|A)
  total<-rep(0,length(block))
  for (j in 1:length(block)){
    subtotal<- rep(0,length(assigned_causes_new[,2]))
    for(i in 1:length(assigned_causes_new[,2])){
      subtotal[i]<-assigned_causes_new[i,2]*as.numeric(probbase[min(block)+j-1,16+i])
    }
    total[j]<- sum(subtotal)
  }
  return(list(Original = missing, Imputed = total))
}

my_function(VA= RandomVA1[38,], block=11:20)
my_function(VA= RandomVA1[38,], block=43:45)
my_function(VA= RandomVA1[38,], block=21:22)
my_function(VA= RandomVA1[38,], block=18:20)
