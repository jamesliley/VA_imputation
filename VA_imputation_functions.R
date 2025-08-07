##**********************************************************************
##  Calibration of symptom frequencies in VA using imputation       ####
##   proposal. 
##   
##  Functions: to be developed into an R package.
##  
##  Nathan Higgins, James Liley, Eilidh Cowan                                                        
##  July 2025                                                              
##**********************************************************************
##
## 

## Functions to be added to this document once working and commented.

##Note for users: The following packages must be loaded for all functions to work
library(openVA)
library(nbc4va)

##' @name logistic
##' @description calculates the logistic function 
##' @param x input
##' @return 1/(1+exp(-x))
##' @export
logistic=function(x) {
  return(1/(1+exp(-x)))
}



##' Remove and impute (e.g. 'recover') a block of answers for InterVA models
##' @name recoverAnswers
##' @description Takes a dataset of (single) verbal autopsies, and a set of questions.
##' Sets the answers to those questions to 'missing', and attempts to impute
##' their values.
##' @param VA VA dataset, in the format of (e.g.) RandomVA1
##' @param block Set of indices for questions, corresponding to rows in VA
##' @param method Method for attaining cause-of-death distribution from VA answers; default 'OpenVA'
##' @param data.type Format of the VA data; default 'WHO2012'
##' @param model VA algorithm used; at the moment supports only 'InterVA'
##' @param version Version of InterVA used; default '4.03'
##' @param HIV Set as 'h','m','l'
##' @param Malaria Set as 'h','m','l'
##' @return A list of two data frames, one called 'Original' and one called 'Imputed'
##' @export
##' @examples
##' ## Impute block of questions 11-20
##' data(RandomVA1)
##' out=recoverAnswers(RandomVA1,block=11:20,method="OpenVA",
##'                    data.type = "WHO2012",model = "InterVA",
##'                    version = "4.03", HIV = "h", Malaria = "h")
##'                   
##' head(out$Original)
##' head(out$Imputed)
recoverAnswers=function(VA,block,method="OpenVA",data.type = "WHO2012",
                        model = "InterVA", version = "4.03", HIV = "h", 
                        Malaria = "h") {
  library(openVA)
  library(nbc4va)
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
  return(list(Original=VA[,block],Imputed=total))
}









##' Remove and impute (e.g. 'recover') a block of answers for InSilicoVA models
##' @name recoverAnswers2
##' @description Takes a dataset of verbal autopsies, a probbase and a set of questions.
##' Sets the answers to those questions to 'missing', and attempts to impute
##' their values.
##' @param VA VA dataset, in the format of (e.g.) RandomVA1
##' @param block Set of indices for questions, corresponding to rows in VA
##' @param method Method for attaining cause-of-death distribution from VA answers; default 'OpenVA'
##' @param data.type Format of the VA data; default 'WHO2012'
##' @param model VA algorithm used; at the moment supports only 'InterVA'
##' @param CondProbNum A  customised probability matrix with 245 row symptomns on 60 column causes in the same order as InterVA4 specification. For example input see condprobnum.
##' @return A list of two data frames, one called 'Original' and one called 'Imputed'
##' @export
##' @examples
##' ## Impute block of questions 11-20
##' data(RandomVA1)
##' out=recoverAnswers2(RandomVA1,block=11:20,method="OpenVA",
##'                    data.type = "WHO2012",model = "InSilicoVA",
##'                    CondProbNum=condprobnum)
##'                   
##' head(out$Original)
##' head(out$Imputed)
recoverAnswers2=function(VA,block,CondProbNum,method="OpenVA",data.type = "WHO2012",
                        model = "InSilicoVA") {
  library(openVA)
  library(nbc4va)
  new_data<- VA
  missing<- VA[,block]
  new_data[, block]<- ""
  
  #Run VA algorithm
  output_new <- insilico(data = new_data, data.type = data.type,Nsim = 10000, auto.length = FALSE,
                         CondProbNum = CondProbNum)
  
  #CSMF
  indivprob<-output_new$indiv.prob
  imputed<- matrix(rep(0, length(indivprob[,1]) * length(block)), 
                     nrow = length(indivprob[,1]), ncol = length(block))
  
                   
  for(k in 1:length(indivprob[,1])){
    CSMFs_new<- as.vector(indivprob[k,])
    #Calc P(Ai|A)
    total<-rep(0,length(block))
    for(j in 1:length(block)){
      subtotal<- rep(0,length(CSMFs_new))
      for(i in 1:length(CSMFs_new)){
        subtotal[i]<-CSMFs_new[i]*as.numeric(CondProbNum[min(block)+j-2,i])
        }
      total[j]<- sum(subtotal)
      }
    imputed[k,]<- total
    }
  return(list(Original=VA[,block],Imputed=imputed))
}






