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



##**********************************************************************
#recoverAnswers function
##**********************************************************************
##

##' Remove and impute (e.g. 'recover') a block of answers for InterVA models
##' @name recoverAnswers
##' @description Takes a dataset of (single) verbal autopsies, and a set of questions.
##' Sets the answers to those questions to 'missing', and attempts to impute
##' their values.
##' @param VA VA dataset, in the format of (e.g.) RandomVA1
##' @param block Set of indices for questions, corresponding to rows in VA
##' @param probBase Default null; uses data(probbase). Any custom probbase must be of exact same format as data(probbase). Only use custom probbase for InterVA model; see insilico function for custom probbase for InSilicoVA model
##' @param letterProb A data frame with first coloumn letter and second coloumn the corresponding probabilities. Should only be specified if using InterVA model and want to change default. If left null then will use the following default: data.frame(grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),value = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
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
##' 
##' 
##' ## Impute block of questions 11-20 with custom probbase and letter probabilities
##' data("probbase")
##' probbase2<- probbase
##' probbase2[,19]<-"A"
##' out=recoverAnswers(RandomVA1,block=11:20,method="OpenVA",
##' data.type = "WHO2012",model = "InterVA",
##' version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))
##' head(out$Original)
##' head(out$Imputed)


recoverAnswers=function(VA,block,method="OpenVA",data.type = "WHO2012",
           model = "InterVA", version = "4.03", HIV = "h", 
           Malaria = "h", probBase=NULL, letterProb=NULL) {
  library(openVA)
  library(nbc4va)
  new_data<- VA
  missing<- VA[,block]
  new_data[, block]<- ""
  
  #NEW SECTION
  if(is.null(probBase)){
    data("probbase", envir = environment())
  }
  if(!is.null(probBase)){
    probbase<- probBase
  }
  if(is.null(letterProb)){
    letterProbs <- data.frame(Letter = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
                                         "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
                              Prob = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,
                                       0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
  }
  if(!is.null(letterProb)){
    letterProbs <- letterProb
  }
  
  #Letter probability conversion
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
  #probbase[1, 1:13] <- rep(0, 13)
  #END NEW SECTION
  
  # Specify the columns to sum
  probbase<-as.data.frame(probbase)
  cols_to_sum <- c("D_LOCATE.C.7", "D_LOGIST.C.7", "D_SERVICE.C.7", "D_PERCEPT.C.7", "D_COST.C.7")
  
  # Create a new column with row sums
  probbase$total <- rowSums(sapply(probbase[, cols_to_sum], as.numeric))
  probbase$total <- probbase$total / 5
  # Remove only the summed columns
  probbase <- probbase[, !names(probbase) %in% cols_to_sum]
  
  #Run VA algorithm
  output_new <- calibrateVA::codeVA2(data = new_data, data.type = data.type, 
                                     model = model, version = version, 
                                     HIV = HIV, Malaria = Malaria, probBase = probBase, letterProb = letterProb)
  
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











##**********************************************************************
#recoverAnswers2 function
##**********************************************************************
##

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




##**********************************************************************
#recoverAnswers3 function
##**********************************************************************
##


##' Remove and impute (e.g. 'recover') a block of answers for InterVA models
##' @name recoverAnswers3
##' @description Takes a dataset of verbal autopsies, and a set of questions.
##' Sets the answers to those questions to 'missing', and attempts to impute
##' their values. The probbase is fixed.
##' @param VA VA dataset, in the format of (e.g.) RandomVA1
##' @param block Set of indices for questions, corresponding to rows in VA
##' @param probBase Default null; uses data(probbase). Any custom probbase must be of exact same format as data(probbase). Only use custom probbase for InterVA model; see insilico function for custom probbase for InSilicoVA model
##' @param letterProb A data frame with first coloumn letter and second coloumn the corresponding probabilities. Should only be specified if using InterVA model and want to change default. If left null then will use the following default: data.frame(grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),value = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
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
##' out=recoverAnswers3(RandomVA1,block=11:20,method="OpenVA",
##'                    data.type = "WHO2012",model = "InterVA",
##'                    version = "4.03", HIV = "h", Malaria = "h")
##'                   
##' head(out$Original)
##' head(out$Imputed)
##' 
##' 
##' 
##' Impute block of questions 11-20 with custom probbase and letter probabilities
##' data("probbase")
##' probbase2<- probbase
##' probbase2[,19]<-"A"
##' 
##' out=y(RandomVA1,block=11:20,method="OpenVA",
##' data.type = "WHO2012",model = "InterVA",
##' version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))
##' 
##' head(out$Original)
##' head(out$Imputed)


recoverAnswers3=function(VA,block,method="OpenVA",data.type = "WHO2012",
                           model = "InterVA", version = "4.03", HIV = "h", 
                           Malaria = "h", probBase=NULL, letterProb=NULL) {
  library(openVA)
  library(nbc4va)
  initial_out<-recoverAnswers(VA[1,],block,method="OpenVA",data.type = "WHO2012",
                              model = "InterVA", version = "4.03", HIV = "h", 
                              Malaria = "h", probBase=probBase, letterProb=letterProb)
  Original<- initial_out$Original
  Imputed<- initial_out$Imputed
  for(i in 2:length(VA[,1])){
    out<-recoverAnswers(VA[i,],block,method="OpenVA",data.type = "WHO2012",
                        model = "InterVA", version = "4.03", HIV = "h", 
                        Malaria = "h", probBase=probBase, letterProb=letterProb)
    Original<- rbind(Original, out$Original)
    Imputed<- rbind(Imputed, out$Imputed)
  }
  return(list(Original=Original,Imputed=Imputed))
}








##**********************************************************************
#evaluateProbbase function
##**********************************************************************
##


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
##' probbase2[,19]<-"A"
##' evaluateProbbase(RandomVA1[1:100,],blocks=list(18:20, 21:24),method="OpenVA",data.type = "WHO2012",
##' model = "InterVA",version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))



                 
evaluateProbbase<-function(VA,blocks, method = "OpenVA",data.type = "WHO2012", model, CondProbNum=NULL, probBase=NULL, letterProb=NULL,
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
                          version = version, HIV = HIV, Malaria = Malaria, probBase=probBase, letterProb=letterProb)
    }
    if(model=="InSilicoVA"){
      out=recoverAnswers2(VA,block,method=method,
                          data.type = data.type,model = model,
                          CondProbNum=CondProbNum)
    }
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








##' Partial Derivative of Imputation Accuracy with respect to a single element of the probbase(for InSilico)
##' @name gradientProbbase
##' @description Takes a dataset of verbal autopsies, a set of multiple blocks of questions and a VA algorithm, a value for epsilon and a single element of the probbase.
##' Sets the answers to those questions to 'missing', and attempts to impute their values for each block in turn.
##' Then calculates overall accuracy of the probbase by averaging imputed values across all answers set to missing. 
##' Repeat this for a modified probbase by changing a single element by epsilon. Partial derivative with respect to this element of the probbase is the difference between accuracy.
##' @param VA VA dataset, in the format of (e.g.) RandomVA1
##' @param blocks A list of blocks of questions. Each block is a set of indices for questions, corresponding to rows in VA
##' @param CondProbNum A customised probability matrix with 245 row symptomns on 60 column causes in the same order as InterVA4 specification. For example input see condprobnum.
##' @param epsilon The amount to adjust the element of the probbase by
##' @param probbase_element The element of the probbase you want to find the partial derivative with respect to
##' @param method Method for attaining cause-of-death distribution from VA answers; default 'OpenVA'
##' @param data.type Format of the VA data; default 'WHO2012'
##' @param model VA algorithm used; at the moment supports 'InterVA' and "InSilicoVA
##' @param version Version of InterVA used; default '4.03'
##' @param HIV Set as 'h','m','l' for InterVA
##' @param Malaria Set as 'h','m','l' for InterVA
##' @return Accuracy of probbase
##' @export
##' @examples
##' ## Impute blocks of questions 18-20 and 21-24, use probbase condprobnum and change element (1,1) of probbase by 0.1
##' data("condprobnum")
##' data("RandomVA1")
##' out<- gradient_probbase(RandomVA1[1:100,],block=list(18:20, 21:24),method="OpenVA",data.type = "WHO2012",
##' model = "InSilicoVA",CondProbNum=condprobnum, probbase_element = list(1,1) , epsilon = 0.1)


gradient_probbase<- function(VA, blocks, epsilon, method = "OpenVA", probbase_element,
                             data.type = "WHO2012", model,CondProbNum, version = NULL, HIV = NULL, 
                             Malaria = NULL){
  #Calculate imputation accuracy with element of probbase set to original value.
  #Then change element of the probbase to original + epsilon
  #Calculate new imputation accuracy and record the change in accuracy to see sensitivity
  #Note: only work with InSilicoVA at the moment since only model which we can have custom probbase at the moment
  
  #Original probbase accuracy
  library(calibrateVA)
  library(openVA)
  library(nbc4va)
  accuracy_original<- evaluateProbbase(VA=VA,
                                  blocks=blocks,
                                  method = method,
                                  data.type = data.type,
                                  model=model,
                                  CondProbNum = CondProbNum,
                                  version = version,
                                  HIV = HIV,
                                  Malaria = Malaria)
  new_probbase<- CondProbNum
  new_probbase[probbase_element[[1]], probbase_element[[2]]]<- new_probbase[probbase_element[[1]],probbase_element[[2]]] + epsilon
  accuracy_new <- evaluateProbbase(VA=VA,
                                   blocks=blocks,
                                   method = method,
                                   data.type = data.type,
                                   model=model,
                                   CondProbNum = new_probbase,
                                   version = version,
                                   HIV = HIV,
                                   Malaria = Malaria)
  gradient<- (accuracy_new - accuracy_original)/epsilon
  return(Gradient = gradient)
}
