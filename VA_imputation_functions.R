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
#Allow custom probbase for InterVA2
##**********************************************************************
##


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
##' @return InterVA output
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
##' probbase2[10,22]<-"A"
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
#Allow change in probbase element for interVA model by epsilon
##**********************************************************************
##


##' Run InterVA4 algorithm whilst changing element of custom probbase/letter probabilities by epsilon
##' @name InterVAepsilon
##' @description This function implements the algorithm in the InterVA4 software. 
##' Allows for customised probbase of the exact form of data(probbase) and also allows custom letter probabilities for the probbase. It specifically changes one element of the probbase by epsilon before running the algorithm
##' @param probBase Default null; uses data(probbase). Any custom probbase must be of exact same format as data(probbase).
##' @param letterProb A data frame with first coloumn letter and second coloumn the corresponding probabilities. If left null then will use the following default: data.frame(grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),value = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
##' @param probbase_element An element of the probbase you want to change by epsilon indexed as a list
##' @param epsilon A small value you want to add to the selected probbase element. If causes probability to go out of range (0,1) then reverse operation(addition or subtraction) will be made so within the range. Absolute value of epsilon should be <=0.1. Epsilon can be positive or negative.
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
##' @return InterVA output
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
##' probbase2[10,22]<-"A"
##' 
##' 
##' #Custom letter probabilities and select probbase2[10,25] to change by +0.4
##' 
##' sample.output1 <- InterVA2(SampleInput, HIV = "h", Malaria = "l", write=FALSE, 
##' replicates = FALSE, probBase=probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)), probbase_element=list(10,25), epsilon = 0.4)
##' 
##' 
##' to get causes of death with group code for further usage
##' sample.output2 <- InterVA2(SampleInput, HIV = "h", Malaria = "l", write=FALSE,
##' replicates = FALSE, groupcode = TRUE, probBase=probbase2, letterProb = NULL,probbase_element=list(10,25), epsilon = 0.4)
##' 
##' 
##' Using default probbase and letter probabilities
##' sample.output3 <- InterVA2(SampleInput, HIV = "h", Malaria = "l", write=FALSE,
##' replicates = FALSE, groupcode = TRUE, probBase=NULL, letterProb = NULL,probbase_element=list(10,25), epsilon = 0.4)


InterVAepsilon<-function (Input, HIV, Malaria, directory = NULL, filename = "VA_result", 
                    output = "classic", append = FALSE, groupcode = FALSE, replicates = FALSE,probBase=NULL, letterProb=NULL,probbase_element, epsilon,
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
  #NEW NEW SECTION- changing single element of probbase by epsilon
  probbase[probbase_element[[1]], probbase_element[[2]]]<- as.numeric(probbase[probbase_element[[1]], probbase_element[[2]]]) + epsilon
  if(probbase[probbase_element[[1]], probbase_element[[2]]] > 1){
    probbase[probbase_element[[1]], probbase_element[[2]]]<- as.numeric(probbase[probbase_element[[1]], probbase_element[[2]]]) - 2*epsilon
  }
  if(probbase[probbase_element[[1]], probbase_element[[2]]] < 0){
    probbase[probbase_element[[1]], probbase_element[[2]]]<- as.numeric(probbase[probbase_element[[1]], probbase_element[[2]]]) - 2*epsilon
  }
  #END NEW NEW SECTION
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
#Allow custom probbase/letter grade for InterVA using codeVA function
##**********************************************************************
##


##' Version of codeVA which provides InterVA4 analysis with custom probbase/letter probability options
##' @name codeVA2
##' @description This function implements the codeVA function but allows for customised probbase of the exact form of data(probbase) and also allows custom letter probabilities for the probbase when using InterVA.
##' @param probBase Default null; uses data(probbase). Any custom probbase must be of exact same format as data(probbase). Only use custom probbase for InterVA model; see insilico function for custom probbase for InSilicoVA model
##' @param letterProb A data frame with first coloumn letter and second coloumn the corresponding probabilities. Should only be specified if using InterVA model and want to change default. If left null then will use the following default: data.frame(grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),value = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
##' @param probbase_element An element of the probbase you want to change by epsilon indexed as a list
##' @param epsilon A small value you want to add to the selected probbase element. If causes probability to go out of range (0,1) then reverse operation(addition or subtraction) will be made so within the range. Absolute value of epsilon should be <=0.1. Epsilon can be positive or negative.
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
##' data(RandomVA3)
##' test <- RandomVA3[1:200, ]
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
##' probbase2[10,22]<-"A"
##' fit2 <- codeVA2(data = test, data.type = "customize", model = "InterVA",
##' data.train = train, causes.train = "cause", write=FALSE,
##' version = "4.02", HIV = "h", Malaria = "l", probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))
##' 
##' 
##' select probbase2[10,25] to change by +0.4
##' fit2.1 <- y(data = test, data.type = "customize", model = "InterVA",
##' data.train = train, causes.train = "cause", write=FALSE,
##' version = "4.02", HIV = "h", Malaria = "l", probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)), probbase_element = list(10,25), epsilon = 0.4)
##' 
##' 
##' #Default probbase and letter probabilities for InterVA
##' fit2 <- codeVA2(data = test, data.type = "customize", model = "InterVA",
##' data.train = train, causes.train = "cause", write=FALSE,
##' version = "4.02", HIV = "h", Malaria = "l")


codeVA2<-function (data, data.type = c("WHO2012", "WHO2016", "PHMRC", 
                              "customize")[2], data.train = NULL, causes.train = NULL, probBase=NULL, letterProb=NULL,probbase_element=NULL, epsilon=NULL,
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
      #NEW NEW SECTION
      if(!is.null(probbase_element) && !is.null(epsilon)){
      fit <- calibrateVA::InterVAepsilon(Input = data, HIV = HIV,probBase=probBase, letterProb = letterProb,
                               Malaria = Malaria, replicates = replicates,probbase_element=probbase_element, epsilon=epsilon, ...)
      }
      else{
      fit <- calibrateVA::InterVA2(Input = data, HIV = HIV,probBase=probBase, letterProb = letterProb,
                               Malaria = Malaria, replicates = replicates, ...)
      }
      #END NEW NEW SECTION
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
##' @param probbase_element An element of the probbase you want to change by epsilon indexed as a list
##' @param epsilon A small value you want to add to the selected probbase element. If causes probability to go out of range (0,1) then reverse operation(addition or subtraction) will be made so within the range. Absolute value of epsilon should be <=0.1. Epsilon can be positive or negative.
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
##' 
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
##' 
##' data("probbase")
##' probbase2<- probbase
##' probbase2[10,22]<-"A"
##' out=recoverAnswers(RandomVA1,block=11:20,method="OpenVA",
##' data.type = "WHO2012",model = "InterVA",
##' version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))
##' head(out$Original)
##' head(out$Imputed)
##' 
##' select probbase2[10,25] to change by +0.4
##' 
##' out=recoverAnswers(RandomVA1,block=11:20,method="OpenVA",
##' data.type = "WHO2012",model = "InterVA",
##' version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)), probbase_element =list(10,25) , epsilon = 0.4)
##' head(out$Original)
##' head(out$Imputed)


recoverAnswers= function(VA,block,method="OpenVA",data.type = "WHO2012",
                           model = "InterVA", version = "4.03", HIV = "h", 
                           Malaria = "h", probBase=NULL, letterProb=NULL, probbase_element=NULL, epsilon=NULL) {
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
  #NEW NEW SECTION- changing single element of probbase by epsilon
  if(!is.null(probbase_element) && !is.null(epsilon)){
    probbase[probbase_element[[1]], probbase_element[[2]]]<- as.numeric(probbase[probbase_element[[1]], probbase_element[[2]]]) + epsilon
    if(probbase[probbase_element[[1]], probbase_element[[2]]] > 1){
      probbase[probbase_element[[1]], probbase_element[[2]]]<- as.numeric(probbase[probbase_element[[1]], probbase_element[[2]]]) - 2*epsilon
    }
    if(probbase[probbase_element[[1]], probbase_element[[2]]] < 0){
      probbase[probbase_element[[1]], probbase_element[[2]]]<- as.numeric(probbase[probbase_element[[1]], probbase_element[[2]]]) - 2*epsilon
    }
  }
  #END NEW NEW SECTION
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
                                     HIV = HIV, Malaria = Malaria, probBase = probBase, letterProb = letterProb, 
                                     probbase_element = probbase_element, epsilon = epsilon)
  
  #CSMF
  #Essentially this function takes the raw probabilities from interVA and sets the answer to ‘undetermined’
  #with probability 1 if no probability exceeds 0.4. So you are usually getting this answer, whereas you 
  #actually just want the raw probabilities. This is why interVA.rule=FALSE
  CSMFs_new <- as.matrix(getCSMF(output_new, interVA.rule=FALSE))
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
##' data(condprobnum)
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
##' their values.
##' @param VA VA dataset, in the format of (e.g.) RandomVA1
##' @param block Set of indices for questions, corresponding to rows in VA
##' @param probBase Default null; uses data(probbase). Any custom probbase must be of exact same format as data(probbase). Only use custom probbase for InterVA model; see insilico function for custom probbase for InSilicoVA model
##' @param letterProb A data frame with first coloumn letter and second coloumn the corresponding probabilities. Should only be specified if using InterVA model and want to change default. If left null then will use the following default: data.frame(grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),value = c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02,0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0, 0))
##' @param probbase_element An element of the probbase you want to change by epsilon indexed as a list
##' @param epsilon A small value you want to add to the selected probbase element. If causes probability to go out of range (0,1) then reverse operation(addition or subtraction) will be made so within the range. Absolute value of epsilon should be <=0.1. Epsilon can be positive or negative.
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
##' probbase2[10,22]<-"A"
##' 
##' out=recoverAnswers3(RandomVA1,block=11:20,method="OpenVA",
##' data.type = "WHO2012",model = "InterVA",
##' version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)))
##' 
##' head(out$Original)
##' head(out$Imputed)
##' 
##' select probbase2[10,25] to change by +0.4
##' out=recoverAnswers3(RandomVA1,block=11:20,method="OpenVA",
##' data.type = "WHO2012",model = "InterVA",
##' version = "4.03", HIV = "h", Malaria = "h",probBase = probbase2, letterProb = data.frame(
##' grade = c("I", "A+", "A", "A-", "B+", "B", "B-", "B -", 
##' "C+", "C", "C-", "D+", "D", "D-", "E", "N", ""),
##' value = c(1, 1, 0.7, 0.2, 0.9, 1, 0.02, 0.4,
##' 0.01, 0.005, 0, 0, 0.2, 0, 0.00001, 0, 0)), probbase_element =list(10,25) , epsilon = 0.4)
##' 
##' head(out$Original)
##' head(out$Imputed)


recoverAnswers3=function(VA,block,method="OpenVA",data.type = "WHO2012",
                           model = "InterVA", version = "4.03", HIV = "h", 
                           Malaria = "h", probBase=NULL, letterProb=NULL, probbase_element=NULL,epsilon=NULL) {
  library(openVA)
  library(nbc4va)
  initial_out<-recoverAnswers(VA[1,],block,method="OpenVA",data.type = "WHO2012",
                              model = "InterVA", version = "4.03", HIV = "h", 
                              Malaria = "h", probBase=probBase, letterProb=letterProb, probbase_element = probbase_element,
                              epsilon = epsilon)
  Original<- initial_out$Original
  Imputed<- initial_out$Imputed
  for(i in 2:length(VA[,1])){
    out<-recoverAnswers(VA[i,],block,method="OpenVA",data.type = "WHO2012",
                        model = "InterVA", version = "4.03", HIV = "h", 
                        Malaria = "h", probBase=probBase, letterProb=letterProb, probbase_element = probbase_element,
                        epsilon = epsilon)
    Original<- rbind(Original, out$Original)
    Imputed<- rbind(Imputed, out$Imputed)
  }
  return(list(Original=Original,Imputed=Imputed))
}










##**********************************************************************
#evaluateProbbase function
##**********************************************************************
##


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






##**********************************************************************
#gradientProbbase function
##**********************************************************************
##



##' Partial Derivative of Imputation Accuracy with respect to a single element of the probbase
##' @name gradientProbbase
##' @description Takes a dataset of verbal autopsies, a set of multiple blocks of questions and a VA algorithm, a value for epsilon and a single element of the probbase.
##' Sets the answers to those questions to 'missing', and attempts to impute their values for each block in turn.
##' Then calculates overall accuracy of the probbase by averaging imputed values across all answers set to missing. 
##' Repeat this for a modified probbase by changing a single element by epsilon. Partial derivative with respect to this element of the probbase is the difference between accuracy divided by epsilon.
##' @param VA VA dataset, in the format of (e.g.) RandomVA1
##' @param blocks A list of blocks of questions. Each block is a set of indices for questions, corresponding to rows in VA
##' @param CondProbNum A customised probability matrix with 245 row symptomns on 60 column causes in the same order as InterVA4 specification. For example input see condprobnum.
##' @param epsilon A small value you want to add to the selected probbase element. If causes probability to go out of range (0,1) then reverse operation(addition or subtraction) will be made so within the range. Absolute value of epsilon should be <=0.1. Epsilon can be positive or negative.
##' @param probbase_element The element of the probbase you want to find the partial derivative with respect to indexed by a list
##' @param method Method for attaining cause-of-death distribution from VA answers; default 'OpenVA'
##' @param data.type Format of the VA data; default 'WHO2012'
##' @param model VA algorithm used; at the moment supports 'InterVA' and "InSilicoVA"
##' @param version Version of InterVA used; default '4.03'
##' @param HIV Set as 'h','m','l' for InterVA
##' @param Malaria Set as 'h','m','l' for InterVA
##' @return Partial Derivative of Imputation Accuracy with respect to a single element of the probbase
##' @export
##' @examples
##' ## Impute blocks of questions 18-20 and 21-24, use probbase condprobnum and change element (1,1) of probbase by 0.1
##' data("condprobnum")
##' data("RandomVA1")
##' out<- gradientProbbase(RandomVA1[1:100,],block=list(18:20, 21:24),method="OpenVA",data.type = "WHO2012",
##' model = "InSilicoVA",CondProbNum=condprobnum, probbase_element = list(1,1) , epsilon = 0.1)
##' 
##' 
##' ## Impute blocks of questions 18-20 and 21-24, use probbase2 and change element (10,25) of probbase by 0.4
##' data("probbase")
##' data("RandomVA1")
##' probbase2<- probbase
##' probbase2[10,22]<-"A"
##' out<- gradientProbbase(RandomVA1[1:100,],block=list(18:20, 21:24),method="OpenVA",data.type = "WHO2012",
##' model = "InterVA",probBase=probbase2, probbase_element = list(10,25) , epsilon = 0.4)
##' 
##' 
##' ## Impute blocks of questions 18-20 and 21-24, use probbase and change element (10,25) of probbase by 0.4
##' data("probbase")
##' data("RandomVA1")
##' out<- gradientProbbase(RandomVA1[1:100,],block=list(18:20, 21:24,50:55),method="OpenVA",data.type = "WHO2012",
##' model = "InterVA",probBase=probbase, probbase_element = list(10,25) , epsilon = -0.7)
##' 



gradientProbbase<- function(VA, blocks, epsilon, method = "OpenVA", probbase_element, probBase=NULL, letterProb=NULL,
                             data.type = "WHO2012", model,CondProbNum=NULL, version = NULL, HIV = NULL, 
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
                                  probBase = probBase,
                                  letterProb = letterProb,
                                  version = version,
                                  HIV = HIV,
                                  Malaria = Malaria)
  if(model == "InSilicoVA"){
    new_probbase<- CondProbNum
    new_probbase[probbase_element[[1]], probbase_element[[2]]]<- new_probbase[probbase_element[[1]],probbase_element[[2]]] + epsilon
    #New part- if out of range (0,1), need to minus epsilon then minus epsilon again to get reverse operation applied to orginal probbase
    if(new_probbase[probbase_element[[1]], probbase_element[[2]]] > 1){
      new_probbase[probbase_element[[1]], probbase_element[[2]]]<- as.numeric(new_probbase[probbase_element[[1]], probbase_element[[2]]]) - 2*epsilon
    }
    if(new_probbase[probbase_element[[1]], probbase_element[[2]]] < 0){
      new_probbase[probbase_element[[1]], probbase_element[[2]]]<- as.numeric(new_probbase[probbase_element[[1]], probbase_element[[2]]]) - 2*epsilon
    }
    accuracy_new <- evaluateProbbase(VA=VA,
                                     blocks=blocks,
                                     method = method,
                                     data.type = data.type,
                                     model=model,
                                     CondProbNum = new_probbase,
                                     version = version,
                                     HIV = HIV,
                                     Malaria = Malaria)
    }
  if(model == "InterVA"){
    accuracy_new <- evaluateProbbase(VA=VA,
                                     blocks=blocks,
                                     method = method,
                                     data.type = data.type,
                                     model=model,
                                     probBase = probBase,
                                     letterProb = letterProb,
                                     version = version,
                                     HIV = HIV,
                                     Malaria = Malaria, probbase_element = probbase_element, epsilon = epsilon)
    }
  gradient<- (accuracy_new - accuracy_original)/epsilon
  return(Gradient = gradient)
}


##**********************************************************************
#Simulation of VA data- first version
##**********************************************************************
##


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


