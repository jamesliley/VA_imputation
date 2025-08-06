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

##' @name logistic
##' @description calculates the logistic function 
##' @param x input
##' @return 1/(1+exp(-x))
##' @export
logistic=function(x) {
  return(1/(1+exp(-x)))
}


