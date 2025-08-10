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








