## Temporary, while package is not constructed
load("~/Research/VA/Git/VA_imputation/Package_temp/imputeVA/data/standard_simulation.RData")


##' Simulate a set of VA answers in WHO20112 format consistent with a given 
##'  probbase, prior on CoDs, conditional dependency structure, and dependency 
##'  matrices. 
##' 
##' @name simulateVA_WHO2012
##' @description This function generates simulated completed verbal autopsy 
##'  questionnaires in WHO2012 format. It does this assuming a given probbase is
##'  'correct'; that is, gives the correct probability of answering a particular
##'  question given a particular cause of death.
##'  
##' It also assumes a dependency structure. We recommend using the preset 
##'  dependency structure, or learning a dependency structure from an existing
##'  set of VA answers using the function learn_VA_dependency. 
##' 
##' Technically, answers b_1,b_2,...,b_n, are simulated
##'  answers according to the hierarchical model:
##'  
##'  c ~ Multinomial(pi)   # Cause of death
##'  z_1,z_2,... z_i|c ~ MultivariateNormal(0,V_{c})  # 'Latent' answers
##'  b_j|c=1{z_j < qnorm(p_{c,j})}
##' 
##' where pi is a prior on cause of death (CoD) frequencies, V_{c} is a 
##'  variance matrix associated with CoD c, and p_{c,j}, coming 
##'  from a probbase, describes the probability of answering question j as 'yes'
##'  given CoD c. 
##'  
##' The input to our simulation is a probbase, a prior, and a list of matrices
##'  V where V[[c]] is the covariance matrix for  CoD c. 
##' 
##' IDs are added at the start of each column.
##' 
##' @param n number of simulated VA answers
##' @param probbase A probbase; defaults to null in which case a standard probbase is used. Must be of the same form as data(probbase).
##' @param letterprob Correspondences between probbase letters and probabilities. Defaults to null, in which case a standard lookup is used.
##' @param prior A prior on cause-of-death types; defaults to null in which case a standardised prior is used.
##' @param V A list of lists of variance matrices as above; defaults to null in which case a standardised list is used.
##' @param agefreq Vector of frequencies of ages (length 6); need not be normalised.
##' @param sexfreq Vector of sex frequencies (length 2).
##' @param va_colnames Column names for VA output
##' @param seed Random seed; set to time if NULL.
##' @return List of two values: one (VA) is simulated VA output, the other is parameters used in simulation.
##' @export
##' @examples
##' 
##' source("~/Research/VA/Git/Package_temp/imputeVA/R/simulations.R")
# 
# data(RandomPhysician)
# code=RandomPhysician$code1
# va=RandomPhysician[,1:246]
# rparam=learn_parameters_from_va(va,code)
# 
# # Get standard simulation settings
# load("~/Research/VA/Git/Package_temp/imputeVA/data/standard_simulation.RData")
# 
# # Set cause-of-death specific covariance matrices.
# load("~/Research/VA/Git/Package_temp/imputeVA/data/cod_subcategories.RData") # Load correspondence between CoDs in RandomPhysician and probbase
# V_rp=list()
# for (i in 1:ncol(standard_simulation$probbase)) {
#   xi=which(cod_subcategories[,1]==colnames(standard_simulation$probbase)[i])
#   subcat=cod_subcategories[xi,2]
#   if (subcat!="NA") V_rp[[i]]=rparam$V[[subcat]] else V_rp[[i]]=rparam$V0
# }
# names(V_rp)=colnames(probbase)
# 
# # Simulate a VA database (10k entries) 
# vax=simulateVA_WHO2012(10000,seed=32212,V=V_rp)
# 
##' 
simulateVA_WHO2012=function(
    n,
    probbase=NULL, 
    letterprob=NULL,
    prior=NULL,
    V=NULL,
    agefreq=NULL,
    sexfreq=NULL,
    va_colnames=NULL,
    seed=NULL) {
  
  if (is.null(seed)) {
    set.seed(as.numeric(Sys.time()))
  } else {
    set.seed(seed)
  }
  
  # Establish defaults for missing inputs
  if (any(is.null(probbase),is.null(letterprob),is.null(prior),is.null(V))) {
    # data("standard_simulation")
    load("~/Research/VA/Git/VA_imputation/Package_temp/imputeVA/data/standard_simulation.RData")
  }
  if (is.null(probbase)) {
    probbase=standard_simulation$probbase
  }
  if (is.null(letterprob)) {
    letterprob=standard_simulation$letterprob
  }
  if (is.null(prior)) {
    prior=standard_simulation$prior
  }
  if (is.null(V)) {
    V=standard_simulation$V
  }
  if (is.null(agefreq)) {
    agefreq=standard_simulation$agefreq
  }
  if (is.null(sexfreq)) {
    sexfreq=standard_simulation$sexfreq
  }
  if (is.null(va_colnames)) {
    va_colnames=standard_simulation$va_colnames
  }
  
  ## ID
  id=paste0("d",1:n)
  
  ## Age
  age=matrix("",n,length(agefreq))
  agegroup=sample(1:length(agefreq),n,prob=agefreq,rep=T)
  age[cbind(1:n,agegroup)]="Y"
  colnames(age)=names(agegroup)
  
  ## Sex
  sex=matrix("",n,length(sexfreq))
  sexgroup=sample(1:length(sexfreq),n,prob=sexfreq,rep=T)
  sex[cbind(1:n,sexgroup)]="Y"
  colnames(sex)=names(sexgroup)
  
  # Cause of death (latent)
  clist=colnames(probbase)
  cod=sample(clist[17:76],n,rep=T,prob=prior[17:76])
  
  # Remainder
  dat=matrix("",n,dim(probbase)[1])
  ex=letterprob$value; names(ex)=letterprob$grade
  
  # Numeric probbase
  probbase_numeric=suppressWarnings(matrix(as.numeric(ex[probbase]),nrow(probbase),ncol(probbase)))
  colnames(probbase_numeric)=colnames(probbase)
  probbase_numeric[which(is.na(probbase_numeric))]=0

  for (i in 1:length(clist)) {
    w=which(cod==clist[i])
    if (length(w)>0) {
      rprob=probbase_numeric[,i] #ex[probbase[,i]]
      #rprob[which(is.na(rprob))]=0
      lq=matrix(rmnorm(length(w),mean=rep(0,length(rprob)),varcov=V[[clist[i]]]),length(w),length(rprob))
      comp=outer(rep(1,length(w)),qnorm(rprob))
      ans=lq<comp
      xans=matrix(c("","Y")[1+ans],nrow(lq),ncol(lq))
      dat[w,]=xans
    }
  }
  dat[,1+(1:(ncol(age)+ncol(sex)))]=cbind(age,sex)
  
  # Put in correct format
  dat=as.data.frame(dat)
  dat[[1]]=id
  for (i in 2:dim(dat)[2]) dat[,i]=factor(dat[[i]],levels=c("Y","","."))
  
  colnames(dat)=va_colnames
  
  return(list(VA=dat,
              parameters=list(
                probbase=probbase,
                letterprob=letterprob,
                probbase_numeric=probbase_numeric,
                prior=prior,
                agefreq=agefreq,
                sexfreq=sexfreq,
                V=V
                )))
}



##' Simulate a set of VA answers in WHO20112 forat consistent with a given 
##'  probbase, prior on CoDs, conditional dependency structure, and dependency 
##'  matrices. 
##' 
##' @name learn_parameters_from_va
##' @description This function takes a set of VA data, and optionally a set of 
##'  assigned disease codes, and learns a covariance structure to simulate it. 
##' 
##' VA answers are tested for correlation, and grouped into approximate 
##'  conditionally-independent blocks. 
##' 
##' @param n number of simulated VA answers
##' @param va A set of VA answers
##' @param code A vector containing (named) disease codes. 
##' @return A list of block-diagonal covariance matrices, one for each disease 
##'  code. 
##' @export
##' @examples
##'  See examples for simulate_va_WHO2012
learn_parameters_from_va=function(
    va,
    code=NULL) {
  
  # Convert VA to integers
  nq=dim(va)[2]
  cva=va; for (i in 1:ncol(cva)) cva[[i]]=as.character(cva[[i]])
  cva[is.na(cva)]="."
  cva[(cva=="")]="N"
  conv=c("Y"=1,"N"=0,"."=NA)
  nva=matrix(conv[as.matrix(cva)],nrow(va),ncol(va))
  
  # If we have codes, use them, otherwise use whole dataset
  if (is.null(code)) code=rep("Unknown",dim(va)[1])
  t0=table(code)
  ncode=names(t0)[which(t0>=10)]
  nc=length(ncode)
  
  # Pairwise tests for each code
  tests=array(1,dim=c(nq,nq,nc))
  mva=colMeans(nva,na.rm=T)
  xva=colSums(is.finite(nva))
  for (cc in 1:nc) {
    for (i in 2:nq) {
      for (j in 1:(i-1)) {
        vi=nva[,i]; vj=nva[,j]
        ni=xva[i]; nj=xva[j]
        res=1
        if (ni>8 & nj>8) {
          w=which(is.finite(vi + vj))
          if (length(w)>8 & var(vi[w]>0) & var(vj[w])>0) {
            res=suppressWarnings(cor.test(vi[w],vj[w])$p.value)
          }
        }
        
        tests[i,j,cc]=res
        tests[j,i,cc]=res
      }
    }
  }
  
  # Cluster for covariance
  pthresh=1-(1e-5) # Cutoff for same-group
  xmin=apply(tests,1:2,min)
  blox=cutree(hclust(as.dist(xmin)),h=pthresh)
  nb=length(unique(blox))
  bloxmat=outer(blox,blox,function(x,y) x==y)
  
  V_out=list()
  
  # General correlation matrix
  getcor=function(xva) {
    
    gen_cor0=suppressWarnings(cor(xva,use="pairwise.complete"))
    gen_cor1=diag(nrow(gen_cor0))
    ofreq=colMeans(xva,na.rm=TRUE); ofreq[which(is.na(ofreq))]=0
    for (i in 1:max(blox)) {
      w=which(blox==i)
      if (length(w)>1) {
        sub=gen_cor0[w,w]
        xy=outer(ofreq[w],rep(1,length(w))); yx=t(xy)
        psub=p2n(sub,xy,yx); psub1=psub
        cx=1
        ex=eigen(psub1)$values
        while(cx>1e-3 & (min(Re(ex))<1e-2 | max(Im(ex)>0))) { # shrink matrix to pos def
          cx=0.99*cx
          psub1=psub*(1-diag(length(w)))*cx + diag(length(w))
          ex=eigen(psub1)$values
        }
        if (cx<= 1e-3) psub1=diag(length(w))
        gen_cor1[w,w]=psub1
      }
    }
    return(gen_cor1) 
  }
  
  g0=getcor(nva)
  # image(gen_cor[order(blox),order(blox)])
  
  for (i in 1:length(unique(ncode))) {
    if (t0[ncode[i]]>10) {
      w=which(code==ncode[i])
      sub_cor=getcor(nva[w,])
      V_out[[i]]=sub_cor
    } else V_out[[i]]=g0
  }
  names(V_out)=ncode
  
  return(list(blocks=blox,V=V_out,V0=g0))
  
}



##' @name p2n
##' @description If (X,Y)~MVN(mean=c(0,0),variance=rbind(c(1,rho),c(rho,1))) 
##'  and cor(X<qnorm(a),Y<qnorm(b))=r, then given r, a, b, solves for rho 
##'  approximately. 
##'  
##' Uses a degree-ten polynomial approximation
##' 
##' @param xcor correlation of indicator variables
##' @param a threshold 1, as above
##' @param b threshold 2, as above
##' @return solution of equation as above. 
##' 
p2n=function(xcor,a,b) {
  
  #   b2n=function(a,b,r) {
  #     qa=qnorm(a); qb=qnorm(b)
  #     target=sqrt(a*b*(1-a)*(1-b))*r + a*b
  #     fr=function(rho) pmnorm(c(qa,qb),mean=c(0,0),varcov=rbind(c(1,rho),c(rho,1)))-target
  #     if (fr(0.99)*fr(-0.99) < 0) {
  #       ux=uniroot(fr,c(-0.99,0.99))
  #       return(ux$root)
  #     } else return(NA)
  #   }
  # 
  # nn=1000000
  # aa=1/(1+exp(-rnorm(nn,sd=2)))
  # bb=1/(1+exp(-rnorm(nn,sd=2)))
  # rr=pmin(1,pmax(-1,rnorm(nn,sd=0.33)))
  # dx=data.frame(aa,bb,rr,rho=rep(NA,nn))
  # for (i in 1:nn) {
  #   dx$rho[i]=b2n(aa[i],bb[i],rr[i])
  #   if ((i %% 10000)==0) print(i)
  # }
  # dx0=dx[which(!is.na(dx$rho)),]
  # deg=10
  # dx1=dx0
  # dx1$aa=qnorm(dx1$aa); dx1$bb=qnorm(dx1$bb);
  # nn=floor(dim(dx1)[1]/2); dx1a=dx1[1:nn,]; dx1b=dx1[(nn+1):dim(dx0)[1],]
  # l1=lm(rho~polym(aa,bb,rr,degree=deg,raw=TRUE),data=dx1a)
  # y2=predict(l1,newdata=dx1b)
  # par(mfrow=c(1,2))
  # plot(y2,dx1b$rho)
  # plot(l1$fitted.values,dx1a$rho)
  # 
  # l1c=l1$coefficients
  # terms=paste0(l1c[1])
  # for (i in 2:length(l1c)) {
  #   ni=names(l1c)[i]
  #   ix=unlist(strsplit(substring(ni,44,nchar(ni)),".",fixed=TRUE))
  #   terms[i]=paste0("(",l1c[i],")*(aa^",ix[1],
  #                   ")*(bb^",ix[2],")*(rr^",ix[3],")")
  # }
  # form=paste0(terms,collapse=" + ")
  # print(form)
  
  
  aa=qnorm(a); bb=qnorm(b); rr=xcor
  out=0.000155437383940755 + 
    (-2.46576649201711e-05)*(aa^1)*(bb^0)*(rr^0) + 
    (-0.000105114391043574)*(aa^2)*(bb^0)*(rr^0) + 
    (0.000246519496885033)*(aa^3)*(bb^0)*(rr^0) + 
    (-0.000303645755712534)*(aa^4)*(bb^0)*(rr^0) + 
    (-0.000220208915671345)*(aa^5)*(bb^0)*(rr^0) + 
    (0.000117885262910898)*(aa^6)*(bb^0)*(rr^0) + 
    (5.94448593330849e-05)*(aa^7)*(bb^0)*(rr^0) + 
    (-8.99732702757768e-06)*(aa^8)*(bb^0)*(rr^0) + 
    (-4.73458070993268e-06)*(aa^9)*(bb^0)*(rr^0) + 
    (-1.08042983611122e-07)*(aa^10)*(bb^0)*(rr^0) + 
    (0.000403221575197639)*(aa^0)*(bb^1)*(rr^0) + 
    (-0.000806773333795039)*(aa^1)*(bb^1)*(rr^0) + 
    (-6.28132736923196e-06)*(aa^2)*(bb^1)*(rr^0) + 
    (0.00526211844991422)*(aa^3)*(bb^1)*(rr^0) + 
    (-6.57676527509364e-05)*(aa^4)*(bb^1)*(rr^0) + 
    (-0.00185279557111216)*(aa^5)*(bb^1)*(rr^0) + 
    (2.20459092119484e-05)*(aa^6)*(bb^1)*(rr^0) + 
    (0.000162384041036247)*(aa^7)*(bb^1)*(rr^0) + 
    (-5.44090053762209e-06)*(aa^8)*(bb^1)*(rr^0) + 
    (-5.95056092072503e-06)*(aa^9)*(bb^1)*(rr^0) + 
    (-0.000483732553803503)*(aa^0)*(bb^2)*(rr^0) + 
    (9.82571676539191e-06)*(aa^1)*(bb^2)*(rr^0) + 
    (0.00148683021622206)*(aa^2)*(bb^2)*(rr^0) + 
    (0.000306078968441182)*(aa^3)*(bb^2)*(rr^0) + 
    (-0.000117173036558581)*(aa^4)*(bb^2)*(rr^0) + 
    (-0.000139974090745832)*(aa^5)*(bb^2)*(rr^0) + 
    (-5.38894487822983e-05)*(aa^6)*(bb^2)*(rr^0) + 
    (1.55848293409663e-05)*(aa^7)*(bb^2)*(rr^0) + 
    (5.56424540899932e-06)*(aa^8)*(bb^2)*(rr^0) + 
    (-0.00086568830029533)*(aa^0)*(bb^3)*(rr^0) + 
    (0.00339796199645789)*(aa^1)*(bb^3)*(rr^0) + 
    (5.69171984295738e-05)*(aa^2)*(bb^3)*(rr^0) + 
    (-0.0039925970172195)*(aa^3)*(bb^3)*(rr^0) + 
    (2.76051493384664e-05)*(aa^4)*(bb^3)*(rr^0) + 
    (0.000512197393567959)*(aa^5)*(bb^3)*(rr^0) + 
    (1.21522400078286e-05)*(aa^6)*(bb^3)*(rr^0) + 
    (-8.98142643196254e-06)*(aa^7)*(bb^3)*(rr^0) + 
    (6.01287086025693e-05)*(aa^0)*(bb^4)*(rr^0) + 
    (-0.000162043984634456)*(aa^1)*(bb^4)*(rr^0) + 
    (-0.000735865205370135)*(aa^2)*(bb^4)*(rr^0) + 
    (2.0617122289632e-05)*(aa^3)*(bb^4)*(rr^0) + 
    (6.52293845661529e-05)*(aa^4)*(bb^4)*(rr^0) + 
    (-8.52148996173494e-06)*(aa^5)*(bb^4)*(rr^0) + 
    (1.93311755640811e-06)*(aa^6)*(bb^4)*(rr^0) + 
    (0.000476139067331829)*(aa^0)*(bb^5)*(rr^0) + 
    (-0.000859765248486759)*(aa^1)*(bb^5)*(rr^0) + 
    (-3.43156862029608e-05)*(aa^2)*(bb^5)*(rr^0) + 
    (0.000526763926472771)*(aa^3)*(bb^5)*(rr^0) + 
    (-1.19574205341895e-05)*(aa^4)*(bb^5)*(rr^0) + 
    (-3.61498676649997e-05)*(aa^5)*(bb^5)*(rr^0) + 
    (6.4554859069407e-05)*(aa^0)*(bb^6)*(rr^0) + 
    (3.89923845095554e-05)*(aa^1)*(bb^6)*(rr^0) + 
    (0.000114828626781942)*(aa^2)*(bb^6)*(rr^0) + 
    (4.14072501796588e-06)*(aa^3)*(bb^6)*(rr^0) + 
    (-9.03067225671331e-06)*(aa^4)*(bb^6)*(rr^0) + 
    (-8.89014172576185e-05)*(aa^0)*(bb^7)*(rr^0) + 
    (-7.81260949566991e-06)*(aa^1)*(bb^7)*(rr^0) + 
    (3.76635288949844e-06)*(aa^2)*(bb^7)*(rr^0) + 
    (-9.26590434815379e-06)*(aa^3)*(bb^7)*(rr^0) + 
    (-1.55227785216941e-05)*(aa^0)*(bb^8)*(rr^0) + 
    (-4.99665817751357e-06)*(aa^1)*(bb^8)*(rr^0) + 
    (-3.84466460282001e-06)*(aa^2)*(bb^8)*(rr^0) + 
    (5.10437403171124e-06)*(aa^0)*(bb^9)*(rr^0) + 
    (1.82636677904366e-06)*(aa^1)*(bb^9)*(rr^0) + 
    (8.21706485051104e-07)*(aa^0)*(bb^10)*(rr^0) + 
    (1.61887813536931)*(aa^0)*(bb^0)*(rr^1) + 
    (-0.0048162947862526)*(aa^1)*(bb^0)*(rr^1) + 
    (0.153993853700454)*(aa^2)*(bb^0)*(rr^1) + 
    (0.0122079574374654)*(aa^3)*(bb^0)*(rr^1) + 
    (0.0719801776061687)*(aa^4)*(bb^0)*(rr^1) + 
    (-0.00745719221761)*(aa^5)*(bb^0)*(rr^1) + 
    (0.00166309230398423)*(aa^6)*(bb^0)*(rr^1) + 
    (0.00139280317247491)*(aa^7)*(bb^0)*(rr^1) + 
    (0.000239159173606566)*(aa^8)*(bb^0)*(rr^1) + 
    (-7.42029427973696e-05)*(aa^9)*(bb^0)*(rr^1) + 
    (-0.00430515341316747)*(aa^0)*(bb^1)*(rr^1) + 
    (0.0064065712250177)*(aa^1)*(bb^1)*(rr^1) + 
    (-0.00828927849304994)*(aa^2)*(bb^1)*(rr^1) + 
    (-0.00928730846928263)*(aa^3)*(bb^1)*(rr^1) + 
    (0.00422386864296986)*(aa^4)*(bb^1)*(rr^1) + 
    (0.0012699042590384)*(aa^5)*(bb^1)*(rr^1) + 
    (-0.000575689902814664)*(aa^6)*(bb^1)*(rr^1) + 
    (3.25824369917826e-05)*(aa^7)*(bb^1)*(rr^1) + 
    (3.19261098824362e-05)*(aa^8)*(bb^1)*(rr^1) + 
    (0.157120979461615)*(aa^0)*(bb^2)*(rr^1) + 
    (-0.000495260641099239)*(aa^1)*(bb^2)*(rr^1) + 
    (0.241983635633936)*(aa^2)*(bb^2)*(rr^1) + 
    (0.00226253199717725)*(aa^3)*(bb^2)*(rr^1) + 
    (-0.00299402179245382)*(aa^4)*(bb^2)*(rr^1) + 
    (-7.94102541101925e-05)*(aa^5)*(bb^2)*(rr^1) + 
    (-0.00140962878365239)*(aa^6)*(bb^2)*(rr^1) + 
    (0.00016097779594253)*(aa^7)*(bb^2)*(rr^1) + 
    (0.0110502315399881)*(aa^0)*(bb^3)*(rr^1) + 
    (-0.0032959626508287)*(aa^1)*(bb^3)*(rr^1) + 
    (0.00690360470526425)*(aa^2)*(bb^3)*(rr^1) + 
    (0.00595055842925371)*(aa^3)*(bb^3)*(rr^1) + 
    (-0.00145665057743269)*(aa^4)*(bb^3)*(rr^1) + 
    (-0.000470847094615287)*(aa^5)*(bb^3)*(rr^1) + 
    (-1.41781220263719e-05)*(aa^6)*(bb^3)*(rr^1) + 
    (0.0713608014611117)*(aa^0)*(bb^4)*(rr^1) + 
    (-0.00221240983280386)*(aa^1)*(bb^4)*(rr^1) + 
    (-0.00588601517816772)*(aa^2)*(bb^4)*(rr^1) + 
    (-0.00217397258386102)*(aa^3)*(bb^4)*(rr^1) + 
    (-0.000867740969689348)*(aa^4)*(bb^4)*(rr^1) + 
    (-0.000295689699628611)*(aa^5)*(bb^4)*(rr^1) + 
    (-0.00651237362645333)*(aa^0)*(bb^5)*(rr^1) + 
    (-0.000689537591585548)*(aa^1)*(bb^5)*(rr^1) + 
    (-0.00055453589095671)*(aa^2)*(bb^5)*(rr^1) + 
    (-0.000287897071640369)*(aa^3)*(bb^5)*(rr^1) + 
    (0.000220288996617952)*(aa^4)*(bb^5)*(rr^1) + 
    (0.00118175506280452)*(aa^0)*(bb^6)*(rr^1) + 
    (0.00155490676957539)*(aa^1)*(bb^6)*(rr^1) + 
    (-0.0011205086344504)*(aa^2)*(bb^6)*(rr^1) + 
    (0.000367631599265092)*(aa^3)*(bb^6)*(rr^1) + 
    (0.00107866200339054)*(aa^0)*(bb^7)*(rr^1) + 
    (8.4925309164347e-05)*(aa^1)*(bb^7)*(rr^1) + 
    (-3.90511827054194e-05)*(aa^2)*(bb^7)*(rr^1) + 
    (0.000347294785004155)*(aa^0)*(bb^8)*(rr^1) + 
    (-0.000130445000787997)*(aa^1)*(bb^8)*(rr^1) + 
    (-3.04605800237333e-05)*(aa^0)*(bb^9)*(rr^1) + 
    (-0.00361690512680858)*(aa^0)*(bb^0)*(rr^2) + 
    (0.012596682945989)*(aa^1)*(bb^0)*(rr^2) + 
    (-0.0283982288012754)*(aa^2)*(bb^0)*(rr^2) + 
    (-0.0117423851580687)*(aa^3)*(bb^0)*(rr^2) + 
    (0.0328919884274164)*(aa^4)*(bb^0)*(rr^2) + 
    (0.00460012466451092)*(aa^5)*(bb^0)*(rr^2) + 
    (-0.00530386839554475)*(aa^6)*(bb^0)*(rr^2) + 
    (-0.000433672834711797)*(aa^7)*(bb^0)*(rr^2) + 
    (-0.000190018303024016)*(aa^8)*(bb^0)*(rr^2) + 
    (-0.00900671242199665)*(aa^0)*(bb^1)*(rr^2) + 
    (0.533809056038243)*(aa^1)*(bb^1)*(rr^2) + 
    (-0.033984641022532)*(aa^2)*(bb^1)*(rr^2) + 
    (-1.96729774672806)*(aa^3)*(bb^1)*(rr^2) + 
    (0.0422817424900405)*(aa^4)*(bb^1)*(rr^2) + 
    (-0.0581882033258064)*(aa^5)*(bb^1)*(rr^2) + 
    (-0.0114680745047052)*(aa^6)*(bb^1)*(rr^2) + 
    (-0.00106185856781095)*(aa^7)*(bb^1)*(rr^2) + 
    (0.0130868348136995)*(aa^0)*(bb^2)*(rr^2) + 
    (-0.0140029402446516)*(aa^1)*(bb^2)*(rr^2) + 
    (0.00918598419449878)*(aa^2)*(bb^2)*(rr^2) + 
    (-0.0323324906940046)*(aa^3)*(bb^2)*(rr^2) + 
    (-0.0110770129554031)*(aa^4)*(bb^2)*(rr^2) + 
    (0.00255762364923108)*(aa^5)*(bb^2)*(rr^2) + 
    (0.00092501798437717)*(aa^6)*(bb^2)*(rr^2) + 
    (0.0269856999800542)*(aa^0)*(bb^3)*(rr^2) + 
    (-1.96071652608428)*(aa^1)*(bb^3)*(rr^2) + 
    (-0.0109531640013649)*(aa^2)*(bb^3)*(rr^2) + 
    (0.267618824751911)*(aa^3)*(bb^3)*(rr^2) + 
    (0.014730957658583)*(aa^4)*(bb^3)*(rr^2) + 
    (0.00311498145961001)*(aa^5)*(bb^3)*(rr^2) + 
    (-0.01100080530212)*(aa^0)*(bb^4)*(rr^2) + 
    (0.0147069082172989)*(aa^1)*(bb^4)*(rr^2) + 
    (-0.00985286381421208)*(aa^2)*(bb^4)*(rr^2) + 
    (0.00519884162684712)*(aa^3)*(bb^4)*(rr^2) + 
    (0.00016287162355905)*(aa^4)*(bb^4)*(rr^2) + 
    (-0.00879449210071235)*(aa^0)*(bb^5)*(rr^2) + 
    (-0.053053379022275)*(aa^1)*(bb^5)*(rr^2) + 
    (-0.00654948098607376)*(aa^2)*(bb^5)*(rr^2) + 
    (0.00383092967876382)*(aa^3)*(bb^5)*(rr^2) + 
    (0.0055881504684442)*(aa^0)*(bb^6)*(rr^2) + 
    (-0.00579254803860827)*(aa^1)*(bb^6)*(rr^2) + 
    (0.000964455235851981)*(aa^2)*(bb^6)*(rr^2) + 
    (0.000338783777715453)*(aa^0)*(bb^7)*(rr^2) + 
    (-0.00185202385799827)*(aa^1)*(bb^7)*(rr^2) + 
    (-0.000608780584592641)*(aa^0)*(bb^8)*(rr^2) + 
    (-0.937682175345588)*(aa^0)*(bb^0)*(rr^3) + 
    (0.0669220694508859)*(aa^1)*(bb^0)*(rr^3) + 
    (-0.965540711064828)*(aa^2)*(bb^0)*(rr^3) + 
    (-0.0738595007870735)*(aa^3)*(bb^0)*(rr^3) + 
    (1.50860891457194)*(aa^4)*(bb^0)*(rr^3) + 
    (-0.0028479061637264)*(aa^5)*(bb^0)*(rr^3) + 
    (0.187350514263808)*(aa^6)*(bb^0)*(rr^3) + 
    (0.0093875879646082)*(aa^7)*(bb^0)*(rr^3) + 
    (0.061374772632771)*(aa^0)*(bb^1)*(rr^3) + 
    (-0.0408679295529855)*(aa^1)*(bb^1)*(rr^3) + 
    (0.0736178312825493)*(aa^2)*(bb^1)*(rr^3) + 
    (-0.0610733230455386)*(aa^3)*(bb^1)*(rr^3) + 
    (0.0165141243012373)*(aa^4)*(bb^1)*(rr^3) + 
    (0.0312029382029666)*(aa^5)*(bb^1)*(rr^3) + 
    (0.00170219815203215)*(aa^6)*(bb^1)*(rr^3) + 
    (-1.02882536683349)*(aa^0)*(bb^2)*(rr^3) + 
    (-0.0323398689893839)*(aa^1)*(bb^2)*(rr^3) + 
    (10.7773429114047)*(aa^2)*(bb^2)*(rr^3) + 
    (-0.0524333513753839)*(aa^3)*(bb^2)*(rr^3) + 
    (-0.355724358840828)*(aa^4)*(bb^2)*(rr^3) + 
    (-0.00592881068681071)*(aa^5)*(bb^2)*(rr^3) + 
    (-0.114686644657706)*(aa^0)*(bb^3)*(rr^3) + 
    (0.0857380859093991)*(aa^1)*(bb^3)*(rr^3) + 
    (0.000515839849995309)*(aa^2)*(bb^3)*(rr^3) + 
    (0.0072798570442386)*(aa^3)*(bb^3)*(rr^3) + 
    (-0.0135700525713031)*(aa^4)*(bb^3)*(rr^3) + 
    (1.54491765218057)*(aa^0)*(bb^4)*(rr^3) + 
    (0.06359212967026)*(aa^1)*(bb^4)*(rr^3) + 
    (-0.378447206670631)*(aa^2)*(bb^4)*(rr^3) + 
    (-0.000660370153280046)*(aa^3)*(bb^4)*(rr^3) + 
    (0.0414728822362899)*(aa^0)*(bb^5)*(rr^3) + 
    (-0.00851447556170987)*(aa^1)*(bb^5)*(rr^3) + 
    (0.00707077170493188)*(aa^2)*(bb^5)*(rr^3) + 
    (0.181444103269122)*(aa^0)*(bb^6)*(rr^3) + 
    (-0.000205954598514967)*(aa^1)*(bb^6)*(rr^3) + 
    (0.00083219030775117)*(aa^0)*(bb^7)*(rr^3) + 
    (0.0582654389626858)*(aa^0)*(bb^0)*(rr^4) + 
    (-0.0937924175618522)*(aa^1)*(bb^0)*(rr^4) + 
    (0.239917993714369)*(aa^2)*(bb^0)*(rr^4) + 
    (0.033358176987036)*(aa^3)*(bb^0)*(rr^4) + 
    (-0.159778841146557)*(aa^4)*(bb^0)*(rr^4) + 
    (-0.0275601257236295)*(aa^5)*(bb^0)*(rr^4) + 
    (0.00461204110933534)*(aa^6)*(bb^0)*(rr^4) + 
    (0.0459916933336156)*(aa^0)*(bb^1)*(rr^4) + 
    (-3.27466321363161)*(aa^1)*(bb^1)*(rr^4) + 
    (0.23273589485233)*(aa^2)*(bb^1)*(rr^4) + 
    (-11.8154439268297)*(aa^3)*(bb^1)*(rr^4) + 
    (-0.0599375928127425)*(aa^4)*(bb^1)*(rr^4) + 
    (-0.399242095382815)*(aa^5)*(bb^1)*(rr^4) + 
    (-0.0473056551899642)*(aa^0)*(bb^2)*(rr^4) + 
    (0.184935659918376)*(aa^1)*(bb^2)*(rr^4) + 
    (-0.0245133563466198)*(aa^2)*(bb^2)*(rr^4) + 
    (0.0634934120714308)*(aa^3)*(bb^2)*(rr^4) + 
    (-0.0618255078730358)*(aa^4)*(bb^2)*(rr^4) + 
    (-0.154069675392178)*(aa^0)*(bb^3)*(rr^4) + 
    (-11.9030070610851)*(aa^1)*(bb^3)*(rr^4) + 
    (-0.0378787217408643)*(aa^2)*(bb^3)*(rr^4) + 
    (0.982142852417649)*(aa^3)*(bb^3)*(rr^4) + 
    (-0.0419680416367472)*(aa^0)*(bb^4)*(rr^4) + 
    (-0.12619191707568)*(aa^1)*(bb^4)*(rr^4) + 
    (0.0346959450252296)*(aa^2)*(bb^4)*(rr^4) + 
    (0.0239907213261728)*(aa^0)*(bb^5)*(rr^4) + 
    (-0.367170824174865)*(aa^1)*(bb^5)*(rr^4) + 
    (-0.00166275949867326)*(aa^0)*(bb^6)*(rr^4) + 
    (1.12250658265507)*(aa^0)*(bb^0)*(rr^5) + 
    (-0.292149593207336)*(aa^1)*(bb^0)*(rr^5) + 
    (4.17918779798415)*(aa^2)*(bb^0)*(rr^5) + 
    (0.200019505739248)*(aa^3)*(bb^0)*(rr^5) + 
    (3.05419545181827)*(aa^4)*(bb^0)*(rr^5) + 
    (-0.0328372138653153)*(aa^5)*(bb^0)*(rr^5) + 
    (-0.215752376447992)*(aa^0)*(bb^1)*(rr^5) + 
    (-0.0925634469452698)*(aa^1)*(bb^1)*(rr^5) + 
    (-0.385282500590179)*(aa^2)*(bb^1)*(rr^5) + 
    (0.473949987292759)*(aa^3)*(bb^1)*(rr^5) + 
    (-0.00660697450435117)*(aa^4)*(bb^1)*(rr^5) + 
    (4.37305660624933)*(aa^0)*(bb^2)*(rr^5) + 
    (0.124074443042935)*(aa^1)*(bb^2)*(rr^5) + 
    (14.6437219363445)*(aa^2)*(bb^2)*(rr^5) + 
    (0.181461168344565)*(aa^3)*(bb^2)*(rr^5) + 
    (0.162371437068663)*(aa^0)*(bb^3)*(rr^5) + 
    (-0.154842185102359)*(aa^1)*(bb^3)*(rr^5) + 
    (0.094551205409844)*(aa^2)*(bb^3)*(rr^5) + 
    (2.99040593892614)*(aa^0)*(bb^4)*(rr^5) + 
    (-0.0987882313437108)*(aa^1)*(bb^4)*(rr^5) + 
    (-0.0175579101743238)*(aa^0)*(bb^5)*(rr^5) + 
    (-0.351508422350644)*(aa^0)*(bb^0)*(rr^6) + 
    (0.192962795780464)*(aa^1)*(bb^0)*(rr^6) + 
    (-0.414251191812918)*(aa^2)*(bb^0)*(rr^6) + 
    (0.101239826687096)*(aa^3)*(bb^0)*(rr^6) + 
    (0.0334484029163793)*(aa^4)*(bb^0)*(rr^6) + 
    (-0.0153535952745056)*(aa^0)*(bb^1)*(rr^6) + 
    (-0.0999025400324779)*(aa^1)*(bb^1)*(rr^6) + 
    (-0.557618141908323)*(aa^2)*(bb^1)*(rr^6) + 
    (-3.68572920101527)*(aa^3)*(bb^1)*(rr^6) + 
    (0.171750116284471)*(aa^0)*(bb^2)*(rr^6) + 
    (-0.092083339347681)*(aa^1)*(bb^2)*(rr^6) + 
    (-0.280125790466323)*(aa^2)*(bb^2)*(rr^6) + 
    (0.164540221773481)*(aa^0)*(bb^3)*(rr^6) + 
    (-3.65087123835667)*(aa^1)*(bb^3)*(rr^6) + 
    (0.141807418974839)*(aa^0)*(bb^4)*(rr^6) + 
    (-1.97155070085856)*(aa^0)*(bb^0)*(rr^7) + 
    (0.393945158070186)*(aa^1)*(bb^0)*(rr^7) + 
    (-3.52450716061689)*(aa^2)*(bb^0)*(rr^7) + 
    (0.0217496195396178)*(aa^3)*(bb^0)*(rr^7) + 
    (0.354976852877355)*(aa^0)*(bb^1)*(rr^7) + 
    (0.259038072214588)*(aa^1)*(bb^1)*(rr^7) + 
    (0.0106036144298769)*(aa^2)*(bb^1)*(rr^7) + 
    (-3.61187170338637)*(aa^0)*(bb^2)*(rr^7) + 
    (0.0604099443115675)*(aa^1)*(bb^2)*(rr^7) + 
    (-0.0182728614046125)*(aa^0)*(bb^3)*(rr^7) + 
    (0.749060382385769)*(aa^0)*(bb^0)*(rr^8) + 
    (-0.100358395697089)*(aa^1)*(bb^0)*(rr^8) + 
    (0.0814626578615307)*(aa^2)*(bb^0)*(rr^8) + 
    (-0.0562266684329912)*(aa^0)*(bb^1)*(rr^8) + 
    (2.89107923004362)*(aa^1)*(bb^1)*(rr^8) + 
    (-0.204453638454221)*(aa^0)*(bb^2)*(rr^8) + 
    (1.35901693700101)*(aa^0)*(bb^0)*(rr^9) + 
    (-0.144076383478568)*(aa^1)*(bb^0)*(rr^9) + 
    (-0.217276749622328)*(aa^0)*(bb^1)*(rr^9) + 
    (-0.500150593430334)*(aa^0)*(bb^0)*(rr^10)
  
  out[which(is.na(out))]=0
  out=pmin(pmax(out,-1),1)
  return(out)
}


##' @name sanitise
##' @description 'Sanitises' a set of WHO2012 VAs by removing 'dontask' columns
##' @param va set of VA data in WHO2012 format, e.g. RandomVA1
##' @param probbase a probbase, e.g. data("probbase")
##' @return a 'sanitised' va
##' 
sanitise=function(va,probbase) {
  for (i in 14:246) {
    x=probbase[i,4:10]
    x=setdiff(x,c("","0"))
    if (length(x)>0) {
      j=match(toupper(x),toupper(colnames(va)))
      j=j[which(is.finite(j))]
      if (length(j)>0) {
        pb=rowSums(va[,j,drop=FALSE]=="Y")
        va[which(pb>0),i]=setdiff(unique(va[,i]),"Y")[1]
      }
    }
    y=probbase[i,12]
    y=setdiff(y,c("","0"))
    if (length(y)>0) {
      j=match(toupper(y),toupper(colnames(va)))
      pb=(va[,i,drop=FALSE]=="Y")
      va[which(pb>0),j]="Y"
    }
  }
  return(va)
}




##' @name fastVA
##' @description Implements an identical function to InterVA, but somewhat faster. 
##' Allows prior and probbase as inputs. 
##' 
##' Does not 'fix' VA elements; use the function ```sanitise``` for this. 
##' 
##' Does not name columns. 
##' 
##' @param va VA data in WHO2012 format
##' @param probbase_numeric a probbase, with numerical entries
##' @param prior a prior
##' @examples 
##' 
##' data("RandomVA1")
# va1=RandomVA1
# ix=InterVAx(va1,HIV="h",Malaria="h")
# iv=t(matrix(unlist(lapply(ix$VA,function(x) x$wholeprob)),63,dim(va1)[1]))
# w1=ix$VA[[1]]$wholeprob
# colnames(iv)=names(w1)
# ivsub=iv[,4:63]
# 
# ## Reproduce output
# 
# # Get standard probbase and causes
# data("probbase3"); probbase=probbase3
# data("causetext")
# 
# 
# # Probabilities for probbase; get numeric version
# lp=data.frame(
#   grade = c("I", "A+", "A", "A-", "B+", "B", "B-","B -", "C+", "C", "C-", "D+", 
#             "D", "D-", "E", "N", ""), 
#   value = c(1,0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.02, 0.01, 0.005, 0.002, 0.001,
#             5e-04, 1e-04, 1e-05, 0, 0))                                                                                                                    -17L))
# ex=lp$value; names(ex)=lp$grade
# rprob=matrix(ex[probbase],nrow(probbase),ncol(probbase))
# colnames(rprob)=colnames(probbase)
# rprob[which(is.na(rprob))]=0
# 
# # Fix 'don't ask if'/'implies' type questions
# xva=sanitise(va1,probbase)
# 
# # Make corrections to prior for high HIV/malaria
# prior=rprob[1,]
# prior[19] = 0.05 # high HIV
# prior[21] = 0.05 # high malaria
# prior[39] = 0.05 # high malaria
# 
# # Reconstruct posterior
# rpost=fastVA(xva,rprob,prior)
# 
# # Get names right
# colnames(rpost)=causetext[4:63,2]
# 
# # Essentially the same
# image(rpost - ivsub)
fastVA=function(va,probbase_numeric,prior) {
  yva=(va=="Y")
  lp=pmax(log(probbase_numeric),-20)
  prior=prior[17:76]
  rpost=exp(yva %*% lp)
  rpost=t(t(rpost[,17:76])*prior)
  mp=rowSums(rpost)
  rpost=rpost/mp
  return(rpost)
}



##' Auxiliary function
logC=function(x) log(pmin(pmax(x,1e-5),1-1e-5))

##' @name imputation_accuracy
##' @description Makes a fast assessment of imputation accuracy using InterVA 
##' @param va VA data
##' @param probbase_numeric Numerical probbase
##' @param prior prior probabilities
##' @param blocks block structure
imputation_accuracy=function(va,probbase_numeric,prior,blocks) {
  bid=unique(blocks) 
  pred=matrix(NA,nrow(va),ncol(va))
  actual=pred
  for (i in 1:length(bid)) {
    w=which(blocks==bid[i])
    va_in=(va[,w,drop=FALSE]=="Y") # True answers
    va_out=va
    va_out[,w]="" # Set to missing
    xpost=fastVA(va_out,probbase_numeric,prior) # Posterior
    va_predicted=(xpost %*% t(probbase_numeric[w,17:76,drop=FALSE]))
    
    # Fill in predictions
    pred[,w]=va_predicted
    actual[,w]=va_in
  }
  return(-sum(actual*logC(pred) + (1-actual)*logC(1-pred)))
}

