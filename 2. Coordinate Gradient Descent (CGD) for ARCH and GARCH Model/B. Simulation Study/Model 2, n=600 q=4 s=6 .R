setwd("C:/Users/Jaffri/Dropbox/Journal Papers drafts/2020/Simulation Studies-Rcode/B. Pure GARCH models/B. Simulation Study")
### LOAD R PACKAGES #######################################
#library(rpgm)
library(robustbase)
library(mvtnorm)
library(MASS)
library(fBasics)
library(matrixcalc)
###########################################################
### LOAD DEFINITIONS      ####################################################
SamSize<-600
maxREP<-1000
### LOAD DATASET GARCH(1,3) ####################################################
GARCH<-function(n){
  #A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
  library(mvtnorm)
  library(MASS)
  library(fBasics)
  library(matrixcalc)
  rCh<- 1
  sCh<- 3 
  n<-n+200  #Obs. below 201 will be deleted 
  Q<- 1 #Number of Variate
  
  ####### Generate r.v of white noise #########################################
  eta<-rmvnorm(n, mean=c(0), diag(1)) 
  eta<-t(eta)
  #############################################################################
  
  ##INTERCEPT
  paraGA0<-matrix(c(0.01),nrow=Q,byrow=TRUE)
  #### MA(p) representative
  paraGAa1<-matrix(c(0.3),nrow=Q,byrow=TRUE)
  #### AR(q) representative
  paraGAb1<-matrix(c(0.1),nrow=Q,byrow=TRUE)
  paraGAb2<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAb3<-matrix(c(0.55),nrow=Q,byrow=TRUE)

  H<-max(rCh,sCh)
  ht<-abs(t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1)))))  #HETEROS. OBS
  err<-t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1)))) #ERR. TERM
  Y <-t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1))))
  ZERO<-matrix(0,nrow=Q,ncol=(dim(eta)[2]-dim(Y)[2]))
  ht<-cbind(ht,ZERO) 
  err<-cbind(err,ZERO)
  Y<-cbind(Y,ZERO) 

GARCHr<-sapply( paste('paraGAa', 1:rCh, sep=''), get ,envir=sys.frame(sys.parent(0)))
GARCHs<-sapply( paste('paraGAb', 1:sCh, sep=''), get ,envir=sys.frame(sys.parent(0)))

    for (i in (H+1):((dim(eta)[2]))) {
      sumRg<-sumSg<-0
       for(j in 1:rCh){
        sumRg<-sumRg + (GARCHr[[j]]%*%((err[,(i-j)])^{2}))
       }
       for(k in 1:sCh){
        sumSg<-sumSg + (GARCHs[[k]]%*%ht[,(i-k)])
       }
     
    ht[,i]<-paraGA0+(sumRg)+(sumSg)
    err[,i]<-eta[,i]%*%sqrt(ht[,i])
    Y[,i]<-err[,i]
  }
Y<-Y[,c(-1:-200)]; #THIS TO REMOVE FIRST 200 OBSERVATIONS TO AVOID SET.SEED EFFECT
Y1<-as.matrix(Y)

coefGARCH<-cbind(paraGA0, paraGAa1, matrix(0,ncol=3),
  paraGAb1,paraGAb2,paraGAb3, matrix(0,ncol=3))

LAGNONZERO<-c(1,2,6,8)
LAGNONZERO2<-c(1,1,0,0,0,1,0,1,0,0,0)  #1: non-zero 0: zero
ALLDATA<-list(Y1,err,coefGARCH,rCh,sCh,LAGNONZERO,LAGNONZERO2,sum(LAGNONZERO2))
names(ALLDATA)<-c("Y1","err","coefGARCH","Gr","Gs","LAGNONZERO","LAGNONZERO2","COUNTNZ")
return(ALLDATA)
}
##############################################################################

### LOAD R FUNCTIONS #########################################################
Mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux)); 
  if(length(x)==length(ux[tab == max(tab)])){print("No mode")}else{ux[tab == max(tab)]}
}
#### BAYESIAN INFORMATION CRITERION (BIC) ###################################
BICFUNC<-function(SIGMASQ0,ERRSQ0,PARA,Ne,Cn){
  return(
    2*(sum(log(SIGMASQ0)+(ERRSQ0/SIGMASQ0)))+(sum(PARA != 0)*log(Ne)*Cn)
  )
}
#LIKELIHOOD IS MINIMIZED (NEGATIVE LIKELIHOOD)
#############################################################################

#### *** LIKELIHOOD FUNC ####################################################
F0<-function(PARA,LAMDDA1,AWEIGHT1,ERR2i){ 
  SUMB<-0
  L0<-0
  SIGMASQ0<-STEPB0(PARA,ERR2i)
  for(t in (h+1):length(ERR2i) ){
    L0 <- L0 +(
      (log(SIGMASQ0[t]) + (ERR2i[t]/(SIGMASQ0[t])))
    )
  }
  for(ii in 1:Dt){
    SUMB <- SUMB + (LAMDDA1*AWEIGHT1[ii]*abs(PARA[ii]))
  }
  L0<-L0+SUMB
  return(L0)	
}
###############################################################################

#### START STEP B. ##########################################################
STEPB0<-function(PARA,ERR2i){   ### RECURSIVE COMPUTATION OF SIGMA^2_{T}
if( (rC>0) && (sC==0) ){     ### ARCH(rC)
SIGMASQb<-numeric()
ERR2A<-SIGMASQA<-rep(0,h)
ERR2A<-c(ERR2A,ERR2i)    #MERGING INITIAL VALUES + PRE-COMPUTED SQUARED ERRORS

 	for( t in 1:length(ERR2i)){
    SUM1<-0
    for(ii in 1:rC){
    SUM1<- SUM1 + (PARA[ii+1]*ERR2A[t-ii+h])
	}
	SIGMASQb[t]<-SIGMASQA[t+h]<-PARA[1]+SUM1
  }
}else if( (rC>0) && (sC>0) ){ ### GARCH(rC,sC)
SIGMASQb<-numeric()
ERR2A<-SIGMASQA<-rep(0,h)
ERR2A<-c(ERR2A,ERR2i)    #MERGING INITIAL VALUES + PRE-COMPUTED SQUARED ERRORS
for( t in 1:length(ERR2i)){
  SUM1<-SUM2<-0
  for(ii in 1:rC){
    SUM1<- SUM1 + (PARA[ii+1]*ERR2A[t-ii+h])
  }
  for(jj in (rC+1):(Dt-1)){
    SUM2<- SUM2 + (PARA[jj+1]*SIGMASQA[t-jj+rC+h])
  }
  SIGMASQb[t]<-SIGMASQA[t+h]<-PARA[1]+SUM1+SUM2
  } #END FOR
}
return(SIGMASQb)
} #END STEPB0
#### END STEP B. ##############################################################

#### START STEP C. ############################################################
STEPC0<-function(PARA,SIGMASQa,ERR2i){   ### RECURSIVE COMPUTATION OF FIRST DERIVATIVE SIGMA^2_{T}
if( (rC>0) && (sC==0) ){     ### ARCH(rC)
	 ERR2A<-rep(ERR2i[1],h)
     ERR2A<-c(ERR2A,ERR2i)    #MERGING INITIAL VALUES + SQUARED ERRORS
     V<-DSIGMA<-list()
  ########################################
     for( t in 1:length(ERR2i)){
     	Ep<-numeric()
        for(ii in 1:rC){
           Ep[ii]<-ERR2A[t+h-ii]
        }
        V[[t]]<-as.matrix(c(1,Ep))
      }
  ########################################
     for( t in 1:length(ERR2i)){
     	DSIGMA[[t]]<-V[[t]]
     }
}else if( (rC>0) && (sC>0) ){ ### GARCH(rC,sC)
      ERR2A<-SIGMASQA<-rep(ERR2i[1],h) #INITIAL VALUES FOLLOWS SUGGESTION FROM FRANQ AND ZAKOIAN (2003)
      ERR2A<-c(ERR2A,ERR2i)    #MERGING INITIAL VALUES + SQUARED ERRORS
      SIGMASQA<-c(SIGMASQA,SIGMASQa) #MERGING INITIAL VALUES + SIGMA^{2}
      V<-DSIGMA<-DSIGMAt<-list()    #STORE VALUES IN LIST
      for(kk in 1:h){
       DSIGMAt[[kk]]<-matrix(c(0),nrow=Dt)
         }
  ########################################
     for( t in 1:length(ERR2i)){
      Ep<-Sigq<-numeric()
        for(ii in 1:rC){
           Ep[ii]<-ERR2A[t+h-ii]
        }
        for(jj in 1:sC){
           Sigq[jj]<-SIGMASQA[t+h-jj]
        }
        V[[t]]<-as.matrix(c(1,Ep,Sigq))  #MERGING VECTORS
      }
  ########################################
      for( t in 1:length(ERR2i)){
         SUM1<-0
        for(kk in (rC+1):(Dt-1)){ SUM1<- SUM1 + (PARA[kk+1]*DSIGMAt[[t-kk+rC+h]])} 
      DSIGMA[[t]]<-DSIGMAt[[t+h]]<-V[[t]] + SUM1
      }
  ########################################    
 }
LIST1<-list(V,DSIGMA)
return(LIST1) 
}
#### END STEP C. ##############################################################

#### START STEP D1. ##########################################################
STEPD1<-function(PARA,DSIGMAa,ERR2i){   ### RECURSIVE COMPUTATION OF SECOND DERIVATIVE SIGMA^2_{T}  ### FOR GARCH MODEL ONLY
DBETA<-D2SIGMAt<-DSIGMAt2<-D2SIGMAt2<-list()
ZEROMATD<-matrix(c(0),ncol=Dt,nrow=Dt)
##########################################################
      for(kk in 1:sC){
        DBETA[[kk]]<-matrix(c(0),ncol=Dt)
        DBETA[[kk]][1,(kk+rC+1)]<-1
      }

      for(kk in 1:h){
       DSIGMAt2[[kk]]<-matrix(c(0),nrow=Dt)   #INITIAL VALUE OF FIRST DERIVATIVE
       D2SIGMAt2[[kk]]<-matrix(c(0),ncol=Dt,nrow=Dt) #INITIAL VALUE OF THE SECOND DERIVATIVE
         }
DSIGMAt2<-c(DSIGMAt2,DSIGMAa)    #DSIGMA FROM STEP C.
##########################################################

for( t in 1:length(ERR2i)){
   A<-ZEROMATD
   for(i in 1:sC){
    A[(i+1+rC),]<-DSIGMAt2[[t+h-i]]
   }

   B<-ZEROMATD
   for(j in 1:(sC)){
     B<-B+(DSIGMAt2[[t+h-j]]%*%DBETA[[j]])
   }

   C<-ZEROMATD
   for(k in (rC+1):(Dt-1)){
     C<- C+(PARA[k+1]*D2SIGMAt2[[t+h-k+rC]])
   }

D2SIGMAt[[t]]<-D2SIGMAt2[[t+h]] <- (A+B+C)
    }
return(D2SIGMAt)
}
#### END STEP C(1). ##############################################################

#### START STEP D. ############################################################
STEPD0<-function(PARA,SIGMASQa,DSIGMAa,ERR2i){  ## COMPUTE THE SCORE VECTOR (Sn), TRUE HESSIAN (HEn), 
#INFOMATION MATRIX (IMn) AND OUTER PRODUCT MATRIX (OMn)
if( (rC>0) && (sC==0) ){     ### ARCH(rC)
Sn<-HEn<-IMn<-OMn<-0 
	for(t in (h+1):length(ERR2i) ){
      Dlike<-((1/SIGMASQa[[t]])- (ERR2i[t]/(SIGMASQa[[t]]^{2})))*DSIGMAa[[t]]
##SCORE MATRIX    
      Sn<-Sn+(Dlike)

##FULL HESSIAN
      HEn<-HEn+(
        ( 
          (((2*ERR2i[t])/(SIGMASQa[[t]]^{3})) - (1/(SIGMASQa[[t]])^2))
               *(DSIGMAa[[t]]%*%t(DSIGMAa[[t]]))
            )
        )

##INFO MATRIX 
      IMn<-IMn+(
               (((1/(SIGMASQa[[t]]^2)) ))*(DSIGMAa[[t]]%*%t(DSIGMAa[[t]])) 
             )

##OUTER PRODUCT MATRIX
      OMn<-OMn+(
        Dlike%*%t(Dlike)
      )
	}
}else if( (rC>0) && (sC>0) ){ ### GARCH(rC,sC)
D2SIGMAa<-STEPD1(PARA,DSIGMAa,ERR2)
Sn<-HEn<-IMn<-OMn<-HEn2<-0 
  for(t in (h+1):length(ERR2i) ){
      Dlike<-((1/SIGMASQa[[t]])- (ERR2i[t]/(SIGMASQa[[t]]^{2})))*DSIGMAa[[t]]
##SCORE MATRIX    
      Sn<-Sn+(Dlike)

##FULL HESSIAN 1
      HEn<-HEn+(
        (
          (  (1/(SIGMASQa[[t]])) - ((ERR2i[t])/(SIGMASQa[[t]]^2))  )*D2SIGMAa[[t]]

        )
          +
        ( 
          (((2*ERR2i[t])/(SIGMASQa[[t]]^{3})) - (1/(SIGMASQa[[t]])^2))
               *(DSIGMAa[[t]]%*%t(DSIGMAa[[t]]))
            )
        )

##INFO MATRIX 
      IMn<-IMn+(
               (((1/(SIGMASQa[[t]]^2)) ))*(DSIGMAa[[t]]%*%t(DSIGMAa[[t]])) 
             )

##OUTER PRODUCT MATRIX
      OMn<-OMn+(
        Dlike%*%t(Dlike)
      )

  }
}
LIST2<-list(Sn,HEn,IMn,OMn)
return(LIST2)
 }
#### END STEP D. ##############################################################

#### STEP E. ARMIJO RULE COMP. ################################################
STEPE0<-function(PARA,LAMBDA1,AWEIGHT1,L,dBi,Qq,ERR2i){  ## DETERMINATION OF \tau^{K-1}_s, THE STEP SIZE
ArM0<- 1 ##SET INITIAL STEP SIZE
ITERArmijo<-TRUE
Alp0<-0.1      #\tilde{\alpha} # default is 0.5 or 0.00001
Del0<-0.5     #\tilde{\delta} #default is 0.8 or 0.5
Gam0<-0        #(DEFAULT AS IN TSENG & YUN (2009), pg 414)
D0<-t(dBi)%*%(Qq%*%dBi)

LiKE1<-F0(PARA,LAMBDA1,AWEIGHT1,ERR2i)
while(ITERArmijo==TRUE){
if( (PARA[L]+(ArM0*(dBi[L]))) < 0 ){
  ArM0<-(Del0*ArM0)
  if( (ArM0 <= 1e-8 ) ){ITERArmijo<-FALSE}
  }else{
LiKE2<-F0((PARA+(ArM0*(dBi))),LAMBDA1,AWEIGHT1,ERR2i)
DIRECTION<-Alp0*ArM0*(Gam0-1)*(D0)
############ COMPUTE ARMIJO RULE #####################
   if( ( LiKE2 > (LiKE1+DIRECTION)) && (ArM0 > 1e-8 )){
    ArM0<-(Del0*ArM0)
     }else{ITERArmijo<-FALSE}
   }
}
return(ArM0)
}
#### END STEP E. ##############################################################

GARCHCGD<-function(LAMBDAi,AWEIGHTi,PARAs1CGDi){
##**** BEGIN: CGD ALGORITHM, FULL QMLE ########################################
PARAs1CGD<-PARAs1CGDi
############## (A) GCD WITH LAMBDA=0 ##########################################
ITERFULL<-TRUE
CYCLEGCD<-0

while(ITERFULL==TRUE){
BEFORE<-PARAs1CGD

for(L in 1:Dt){ ###################### CYCLE FORM 1 to D
####################### STEP B & C
SIGMASQ<-STEPB0(PARAs1CGD,ERR2)
STEPCresult<-STEPC0(PARAs1CGD,SIGMASQ,ERR2)
DSIGMA<-STEPCresult[[2]]
####################### STEP D
STEPDresult<-STEPD0(PARAs1CGD,SIGMASQ,DSIGMA,ERR2)
Sn<- STEPDresult[[1]]
HEn<-STEPDresult[[2]]
IMn<-STEPDresult[[3]]
OMn<-STEPDresult[[4]]/2
### SET THE SECOND DERIVATIVE OF LOG-LIKELIHOOD = Information Matrix
Qq<- IMn
##################################
dB<-rep(0,Dt)
  fA<-as.matrix(-Sn[L]-(Qq[L,L]*(-PARAs1CGD[L])))
    if(fA < (LAMBDAi*AWEIGHTi[L])){
      dB[L]<- -PARAs1CGD[L]
    }else{ 
      dB[L]<- (-(LAMBDAi*AWEIGHTi[L])-Sn[L])/Qq[L,L]
    }
####################### STEP D
####################### UPDATE PARAMETER
if(norm(as.matrix(dB),"1")!=0){
ArM<-STEPE0(PARAs1CGD,LAMBDAi,AWEIGHTi,L,dB,Qq,ERR2)
PARAs1CGD[L]<-max((PARAs1CGD[L]+(ArM*(dB[L]))),0)
  }
}##################################### END CYCLE
if( (norm(matrix(BEFORE-PARAs1CGD),type = "1") < 0.0001) || (CYCLEGCD==10000) ){
  ITERFULL<-FALSE
}
CYCLEGCD<-CYCLEGCD+1
  }
return(PARAs1CGD)
}
###############################################################################


### START : SIMULATION ########################################################
start_time <- Sys.time()
for(REPLICATE in 1:maxREP){
cat("REPLICATION=", REPLICATE,"\n")
rC<-4
sC<-6
Dt<-(1+rC+sC)
### STEP 1: generate dataset and other values. #################################
GARCH1<-GARCH(SamSize)
ERR<-as.vector(GARCH1$Y1)
ERR2<-as.vector(ERR)^2                 #COMPUTE SQUARE ERRORS
h<-max(rC,sC)
if(REPLICATE==1){
### COLLECT SERIES #############################################################
listSERIES1<-listERR2<-list()
#### TRUE VALUES #########
  TRUEPARA<-GARCH1$coefGARCH
  NONZERO<-GARCH1$LAGNONZERO
  CTRUENONZERO<-GARCH1$COUNTNZ
  propoTrue<-GARCH1$LAGNONZERO2
#### COUNT ESTIMATED NONZERO VARIABLES
  cNonZerofull<-cNonZeroadaplasso<-integer()
  cBIASNonZerofull<-cBIASNonZeroadaplasso<-integer()
  cEVARNonZerofull<-cEVARNonZeroadaplasso<-integer()
#### CHECK TRUE OR FALSE MODEL
  pickAllSigfull<-pickAllSigadaplasso<-0
  isTRUEinfull<-isTRUEinadaplasso<-0
#### VARIABLES FOR ESTIMATED PARAMETERS
  colESTfull<-colESTadaplasso<-list()
  colBIASfull<-colBIASadaplasso<-list()
  colEVARfull<-colEVARadaplasso<-list()
  cpropoESTfull<-cpropoESTadaplasso<-list()
               }
listSERIES1[[REPLICATE]]<-ERR
listERR2[[REPLICATE]]<-ERR2
### BEGIN : FULL CONSTRAINT QMLE ############################################
PARAinit<-c(mean(ERR2),rep(0, (Dt-1)))
LAMBDAa<-0
AWEIGHTa<-c(0,rep(1,rC),rep(1,sC))
PARAFULL<-GARCHCGD(LAMBDAa,AWEIGHTa,PARAinit)
### COLLECT RESULT: FULL QMLE
colESTfull[[REPLICATE]]<-PARAFULL
colBIASfull[[REPLICATE]]<-(PARAFULL-TRUEPARA)
colEVARfull[[REPLICATE]]<-(PARAFULL-TRUEPARA)^2
pickAllSigfull[REPLICATE]<-all(NONZERO%in%which(PARAFULL!=0))
isTRUEinfull[REPLICATE]<-isTRUE(all.equal(NONZERO,which(PARAFULL!=0)))
cNonZerofull[REPLICATE]<-sum(as.integer(PARAFULL!=0))
cBIASNonZerofull[REPLICATE]<- CTRUENONZERO-sum(as.integer(PARAFULL!=0))
cEVARNonZerofull[REPLICATE]<- (CTRUENONZERO-sum(as.integer(PARAFULL!=0)))^2
cpropoESTfull[[REPLICATE]]<-as.integer(as.integer(PARAFULL!=0))
### END : FULL CONSTRAINT QMLE ##############################################

#########** BEGIN ADAPTIVE LASSO ############################################
PARAWEIGHT<-PARAFULL
PARAWEIGHT[PARAWEIGHT==0]<-1e-20
BICadaplasso<-numeric()
PARAadlasso<-c(mean(ERR2),rep(0, (Dt-1)))
cPARAadaplasso<-list()
AWEIGHTb<-c(0,(PARAWEIGHT[2:(Dt)]^(-1)))
LAMBDAlist<-c(seq(20, 1, length.out= 20), seq(0.9, 0, length.out= 20))

for(KK in 1:length(LAMBDAlist)){
LAMBDAb<-LAMBDAlist[KK]
if(KK==1){
PARAinit2<-PARAadlasso<-cPARAadaplasso[[KK]]<-GARCHCGD(LAMBDAb,AWEIGHTb,PARAinit)
}else{
PARAinit2<-PARAadlasso<-cPARAadaplasso[[KK]]<-GARCHCGD(LAMBDAb,AWEIGHTb,PARAinit2)	
}
### COMPUTE BIC ################################################################
BICadaplasso[KK]<-BICFUNC(STEPB0(cPARAadaplasso[[KK]],ERR2),ERR2,cPARAadaplasso[[KK]],length(ERR2),5)
### END STEP 3B ################################################################
}
RePARAadaplassoBIC<-cPARAadaplasso[[which(BICadaplasso==min(BICadaplasso))[1]]]
#### COLLECT RESULT: ADAPTIVE LASSO
colESTadaplasso[[REPLICATE]]<-RePARAadaplassoBIC
colBIASadaplasso[[REPLICATE]]<-(RePARAadaplassoBIC-TRUEPARA)
colEVARadaplasso[[REPLICATE]]<-(RePARAadaplassoBIC-TRUEPARA)^2
pickAllSigadaplasso[REPLICATE]<-all(NONZERO%in%which(RePARAadaplassoBIC!=0))
isTRUEinadaplasso[REPLICATE]<-isTRUE(all.equal(NONZERO,which(RePARAadaplassoBIC!=0)))
cNonZeroadaplasso[REPLICATE]<-sum(as.integer(RePARAadaplassoBIC!=0))
cBIASNonZeroadaplasso[REPLICATE]<- CTRUENONZERO-sum(as.integer(RePARAadaplassoBIC!=0))
cEVARNonZeroadaplasso[REPLICATE]<- (CTRUENONZERO-sum(as.integer(RePARAadaplassoBIC!=0)))^2
cpropoESTadaplasso[[REPLICATE]]<-as.integer(as.integer(RePARAadaplassoBIC!=0))
### END : ADAPTIVE LASSO  ##############################################
}
end_time <- Sys.time()

#### COLLECT RESULTS ########################################################################
for(REPLICATE in 1:maxREP){
if(REPLICATE==1){
### RESULT: FULL QMLE
colPARAmeanQMLE<-cbind(colESTfull[[REPLICATE]])
colPARAbiasQMLE<-cbind(t(colBIASfull[[REPLICATE]]))
colPARAVARQMLE<-cbind(t(colEVARfull[[REPLICATE]]))
colPropofull<-cbind(cpropoESTfull[[REPLICATE]])
### RESULT: ADAPTIVE LASSO QMLE
colPARAmeanadapLASSO<-cbind(colESTadaplasso[[REPLICATE]])
colPARAbiasadapLASSO<-cbind(t(colBIASadaplasso[[REPLICATE]]))
colPARAVARadapLASSO<-cbind(t(colEVARadaplasso[[REPLICATE]]))
colPropoadapLASSO<-cbind(cpropoESTadaplasso[[REPLICATE]])
 }else{
### RESULT: FULL QMLE
colPARAmeanQMLE<-cbind(colPARAmeanQMLE,as.matrix(colESTfull[[REPLICATE]]))
colPARAbiasQMLE<-cbind(colPARAbiasQMLE,t(colBIASfull[[REPLICATE]]))
colPARAVARQMLE<-cbind(colPARAVARQMLE,t(colEVARfull[[REPLICATE]]))
colPropofull<-cbind(colPropofull,as.matrix(cpropoESTfull[[REPLICATE]]))
### RESULT: ADAPTIVE LASSO QMLE
colPARAmeanadapLASSO<-cbind(colPARAmeanadapLASSO,as.matrix(colESTadaplasso[[REPLICATE]]))
colPARAbiasadapLASSO<-cbind(colPARAbiasadapLASSO,t(colBIASadaplasso[[REPLICATE]]))
colPARAVARadapLASSO<-cbind(colPARAVARadapLASSO,t(colEVARadaplasso[[REPLICATE]]))
colPropoadapLASSO<-cbind(colPropoadapLASSO,as.matrix(cpropoESTadaplasso[[REPLICATE]]))
 }
}

### RESULT: FULL QMLE ##################################################
## COUNT NON-ZERO FOR FULL QMLE
min(cNonZerofull)                #MINIMUM
max(cNonZerofull)                #MAXIMUM
round(mean(cNonZerofull))        #MEAN (rounded)
Mode(cNonZerofull)               #MODE
sum(cBIASNonZerofull)/REPLICATE  #BIAS
sqrt(sum(cEVARNonZerofull)/REPLICATE)  #ESD

table(pickAllSigfull)
table(isTRUEinfull)

PARAmeanQMLE<-rowSums(colPARAmeanQMLE)/REPLICATE; 
rowMins(colPARAmeanQMLE)
rowMedians(colPARAmeanQMLE)
rowMaxs(colPARAmeanQMLE)
PARAbiasQMLE<-rowSums(colPARAbiasQMLE)/REPLICATE; 
PARAesdQMLE<-sqrt(rowSums(colPARAVARQMLE)/REPLICATE); 
PARApropoQMLE<-rowSums(colPropofull)/REPLICATE; 
########################################################################

### RESULT: ADAPTIVE LASSO #############################################
## COUNT NON-ZERO FOR ADAPTIVE LASSO
min(cNonZeroadaplasso)                #MINIMUM
max(cNonZeroadaplasso)                #MAXIMUM
round(mean(cNonZeroadaplasso))        #MEAN (rounded)
Mode(cNonZeroadaplasso)               #MODE
sum(cBIASNonZeroadaplasso)/REPLICATE  #BIAS
sqrt(sum(cEVARNonZeroadaplasso)/REPLICATE)  #ESD

table(pickAllSigadaplasso)
table(isTRUEinadaplasso)

PARAmeanadapLASSO<-rowSums(colPARAmeanadapLASSO)/REPLICATE; 
rowMins(colPARAmeanadapLASSO)
rowMedians(colPARAmeanadapLASSO)
rowMaxs(colPARAmeanadapLASSO)
PARAbiasadapLASSO<-rowSums(colPARAbiasadapLASSO)/REPLICATE; 
PARAesdadapLASSO<-sqrt(rowSums(colPARAVARadapLASSO)/REPLICATE); 
PARApropadapLASSO<-rowSums(colPropoadapLASSO)/REPLICATE;
##########################################################################

  simuTIME<- end_time - start_time
  TEXTtime<-"#### TIME TO FINISH SIMULATION STUDY ####"
  TEXT01<-"#### RESULT: FULL QMLE ####" 
  TEXT02<-"## COUNT NON-ZERO FOR FULL QMLE ##"
  TEXT03<-"## PARAMETER ESTIMATES: FULL QMLE ##"
  TEXT04<-"#### RESULT: ADAPTIVE LASSO ####"  
  TEXT05<-"## COUNT NON-ZERO FOR ADAPTIVE LASSO ##"
  TEXT06<-"## PARAMETER ESTIMATES: ADAPTIVE LASSO ##"

LISTFULL0<-list(TEXTtime, simuTIME , 
	            TEXT01, TEXT02, min(cNonZerofull), max(cNonZerofull), round(mean(cNonZerofull)), Mode(cNonZerofull),
	            sum(cBIASNonZerofull)/REPLICATE, sqrt(sum(cEVARNonZerofull)/REPLICATE),
	            table(pickAllSigfull), table(isTRUEinfull), 
	            TEXT03, PARAmeanQMLE, rowMins(colPARAmeanQMLE), rowMedians(colPARAmeanQMLE),
	            rowMaxs(colPARAmeanQMLE), PARAbiasQMLE, PARAesdQMLE, PARApropoQMLE,
                TEXT04, TEXT05, min(cNonZeroadaplasso), max(cNonZeroadaplasso), round(mean(cNonZeroadaplasso)), Mode(cNonZeroadaplasso),
                sum(cBIASNonZeroadaplasso)/REPLICATE, sqrt(sum(cEVARNonZeroadaplasso)/REPLICATE),
                table(pickAllSigadaplasso), table(isTRUEinadaplasso),
                TEXT06, PARAmeanadapLASSO, rowMins(colPARAmeanadapLASSO), rowMedians(colPARAmeanadapLASSO), 
                rowMaxs(colPARAmeanadapLASSO), PARAbiasadapLASSO,  PARAesdadapLASSO, PARApropadapLASSO
             )

names(LISTFULL0)<-c("", "Simulation time", 
              "", "", "Minimum", "Maximum", "Average", "Mode",
              "Bias", "ESD",
              "Pick All Relevant", "Is Correct model??", 
              "", "Est. Avg. Full", "Est. Min. Full", "Est. Med. Full",
              "Est. Max. Full", "Est. Bias Full", "Est. ESD Full", "Est. propo Full",
                "", "", "Minimum", "Maximum", "Average", "Mode",
                "Bias", "ESD",
                "Pick All Relevant", "Is Correct model??",
                "", "Est. Avg. A.LASSO", "Est. Min. A.LASSO", "Est. Med. A.LASSO", 
                "Est. Max. A.LASSO", "Est. Bias A.LASSO",  "Est. ESD A.LASSO", "Est. propo A.LASSO"
            )

OUTA<-capture.output(list(LISTFULL0)); 
cat("n=",SamSize,"Replication=",maxREP,OUTA,file="./Results/Ver1/Model 2/Model2n600q4s6.txt",sep="\n", append = FALSE)
