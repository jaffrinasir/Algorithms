library(mvtnorm)
library(MASS)
library(fBasics)
library(matrixcalc)

library(mvtnorm)
library(MASS)
library(fBasics)
library(matrixcalc)

### LOAD DATASET ARCH(12) ####################################################
ARCH<-function(n,SEED){
  #A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
  rCh<- 12
  n<-n+200  #Obs. below 201 will be deleted 
  Q<- 1 #Number of Variate
  
  #######  Generate r.v of white noise ########################################
  eta<-rmvnorm(n, mean=c(0), diag(1)) 
  eta<-t(eta)
  #############################################################################
  
  ##INTERCEPT
  paraGA0<-matrix(c(0.01),nrow=Q,byrow=TRUE)
  #### MA(p) representative
  paraGAa1<-matrix(c(0.15),nrow=Q,byrow=TRUE)
  paraGAa2<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa3<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa4<-matrix(c(0.3),nrow=Q,byrow=TRUE)
  paraGAa5<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa6<-matrix(c(0.2),nrow=Q,byrow=TRUE)
  paraGAa7<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa8<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa9<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa10<-matrix(c(0.15),nrow=Q,byrow=TRUE)
  paraGAa11<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa12<-matrix(c(0.19),nrow=Q,byrow=TRUE)

  H<-max(rCh)
  ht<-abs(t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1)))))  #HETEROS. OBS
  err<-t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1)))) #ERR. TERM
  Y <-t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1))))
  ZERO<-matrix(0,nrow=Q,ncol=(dim(eta)[2]-dim(Y)[2]))
  ht<-cbind(ht,ZERO) 
  err<-cbind(err,ZERO)
  Y<-cbind(Y,ZERO) 

ARCHr<-sapply( paste('paraGAa', 1:rCh, sep=''), get ,envir=sys.frame(sys.parent(0)))

    for (i in (H+1):((dim(eta)[2]))) {
      sumRg<-0
       for(j in 1:rCh){
        sumRg<-sumRg + (ARCHr[[j]]%*%(err[,(i-j)]^{2}))
       }
    ht[,i]<-paraGA0+(sumRg)
    err[,i]<-eta[,i]%*%sqrt(ht[,i])
    Y[,i]<-err[,i]
  }
Y<-Y[,c(-1:-200)]; #THIS TO REMOVE FIRST 200 OBSERVATIONS TO AVOID SET.SEED EFFECT
Y1<-as.matrix(Y)

coefARCH<-cbind(paraGA0, paraGAa1, paraGAa2, paraGAa3,paraGAa4,
                         paraGAa5, paraGAa6, paraGAa7,paraGAa8,
                         paraGAa9, paraGAa10, paraGAa11,paraGAa12)

LAGNONZERO<-c(1,2,5,7,11,13)   #LOCATION OF NONZERO PARAMETERS
LAGNONZERO2<-c(1,1,0,0,1,0,1,0,0,0,1,0,1)  #1: non-zero 0: zero
ALLDATA<-list(Y1,err,coefARCH,rCh,LAGNONZERO,LAGNONZERO2,sum(LAGNONZERO2))
names(ALLDATA)<-c("Y1","err","coefARCH","Gr","LAGNONZERO","LAGNONZERO2","COUNTNZ")
return(ALLDATA)
}
##################################################################################################
ARCH1<-ARCH(1200)
ERR<-as.vector(ARCH1$Y1)
#### START STEP A. INITIALIZATION ###########################################
rC<-1
sC<-3
Dt<-(1+rC+sC)
ERR2<-as.vector(ERR)^2                 #COMPUTE SQUARE ERRORS
h<-max(rC,sC)
#### END STEP A. ############################################################

#### BAYESIAN INFORMATION CRITERION (BIC) ###################################
BICFUNC<-function(SIGMASQ0,PARA,Ne,DELTA){
  return(
    (sum(log(SIGMASQ0)))+(sum(PARA != 0)*log(Ne)*DELTA)
  )
}
#LIKELIHOOD IS MAXIMIZED, NOT MINIMIZED
#############################################################################

#### START STEP B. ##########################################################
STEPB0<-function(PARA){   ### RECURSIVE COMPUTATION OF SIGMA^2_{T}

if( (rC>0) && (sC==0) ){     ### ARCH(rC)
SIGMASQb<-numeric()
ERR2A<-SIGMASQA<-rep(0,h)
ERR2A<-c(ERR2A,ERR2)    #MERGING INITIAL VALUES + PRE-COMPUTED SQUARED ERRORS

 	for( t in 1:length(ERR2)){
    SUM1<-0
    for(ii in 1:rC){
    SUM1<- SUM1 + (PARA[ii+1]*ERR2A[t-ii+h])
	}
	SIGMASQb[t]<-SIGMASQA[t+h]<-PARA[1]+SUM1
  }
}else if( (rC>0) && (sC>0) ){ ### GARCH(rC,sC)
SIGMASQb<-numeric()
ERR2A<-SIGMASQA<-rep(0,h)
ERR2A<-c(ERR2A,ERR2)    #MERGING INITIAL VALUES + PRE-COMPUTED SQUARED ERRORS
for( t in 1:length(ERR2)){
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

#### *** LIKELIHOOD FUNC ######################################################
LIKELIHOOD<-function(PARA,LAMDDA1,AWEIGHT1){ 
  SUMB<-0
  L0<-0
  SIGMASQ0<-STEPB0(PARA)
  for(t in (h+1):length(ERR2) ){
    L0 <- L0 +(
      -(log(SIGMASQ0[t]) + (ERR2[t]/(SIGMASQ0[t])))
    )
  }
  for(ii in 1:Dt){
    SUMB <- SUMB - (LAMDDA1*AWEIGHT1[ii]*abs(PARA[ii]))
  }
  L0<-L0+SUMB
  return(L0)	
}
###############################################################################

#### START STEP C. ############################################################
STEPC0<-function(PARA,SIGMASQa){   ### RECURSIVE COMPUTATION OF FIRST DERIVATIVE SIGMA^2_{T}
if( (rC>0) && (sC==0) ){     ### ARCH(rC)
	 ERR2A<-rep(ERR2[1],h)
     ERR2A<-c(ERR2A,ERR2)    #MERGING INITIAL VALUES + SQUARED ERRORS
     V<-DSIGMA<-list()
  ########################################
     for( t in 1:length(ERR2)){
     	Ep<-numeric()
        for(ii in 1:rC){
           Ep[ii]<-ERR2A[t+h-ii]
        }
        V[[t]]<-as.matrix(c(1,Ep))
      }
  ########################################
     for( t in 1:length(ERR2)){
     	DSIGMA[[t]]<-V[[t]]
     }
}else if( (rC>0) && (sC>0) ){ ### GARCH(rC,sC)
      ERR2A<-SIGMASQA<-rep(ERR2[1],h) #INITIAL VALUES FOLLOWS SUGGESTION FROM FRANQ AND ZAKOIAN (2003)
      ERR2A<-c(ERR2A,ERR2)    #MERGING INITIAL VALUES + SQUARED ERRORS
      SIGMASQA<-c(SIGMASQA,SIGMASQa) #MERGING INITIAL VALUES + SIGMA^{2}
      V<-DSIGMA<-DSIGMAt<-list()    #STORE VALUES IN LIST
      for(kk in 1:h){
       DSIGMAt[[kk]]<-matrix(c(0),nrow=Dt)
         }
  ########################################
     for( t in 1:length(ERR2)){
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
      for( t in 1:length(ERR2)){
         SUM1<-0
        for(kk in (rC+1):(Dt-1)){ SUM1<- SUM1 + (PARA[kk+1]*DSIGMAt[[t-kk+rC+h]])}    #THIS IS MY FIRST MISTAKE *SUM1<- (----) is now corrected to SUM1<- SUM1 + (-----)
      DSIGMA[[t]]<-DSIGMAt[[t+h]]<-V[[t]] + SUM1
      }
  ########################################    
 }
LIST1<-list(V,DSIGMA)
return(LIST1) 
}
#### END STEP C. ##############################################################

#### START STEP D1. ##########################################################
STEPD1<-function(PARA,DSIGMAa){   ### RECURSIVE COMPUTATION OF SECOND DERIVATIVE SIGMA^2_{T}  ### FOR GARCH MODEL ONLY
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

for( t in 1:length(ERR2)){
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
STEPD0<-function(PARA,SIGMASQa,DSIGMAa){  ## COMPUTE THE SCORE VECTOR (Sn), TRUE HESSIAN (HEn), INFOMATION MATRIX (IMn) AND OUTER PRODUCT MATRIX (OMn)
if( (rC>0) && (sC==0) ){     ### ARCH(rC)
Sn<-HEn<-IMn<-OMn<-0 
	for(t in (h+1):length(ERR2) ){
      Dlike<-(1/2)*((1/SIGMASQa[[t]])*( (ERR2[t]/SIGMASQa[[t]]) - 1 )*DSIGMAa[[t]])
##SCORE MATRIX    
      Sn<-Sn+(Dlike)

##FULL HESSIAN
      HEn<-HEn+(
        ( 
          ( (1/(SIGMASQa[[t]])^2)*(1- ((2*ERR2[t])/(SIGMASQa[[t]]))) )
               *(DSIGMAa[[t]]%*%t(DSIGMAa[[t]]))
            )
        )

##INFO MATRIX 
      IMn<-IMn-(
               (((1/(SIGMASQa[[t]]^2)) ))*(DSIGMAa[[t]]%*%t(DSIGMAa[[t]])) 
             )

##OUTER PRODUCT MATRIX
      OMn<-OMn-(
        Dlike%*%t(Dlike)
      )
	}
}else if( (rC>0) && (sC>0) ){ ### GARCH(rC,sC)
D2SIGMAa<-STEPD1(PARA,DSIGMAa)
Sn<-HEn<-IMn<-OMn<-HEn2<-0 
  for(t in (h+1):length(ERR2) ){
      Dlike<-((1/SIGMASQa[[t]])*( (ERR2[t]/SIGMASQa[[t]]) - 1 )*DSIGMAa[[t]])
##SCORE MATRIX    
      Sn<-Sn+(Dlike)

##FULL HESSIAN 1
      HEn<-HEn+(
        (
          (  ((ERR2[t])/(SIGMASQa[[t]])^2) -  (1/(SIGMASQa[[t]]))  )*D2SIGMAa[[t]]

        )
          +
        ( 
          ( (1/(SIGMASQa[[t]])^2)*(1- ((2*ERR2[t])/(SIGMASQa[[t]]))) )
               *(DSIGMAa[[t]]%*%t(DSIGMAa[[t]]))
            )
        )

##INFO MATRIX 
      IMn<-IMn-(
               (((1/(SIGMASQa[[t]]^2)) ))*(DSIGMAa[[t]]%*%t(DSIGMAa[[t]])) 
             )

##OUTER PRODUCT MATRIX
      OMn<-OMn-(
        Dlike%*%t(Dlike)
      )

  }
}
LIST2<-list(Sn,HEn,IMn,OMn)
return(LIST2)
 }
#### END STEP D. ##############################################################

#### STEP E. ARMIJO RULE COMP. ################################################
STEPE0<-function(PARA,LAMBDA1,AWEIGHT1,L){  ## DETERMINATION OF \tau^{K-1}_s, THE STEP SIZE
ArM0<- 1 ##SET INITIAL STEP SIZE
ITERArmijo<-TRUE
Alp0<-0.1 #\tilde{\alpha} # default is 0.5 or 0.00001
Del0<-0.01     #\tilde{\delta} #default is 0.8 or 0.5
Gam0<-1                                                  #(DEFAULT AS IN TSENG & YUN (2009), pg 414)
D0<-t(dB)%*%(In%*%dB)

LiKE1<-LIKELIHOOD(PARA,LAMBDA1,AWEIGHT1)
while(ITERArmijo==TRUE){
if( (PARA[L]+(ArM0*(dB[L]))) < 0 ){
  ArM0<-(Del0*ArM0)
  if( (ArM0 <= 1e-8 ) ){ITERArmijo<-FALSE}
  }else{
LiKE2<-LIKELIHOOD((PARA+(ArM0*(dB))),LAMBDA1,AWEIGHT1)
DIRECTION<-Alp0*ArM0*(Gam0-1)*(D0)
############ COMPUTE ARMIJO RULE #####################
   if( ((LiKE1+DIRECTION) > LiKE2) && (ArM0 > 1e-8 )){
    ArM0<-(Del0*ArM0)
     }else{ITERArmijo<-FALSE}
   }
}
return(ArM0)
}
#### END STEP E. ##############################################################

##**** BEGIN: CGD ALGORITHM, FULL QMLE ########################################
PARAs1<-c(mean(ERR2),rep(0, (Dt-1)))
LAMBDA<-0
AWEIGHT<-c(0,rep(1,rC),rep(1,sC))
############## (A) GCD WITH LAMBDA=0 ##########################################
ITERFULL<-TRUE
CYCLEGCD<-0

while(ITERFULL==TRUE){
BEFORE<-PARAs1

for(L in 1:Dt){ ###################### CYCLE FORM 1 to D
####################### STEP B & C
SIGMASQ<-STEPB0(PARAs1)
STEPCresult<-STEPC0(PARAs1,SIGMASQ)
DSIGMA<-STEPCresult[[2]]
####################### STEP D
STEPDresult<-STEPD0(PARAs1,SIGMASQ,DSIGMA)
Sn<- STEPDresult[[1]]
HEn<-STEPDresult[[2]]
### Q matrix 
In<-HEn
##################################

dB<-rep(0,Dt)
  fA<-as.matrix(Sn[L]+(In[L,L]*(-PARAs1[L])))
    if(norm(fA,"1")<=(LAMBDA*AWEIGHT[L])){
      dB[L]<- -PARAs1[L]
    }else{ 
      if(fA > (LAMBDA*AWEIGHT[L])){
      dB[L]<- ((LAMBDA*AWEIGHT[L])-Sn[L])/In[L,L]
       }else if(fA < -(LAMBDA*AWEIGHT[L])){
      dB[L]<- (-(LAMBDA*AWEIGHT[L])-Sn[L])/In[L,L]
      }
    }
####################### STEP D
####################### UPDATE PARAMETER
if(norm(as.matrix(dB),"1")!=0){
ArM<-STEPE0(PARAs1,LAMBDA,AWEIGHT,L)
PARAs1[L]<-max((PARAs1[L]+(ArM*(dB[L]))),0)
  }
}##################################### END CYCLE
print(norm(matrix(BEFORE-PARAs1),type = "1"))
if( (norm(matrix(BEFORE-PARAs1),type = "1") < 0.0001) || (CYCLEGCD==1000) ){
  ITERFULL<-FALSE
}
CYCLEGCD<-CYCLEGCD+1
}