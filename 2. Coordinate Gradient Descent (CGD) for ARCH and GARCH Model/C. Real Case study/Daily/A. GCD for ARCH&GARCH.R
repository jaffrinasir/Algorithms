#### LOAD R LIBRARY  ########################################################
library(mvtnorm)
library(MASS)
library(fBasics)
library(matrixcalc)
#############################################################################

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
D2SIGMAa<-STEPD1(PARA,DSIGMAa,ERR2i)
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
Qq<- HEn
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
print(norm(matrix(BEFORE-PARAs1CGD),type = "1"))
if( (norm(matrix(BEFORE-PARAs1CGD),type = "1") < 0.0001) || (CYCLEGCD==10000) ){
  ITERFULL<-FALSE
}
CYCLEGCD<-CYCLEGCD+1
  }
return(PARAs1CGD)
}

#### START STEP A. INITIALIZATION ###########################################
rC<-15
sC<-0
Dt<-(1+rC+sC)
ERR2<-as.vector(Return)^2                 #COMPUTE SQUARE ERRORS
h<-max(rC,sC)
#### END STEP A. ############################################################
PARAinit<-c(mean(ERR2),rep(0, (Dt-1)))
LAMBDAa<-0
AWEIGHTa<-c(0,rep(1,rC),rep(1,sC))
PARAFULL<-GARCHCGD(LAMBDAa,AWEIGHTa,PARAinit)
### END : FULL CONSTRAINT QMLE ##############################################

#########** BEGIN LASSO & ADAPTIVE LASSO ***#################################
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
print(KK)
}
print(cPARAadaplasso[[which(BICadaplasso==min(BICadaplasso))[1]]])