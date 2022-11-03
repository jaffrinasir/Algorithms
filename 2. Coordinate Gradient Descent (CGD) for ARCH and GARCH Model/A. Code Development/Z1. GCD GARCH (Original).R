GCDGARCH<-function(rCi,sCi,ERR2i,PARAINIT){
  Dt<-rCi+sCi+1
  h<-max(rCi,sCi)
#### *** LIKELIHOOD FUNC ####################################################
F0<-function(PARA,ERR2i){
  L0<-0
  SIGMASQ0<-STEPB0(rCi,sCi,PARA,ERR2i)
  for(t in (h+1):length(ERR2i) ){
    L0 <- L0 +(
      (log(SIGMASQ0[t]) + (ERR2i[t]/(SIGMASQ0[t])))
    )
  }
  return(L0)	
}
###############################################################################

#### START STEP B. ##########################################################
STEPB0<-function(rCi,sCi,PARA,ERR2i){   ### RECURSIVE COMPUTATION OF SIGMA^2_{T}
if( (rCi>0) && (sCi==0) ){     ### ARCH(rCi)
SIGMASQb<-numeric()
ERR2A<-SIGMASQA<-rep(0,h)
ERR2A<-c(ERR2A,ERR2i)    #MERGING INITIAL VALUES + PRE-COMPUTED SQUARED ERRORS

  for( t in 1:length(ERR2i)){
    SUM1<-0
    for(ii in 1:rCi){
    SUM1<- SUM1 + (PARA[ii+1]*ERR2A[t-ii+h])
  }
  SIGMASQb[t]<-SIGMASQA[t+h]<-PARA[1]+SUM1
  }
}else if( (rCi>0) && (sCi>0) ){ ### GARCH(rCi,sCi)
SIGMASQb<-numeric()
ERR2A<-SIGMASQA<-rep(0,h)
ERR2A<-c(ERR2A,ERR2i)    #MERGING INITIAL VALUES + PRE-COMPUTED SQUARED ERRORS
for( t in 1:length(ERR2i)){
  SUM1<-SUM2<-0
  for(ii in 1:rCi){
    SUM1<- SUM1 + (PARA[ii+1]*ERR2A[t-ii+h])
  }
  for(jj in (rCi+1):(Dt-1)){
    SUM2<- SUM2 + (PARA[jj+1]*SIGMASQA[t-jj+rCi+h])
  }
  SIGMASQb[t]<-SIGMASQA[t+h]<-PARA[1]+SUM1+SUM2
  } #END FOR
}
return(SIGMASQb)
} #END STEPB0
#### END STEP B. ##############################################################

#### START STEP C. ############################################################
STEPC0<-function(rCi,sCi,PARA,SIGMASQa,ERR2i){   ### RECURSIVE COMPUTATION OF FIRST DERIVATIVE SIGMA^2_{T}
if( (rCi>0) && (sCi==0) ){     ### ARCH(rCi)
   ERR2A<-rep(ERR2i[1],h)
     ERR2A<-c(ERR2A,ERR2i)    #MERGING INITIAL VALUES + SQUARED ERRORS
     V<-DSIGMA<-list()
  ########################################
     for( t in 1:length(ERR2i)){
      Ep<-numeric()
        for(ii in 1:rCi){
           Ep[ii]<-ERR2A[t+h-ii]
        }
        V[[t]]<-as.matrix(c(1,Ep))
      }
  ########################################
     for( t in 1:length(ERR2i)){
      DSIGMA[[t]]<-V[[t]]
     }
}else if( (rCi>0) && (sCi>0) ){ ### GARCH(rCi,sCi)
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
        for(ii in 1:rCi){
           Ep[ii]<-ERR2A[t+h-ii]
        }
        for(jj in 1:sCi){
           Sigq[jj]<-SIGMASQA[t+h-jj]
        }
        V[[t]]<-as.matrix(c(1,Ep,Sigq))  #MERGING VECTORS
      }
  ########################################
      for( t in 1:length(ERR2i)){
         SUM1<-0
        for(kk in (rCi+1):(Dt-1)){ SUM1<- SUM1 + (PARA[kk+1]*DSIGMAt[[t-kk+rCi+h]])} 
      DSIGMA[[t]]<-DSIGMAt[[t+h]]<-V[[t]] + SUM1
      }
  ########################################    
 }
LIST1<-list(V,DSIGMA)
return(LIST1) 
}
#### END STEP C. ##############################################################

#### START STEP D1. ##########################################################
STEPD1<-function(rCi,sCi,PARA,DSIGMAa,ERR2i){   ### RECURSIVE COMPUTATION OF SECOND DERIVATIVE SIGMA^2_{T}  ### FOR GARCH MODEL ONLY
DBETA<-D2SIGMAt<-DSIGMAt2<-D2SIGMAt2<-list()
ZEROMATD<-matrix(c(0),ncol=Dt,nrow=Dt)
##########################################################
      for(kk in 1:sCi){
        DBETA[[kk]]<-matrix(c(0),ncol=Dt)
        DBETA[[kk]][1,(kk+rCi+1)]<-1
      }

      for(kk in 1:h){
       DSIGMAt2[[kk]]<-matrix(c(0),nrow=Dt)   #INITIAL VALUE OF FIRST DERIVATIVE
       D2SIGMAt2[[kk]]<-matrix(c(0),ncol=Dt,nrow=Dt) #INITIAL VALUE OF THE SECOND DERIVATIVE
         }
DSIGMAt2<-c(DSIGMAt2,DSIGMAa)    #DSIGMA FROM STEP C.
##########################################################
for( t in 1:length(ERR2i)){
   A<-ZEROMATD
   for(i in 1:sCi){
    A[(i+1+rCi),]<-DSIGMAt2[[t+h-i]]
   }

   B<-ZEROMATD
   for(j in 1:(sCi)){
     B<-B+(DSIGMAt2[[t+h-j]]%*%DBETA[[j]])
   }

   C<-ZEROMATD
   for(k in (rCi+1):(Dt-1)){
     C<- C+(PARA[k+1]*D2SIGMAt2[[t+h-k+rCi]])
   }
D2SIGMAt[[t]]<-D2SIGMAt2[[t+h]] <- (A+B+C)
    }
return(D2SIGMAt)
}
#### END STEP C(1). ##############################################################

#### START STEP D. ############################################################
STEPD0<-function(rCi,sCi,PARA,SIGMASQa,DSIGMAa,ERR2i){  ## COMPUTE THE SCORE VECTOR (Sn), TRUE HESSIAN (HEn), 
#INFOMATION MATRIX (IMn) AND OUTER PRODUCT MATRIX (OMn)
if( (rCi>0) && (sCi==0) ){     ### ARCH(rCi)
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
}else if( (rCi>0) && (sCi>0) ){ ### GARCH(rCi,sCi)
D2SIGMAa<-STEPD1(rCi,sCi,PARA,DSIGMAa,ERR2i)
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
STEPE0<-function(PARA,L,dBi,Qq,ERR2i){  ## DETERMINATION OF \tau^{K-1}_s, THE STEP SIZE
ArM0<- 1 ##SET INITIAL STEP SIZE
ITERArmijo<-TRUE
Alp0<-0.1      #\tilde{\alpha} # default is 0.5 or 0.00001
Del0<-0.5     #\tilde{\delta} #default is 0.8 or 0.5
Gam0<-0        #(DEFAULT AS IN TSENG & YUN (2009), pg 414)
D0<-t(dBi)%*%(Qq%*%dBi)

LiKE1<-F0(PARA,ERR2i)
while(ITERArmijo==TRUE){
if( (PARA[L]+(ArM0*(dBi[L]))) < 0 ){
  ArM0<-(Del0*ArM0)
  if( (ArM0 <= 1e-8 ) ){ITERArmijo<-FALSE}
  }else{
LiKE2<-F0((PARA+(ArM0*(dBi))),ERR2i)
DIRECTION<-Alp0*ArM0*(Gam0-1)*(D0)
############ COMPUTE ARMIJO RULE #####################
   if((LiKE2 > (LiKE1+DIRECTION)) && (ArM0 > 1e-8 )){
    ArM0<-(Del0*ArM0)
     }else{ITERArmijo<-FALSE}
   }
}
return(ArM0)
}
#### END STEP E. ##############################################################

GARCHCGD<-function(rCi,sCi,PARAs1CGDi,ERR2i){
##**** BEGIN: CGD ALGORITHM, FULL QMLE ########################################
PARAs1CGD<-PARAs1CGDi
############## (A) GCD WITH LAMBDA=0 ##########################################
ITERFULL<-TRUE
CYCLEGCD<-0

while(ITERFULL==TRUE){
BEFORE<-PARAs1CGD
##################################
for(L in 1:Dt){ ###################### CYCLE FORM 1 to D
####################### STEP B & C
SIGMASQ<-STEPB0(rCi,sCi,PARAs1CGD,ERR2i)
STEPCresult<-STEPC0(rCi,sCi,PARAs1CGD,SIGMASQ,ERR2i)
DSIGMA<-STEPCresult[[2]]
####################### STEP D
STEPDresult<-STEPD0(rCi,sCi,PARAs1CGD,SIGMASQ,DSIGMA,ERR2i)
Sn<- STEPDresult[[1]]
HEn<-STEPDresult[[2]]
IMn<-STEPDresult[[3]]
OMn<-STEPDresult[[4]]/2
### SET THE SECOND DERIVATIVE OF LOG-LIKELIHOOD = Information Matrix
Qq<- HEn
dB<-rep(0,Dt)
dB[L]<- (-Sn[L])/Qq[L,L]
####################### STEP D
####################### UPDATE PARAMETER
if(norm(as.matrix(dB),"1")!=0){
ArM<-STEPE0(PARAs1CGD,L,dB,Qq,ERR2i)
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

COLLRESULT<-list()
FiNALESTPARA<-GARCHCGD(rCi,sCi,PARAINIT,ERR2i)
FiNALESTSIGMASQ<-STEPB0(rCi,sCi,FiNALESTPARA,ERR2i)
COLLRESULT[[1]]<-FiNALESTPARA
COLLRESULT[[2]]<-FiNALESTSIGMASQ
return(COLLRESULT)
}
#ERR<-as.vector(ERR)
ERR2<-ERR^2
rC<-1
sC<-2
Dt<-(1+rC+sC)
PARAinit<-c(mean(ERR2),rep(0, (Dt-1)))
GCDGARCH(rC,sC,ERR2,PARAinit)