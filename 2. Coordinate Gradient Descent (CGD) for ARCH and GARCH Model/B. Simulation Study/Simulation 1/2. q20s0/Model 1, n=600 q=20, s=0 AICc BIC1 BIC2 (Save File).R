setwd("/home/jaffri/Documents/Simulation")
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
minREP<-1
maxREP<-1009
TRAIN<-90 #IN PERCENTAGE
TEST<-10 #IN PERCENTAGE
EXCLUDE<-c(10,29,237, 356,465,525,826,905,945)
### LOAD DATASET ARCH(12) ####################################################
ARCH<-function(n, REPLI){
  #A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
  rCh<- 12
  n<-n+200  #Obs. below 201 will be deleted 
  Q<- 1 #Number of Variate
  
  #######  Generate r.v of white noise ########################################
  set.seed(REPLI)
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
                         paraGAa9, paraGAa10, paraGAa11,paraGAa12,matrix(0,ncol=8))

LAGNONZERO<-c(1,2,5,7,11,13)   #LOCATION OF NONZERO PARAMETERS
LAGNONZERO2<-c(1,1,0,0,1,0,1,0,0,0,1,0,1,rep(0,8))  #1: non-zero 0: zero
ALLDATA<-list(Y1,err,coefARCH,rCh,LAGNONZERO,LAGNONZERO2,sum(LAGNONZERO2))
names(ALLDATA)<-c("Y1","err","coefARCH","Gr","LAGNONZERO","LAGNONZERO2","COUNTNZ")
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
ICFUNC<-function(SIGMASQ0,ERRSQ0,PARA,Ne,OPTION){
  return(
    if(OPTION==1){
      2*(sum(log(SIGMASQ0)+(ERRSQ0/SIGMASQ0)))+(sum(PARA != 0)*2*Ne/(Ne-sum(PARA != 0)-1))
    }
    else if(OPTION==2){
      2*(sum(log(SIGMASQ0)+(ERRSQ0/SIGMASQ0)))+(sum(PARA != 0)*log(Ne)*1)
      }else if(OPTION==3){
      2*(sum(log(SIGMASQ0)+(ERRSQ0/SIGMASQ0)))+(sum(PARA != 0)*log(Ne)*2)
      }
  )
}
#LIKELIHOOD IS MINIMIZED (NEGATIVE LIKELIHOOD)

#PREDICTION
RMSE0<-function(ERR2M,SIGMASQM){
  if(length(ERR2M)==length(SIGMASQM)){
    return(
        sqrt(sum((ERR2M-SIGMASQM)^2)/length(ERR2M))
      )
  }
}

MAPE0<-function(ERR2M,SIGMASQM){
  if(length(ERR2M)==length(SIGMASQM)){
    return(
        sum(abs((ERR2M-SIGMASQM)/ERR2M))/length(ERR2M)
      )
  }
}

MAE0<-function(ERR2M,SIGMASQM){
  if(length(ERR2M)==length(SIGMASQM)){
    return(
        sum(abs(ERR2M-SIGMASQM))/length(ERR2M)
      )
  }
}

STEPPRED<-function(PARA,ERR2i){   ### RECURSIVE COMPUTATION OF SIGMA^2_{T}
if( (rC>0) && (sC==0) ){     ### ARCH(rC)
ERR2A<-SIGMASQA<-rep(0,h)
ERR2A<-c(ERR2A,ERR2i)    #MERGING INITIAL VALUES + PRE-COMPUTED SQUARED ERRORS

  for( t in 1:length(ERR2i)){
    SUM1<-0
    for(ii in 1:rC){
    SUM1<- SUM1 + (PARA[ii+1]*ERR2A[t-ii+h])
  }
  SIGMASQA[t+h]<-PARA[1]+SUM1
  }
}else if( (rC>0) && (sC>0) ){ ### GARCH(rC,sC)
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
  SIGMASQA[t+h]<-PARA[1]+SUM1+SUM2
  } #END FOR
}
return(list(ERR2A[-1:-h],SIGMASQA[-1:-h]))
} #END STEPB0

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
if( (norm(matrix(BEFORE-PARAs1CGD),type = "1") < 0.001) || (CYCLEGCD==10000) ){
  ITERFULL<-FALSE
}
CYCLEGCD<-CYCLEGCD+1
  }
return(PARAs1CGD)
}
###############################################################################


### START : SIMULATION ########################################################
LISTNO<-c(minREP:maxREP)
LISTNO<-LISTNO[! LISTNO %in% EXCLUDE]
for(REPLICATE in LISTNO){
cat("REPLICATION=", REPLICATE,"\n")
rC<-20
sC<-0
Dt<-(1+rC+sC)
### STEP 1: generate dataset and other values. #################################
ARCH1<-ARCH(SamSize,REPLICATE)
ERR<-as.vector(ARCH1$Y1)[1:(TRAIN/100*SamSize)]
ERRa<-as.vector(ARCH1$Y1)[((TRAIN/100*SamSize)+1):SamSize]
ERR2<-as.vector(ERR)^2                 #COMPUTE SQUARE ERRORS
ERR2a<-as.vector(ERRa)^2 
h<-max(rC,sC)
  TRUEPARA<-ARCH1$coefARCH
  NONZERO<-ARCH1$LAGNONZERO
  CTRUENONZERO<-ARCH1$COUNTNZ
  propoTrue<-ARCH1$LAGNONZERO2
### BEGIN : FULL CONSTRAINT QMLE ############################################
PARAinit<-c(mean(ERR2),rep(0, (Dt-1)))
LAMBDAa<-0
AWEIGHTa<-c(0,rep(1,rC),rep(1,sC))
start_time <- Sys.time()
PARAFULL<-GARCHCGD(LAMBDAa,AWEIGHTa,PARAinit)
end_time <- Sys.time()
### COLLECT RESULT: FULL QMLE
A1<-PARAFULL
A2<-(PARAFULL-TRUEPARA)
A3<-(PARAFULL-TRUEPARA)^2
A4<-all(NONZERO%in%which(PARAFULL!=0))
A5<-isTRUE(all.equal(NONZERO,which(PARAFULL!=0)))
A6<-sum(as.integer(PARAFULL!=0))
A7<- CTRUENONZERO-sum(as.integer(PARAFULL!=0))
A8<- (CTRUENONZERO-sum(as.integer(PARAFULL!=0)))^2
A9<-as.integer(as.integer(PARAFULL!=0))
A10<-RMSE0(STEPPRED(PARAFULL,ERR2a)[[1]],STEPPRED(PARAFULL,ERR2a)[[2]])
A11<-MAPE0(STEPPRED(PARAFULL,ERR2a)[[1]],STEPPRED(PARAFULL,ERR2a)[[2]])
A12<-MAE0(STEPPRED(PARAFULL,ERR2a)[[1]],STEPPRED(PARAFULL,ERR2a)[[2]])
### END : FULL CONSTRAINT QMLE ##############################################

#########** BEGIN ADAPTIVE LASSO ############################################
PARAWEIGHT<-PARAFULL
PARAWEIGHT[PARAWEIGHT==0]<-1e-20
AICcadaplasso<-BIC1adaplasso<-BIC2adaplasso<-numeric()
PARAadlasso<-c(mean(ERR2),rep(0, (Dt-1)))
cPARAadaplasso<-list()
AWEIGHTb<-c(0,(PARAWEIGHT[2:(Dt)]^(-1)))
LAMBDAlist<-c(seq(20, 1, length.out= 20), seq(0.9, 0, length.out= 20))

start_timeA <- Sys.time()
for(KK in 1:length(LAMBDAlist)){
LAMBDAb<-LAMBDAlist[KK]
if(KK==1){
PARAinit2<-PARAadlasso<-cPARAadaplasso[[KK]]<-GARCHCGD(LAMBDAb,AWEIGHTb,PARAinit)
}else{
PARAinit2<-PARAadlasso<-cPARAadaplasso[[KK]]<-GARCHCGD(LAMBDAb,AWEIGHTb,PARAinit2)	
}
### COMPUTE IC ################################################################
AICcadaplasso[KK]<-ICFUNC(STEPB0(cPARAadaplasso[[KK]],ERR2),ERR2,cPARAadaplasso[[KK]],length(ERR2),1)
BIC1adaplasso[KK]<-ICFUNC(STEPB0(cPARAadaplasso[[KK]],ERR2),ERR2,cPARAadaplasso[[KK]],length(ERR2),2)
BIC2adaplasso[KK]<-ICFUNC(STEPB0(cPARAadaplasso[[KK]],ERR2),ERR2,cPARAadaplasso[[KK]],length(ERR2),3)
### END STEP 3B ################################################################
}
end_timeA <- Sys.time()
##### AICc
RePARAadaplassoAICc<-cPARAadaplasso[[which(AICcadaplasso==min(AICcadaplasso))[1]]]
#### COLLECT RESULT: ADAPTIVE LASSO
B1<-RePARAadaplassoAICc
B2<-(RePARAadaplassoAICc-TRUEPARA)
B3<-(RePARAadaplassoAICc-TRUEPARA)^2
B4<-all(NONZERO%in%which(RePARAadaplassoAICc!=0))
B5<-isTRUE(all.equal(NONZERO,which(RePARAadaplassoAICc!=0)))
B6<-sum(as.integer(RePARAadaplassoAICc!=0))
B7<- CTRUENONZERO-sum(as.integer(RePARAadaplassoAICc!=0))
B8<- (CTRUENONZERO-sum(as.integer(RePARAadaplassoAICc!=0)))^2
B9<-as.integer(as.integer(RePARAadaplassoAICc!=0))
B10<-RMSE0(STEPPRED(RePARAadaplassoAICc,ERR2a)[[1]],STEPPRED(RePARAadaplassoAICc,ERR2a)[[2]])
B11<-MAPE0(STEPPRED(RePARAadaplassoAICc,ERR2a)[[1]],STEPPRED(RePARAadaplassoAICc,ERR2a)[[2]])
B12<-MAE0(STEPPRED(RePARAadaplassoAICc,ERR2a)[[1]],STEPPRED(RePARAadaplassoAICc,ERR2a)[[2]])

##### BIC1
RePARAadaplassoBIC1<-cPARAadaplasso[[which(BIC1adaplasso==min(BIC1adaplasso))[1]]]
#### COLLECT RESULT: ADAPTIVE LASSO
C1<-RePARAadaplassoBIC1
C2<-(RePARAadaplassoBIC1-TRUEPARA)
C3<-(RePARAadaplassoBIC1-TRUEPARA)^2
C4<-all(NONZERO%in%which(RePARAadaplassoBIC1!=0))
C5<-isTRUE(all.equal(NONZERO,which(RePARAadaplassoBIC1!=0)))
C6<-sum(as.integer(RePARAadaplassoBIC1!=0))
C7<- CTRUENONZERO-sum(as.integer(RePARAadaplassoBIC1!=0))
C8<- (CTRUENONZERO-sum(as.integer(RePARAadaplassoBIC1!=0)))^2
C9<-as.integer(as.integer(RePARAadaplassoBIC1!=0))
C10<-RMSE0(STEPPRED(RePARAadaplassoBIC1,ERR2a)[[1]],STEPPRED(RePARAadaplassoBIC1,ERR2a)[[2]])
C11<-MAPE0(STEPPRED(RePARAadaplassoBIC1,ERR2a)[[1]],STEPPRED(RePARAadaplassoBIC1,ERR2a)[[2]])
C12<-MAE0(STEPPRED(RePARAadaplassoBIC1,ERR2a)[[1]],STEPPRED(RePARAadaplassoBIC1,ERR2a)[[2]])

##### BIC2
RePARAadaplassoBIC2<-cPARAadaplasso[[which(BIC2adaplasso==min(BIC2adaplasso))[1]]]
#### COLLECT RESULT: ADAPTIVE LASSO
D1<-RePARAadaplassoBIC2
D2<-(RePARAadaplassoBIC2-TRUEPARA)
D3<-(RePARAadaplassoBIC2-TRUEPARA)^2
D4<-all(NONZERO%in%which(RePARAadaplassoBIC2!=0))
D5<-isTRUE(all.equal(NONZERO,which(RePARAadaplassoBIC2!=0)))
D6<-sum(as.integer(RePARAadaplassoBIC2!=0))
D7<- CTRUENONZERO-sum(as.integer(RePARAadaplassoBIC2!=0))
D8<- (CTRUENONZERO-sum(as.integer(RePARAadaplassoBIC2!=0)))^2
D9<-as.integer(as.integer(RePARAadaplassoBIC2!=0))
D10<-RMSE0(STEPPRED(RePARAadaplassoBIC2,ERR2a)[[1]],STEPPRED(RePARAadaplassoBIC2,ERR2a)[[2]])
D11<-MAPE0(STEPPRED(RePARAadaplassoBIC2,ERR2a)[[1]],STEPPRED(RePARAadaplassoBIC2,ERR2a)[[2]])
D12<-MAE0(STEPPRED(RePARAadaplassoBIC2,ERR2a)[[1]],STEPPRED(RePARAadaplassoBIC2,ERR2a)[[2]])

### END : ADAPTIVE LASSO  ##############################################
LISTsimu<-list(
  list(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12),
  list(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12),
  list(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12),
  list(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12),
  list(end_time-start_time, end_timeA-start_timeA)
    )
LISTA<-list(assign(paste("LIST",REPLICATE,sep=""),LISTsimu))
save(LISTA,file=paste("Rlistq20s0/600/SIMUAIC600MODEL1.(",REPLICATE,").RData",sep = ""))
}
