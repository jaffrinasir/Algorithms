rm(list = ls())
#### SET SYSTEM ENVIRONMENT (Cpp) ######################
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -std=c++11")
setwd("C:/Users/Jaffri/Dropbox/Journal Papers drafts/2020-2023/1.Paper (SETAR)/Springer Nature/Simulation Files")
#Set or change the workspace directory first - otherwise the code will not execute properly
########################################################

####  LOAD REQUIRED PACKAGES ###########################
library(Matrix)
library(fBasics)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(rootSolve)
library(robustbase)
library(mvtnorm)
library(ks)
library(matrixcalc)

### LOAD DEFINITIONS      ####################################################
minREP<-1
maxREP<-250
LAMBDAlength<-1
LAMBDAlist<-seq(0.5,0.01,length.out=LAMBDAlength)  #GENERATE LIST OF LAMBDAS

DATA3A<-function(n,d,r1,r2,SEED){
#A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
library(mvtnorm)
library(MASS)
library(fBasics)
p<-2      #AR(2)
n<-n+200  #Below obs 201 will be deleted (See line 37)
Q<- 1 #Number of Variate

#############################################################################
set.seed(SEED)
ETA<-rmvnorm(n, mean=c(0), diag(1)) #Generate r.v
ETA<-t(ETA)
#############################################################################

d<- d   #d is delay parameter
r1<- r1  #r1-r2 is threshold variable
r2<- r2

paraC1<-matrix(c(0),nrow=Q,byrow=TRUE)
paraC2<-matrix(c(0),nrow=Q,byrow=TRUE)
paraC3<-matrix(c(0),nrow=Q,byrow=TRUE)
#AR parameter for FIRST regime 
paraT11<-matrix(c(0.8),nrow=Q,byrow=TRUE)
paraT12<-matrix(c(-0.2),nrow=Q,byrow=TRUE)
#AR parameter for SECOND regime 
paraT21<-matrix(c(1.9),nrow=Q,byrow=TRUE)
paraT22<-matrix(c(-0.81),nrow=Q,byrow=TRUE)
#AR parameter for THIRD regime 
paraT31<-matrix(c(0.6),nrow=Q,byrow=TRUE)
paraT32<-matrix(c(-1.0),nrow=Q,byrow=TRUE)

#HERE, SIGMA1-SIGMA3 and SD1-SD3 matrices can have different diagonal values for each regime
SIGMA1<-matrix(c(1),ncol=Q,byrow=TRUE); SD<-chol(SIGMA1); #Var-Cov matrix (Variance for univariate)
sigmaY<-matrix(c(1),ncol=Q,byrow=TRUE) 

Y <- t(rmvnorm(n=p, mean=c(0), sigma=sigmaY))   #GENERATE INITAL SERIES OF Y, #MUST TRANSPOSE
ZERO<-matrix(0,nrow=Q,ncol=(dim(ETA)[2]-dim(Y)[2]))
Y<-cbind(Y,ZERO); 

for (i in (p+1):((dim(ETA)[2]))) { 
  if(Y[,i-d]<=r1){                                       
  Y[,i]<-paraC1+(paraT11%*%(as.matrix(Y[,i-1])))+(paraT12%*%(as.matrix(Y[,i-2])))+(SD%*%ETA[,i])  
  }
  else if(Y[,i-d]> r1 && Y[,i-d]<=r2){
  Y[,i]<-paraC2+(paraT21%*%(as.matrix(Y[,i-1])))+(paraT22%*%(as.matrix(Y[,i-2])))+(SD%*%ETA[,i])  
  }else{
  Y[,i]<-paraC3+(paraT31%*%(as.matrix(Y[,i-1])))+(paraT32%*%(as.matrix(Y[,i-2])))+(SD%*%ETA[,i])  
  }
}
Y<-Y[,c(-1:-200)]; #REMOVE SET.SEED EFFECT
Y1<-as.matrix(Y)

#B.] CONSTRUCT ORDERED AUTOREGRESSION FOR STANDARD FORM 
for(par in 1:p){  #DEVELOP MODIFIED DESIGN MATRIX
  if(par==1){
    DATAY<-matrix(c(Y1[((p+1):length(Y1[,1])),]),ncol=1)
    DATAX<-matrix(Y1[(((p+1)-par):(length(Y1[,1])-par)),])
  }else{
    DATAXt<-Y1[(((p+1)-par):(length(Y1[,1])-par)),]
    DATAX<-cbind(DATAX,DATAXt)
  }
}

DATA1<-cbind(DATAY,DATAX) 

##CREATE TEMPORARILY THRESHOLD VARIABLE FOR SORTING, APPEND TO DATA1
###CREATE OPTIONS FOR SELECTING TYPES OF THRESHOLD VARIABLE
d=d   #RECALL DELAY PARAMETER
  TV1<-as.matrix(DATAX[, d])
  DATA2<-cbind(TV1,DATA1)

SSTEP<-DATA2;

#APPLYING ORDERED AUTOREGRESSION 
DATA3<-DATA2[order(DATA2[,1]),] 

tCOFF<-cbind(paraC1,paraT11,paraT12,paraC2,paraT21,paraT22,paraC3,paraT31,paraT32)

ALLDATA<-list(DATA3,Y1,d,tCOFF,SSTEP)
names(ALLDATA)<-c("DATA3","Y1","delay","tCOFF","SSTEP")
return(ALLDATA)
}

########################################################
sourceCpp("./CppFiles/1.YBySGRAM.cpp") 
# This will produce SGRAM() and YBy() C++ function
sourceCpp("./CppFiles/B1.COMPSSE.cpp")                    #FOR "COMPSSE" FUNCTIONS
sourceCpp("./CppFiles/B2.ARCOEF.cpp")                     #FOR "ARCOEF"  FUNCTION
#############################################################################################
#### LOAD R FUNCTIONS  #################################
OneSHAUS<-function(A,B){
  #ONE SIDED A -> B
  VecA<-numeric()
  for(i in 1:length(A)){
    VecB<-numeric()
    for(ii in 1:length(B)){
      VecB[ii]<-abs(A[i]-B[ii])
    }
    VecA[i]<-min(VecB)
  }
  return(max(VecA)) 
}

ICb<-function(jCard,No,SSEa,Cn){
    IC1b<-(No*log(SSEa/No))+(jCard*Cn*log(No))
  return(IC1b)
}

BICSETAR1S<-function(SAMPLEa,SSEa,Aset,Pa,Cn){   #BIC FOR SETAR, FIRST STAGE ESTIMATE
BIC1a<-((SAMPLEa)*log(SSEa/SAMPLEa))+(length(Aset)*Cn*log(SAMPLEa))
return(BIC1a)
}

SSESETAR1S<-function(Aset,BETAa){                   #SSE FOR SETAR, FIRST STAGE ESTIMATE
SSEvec<-numeric()
for(II in 1:SAMPLE){
SSEvec[II]<-(Y0[II]-(t(as.matrix(DATAC[II,]))%*%
as.matrix(rowSums(as.matrix(BETAa[,c(Aset[which(Aset<=II)])])))
))^{2}
}
return(sum(SSEvec))
}

SLSalgo<-function(LAMBDA,ActiVE,BETAab){
PRECISION<-1e-3                               #PRESSION CONSTANT (smaller constant produce slower convergence)
BETA<-BETAab            
j<-sort(ActiVE)                                    #INITIALIZATION FOR BETA

Fj<-function(INDEX,jA,B){
  if(length(jA)<=1){AC<-0}else{
      jB<-jA[! jA %in% INDEX]
      for(i3A in 1:length(jB)){
        if(i3A==1){
          AC<- XlXl[[max(jB[i3A],INDEX)]]%*%as.matrix(B[,jB[i3A]])
        }else{
          AC<-AC+(XlXl[[max(jB[i3A],INDEX)]]%*%as.matrix(B[,jB[i3A]]))
        }
      }
  }
  return(AC)
}
   
XlXl<-list()
Xlyl<-list()
vl<-list()
dl<-list()

for(i in 1:SAMPLE){
XlXl[[i]]<-SGRAM(DATAC,i,SAMPLE)     #\sum_{k=i}^{n*}x_{\pi(k)}x_{\pi(k)}^{T}
Xlyl[[i]]<-YBy(DATAC,Y0,i,SAMPLE)    #\sum_{k=i}^{n*}x_{\pi(k)}y_{\pi(k)+d}
SVDSt<-svd(XlXl[[i]])
vl[[i]]<-SVDSt$v
dl[[i]]<-SVDSt$d
}

CONDITION1<-TRUE
KK<-0
while(CONDITION1==TRUE){
  BETATempo<-BETA
  for(i4 in 1:SAMPLE){
         #j<-numeric()                                                #j is an active set
           #for(iii in 1:SAMPLE){
               #  if(norm(as.matrix(BETA[,iii]),type = c("m"))!=0){
               # j<-sort(c(j,iii))
              # }
              #}
    W<-Xlyl[[i4]]-Fj(i4,j,BETA)
    V<-vl[[i4]]%*%W
    D<-dl[[i4]]
    if(i4==1){
      BETA[,i4]<-solve(XlXl[[i4]])%*%W
      j<-sort(unique(c(j,i4)))
    }else{
      if((2*norm(W,type="2")/SAMPLE) < LAMBDA){
        BETA[,i4]<-0
        j<-j[! j %in% i4]
        }else{
        f2a<-function(u1){
             for(ij in 1:length(D)){
              if(ij==1){
                uu<-(V[ij])^2/((D[ij]*u1)+((2)^(-1)*SAMPLE*LAMBDA))^2
                }else{
                uu<-uu+((V[ij])^2/((D[ij]*u1)+((2)^(-1)*SAMPLE*LAMBDA))^2)
             }
          }
          return(-1+uu)
        }
        uvect<-uniroot.all(f2a, c(1e-5, 1e5))
        if(length(uvect)==0){u=1}else if(length(uvect)>0){u<-uvect[1]}
        BETA[,i4]<-solve(XlXl[[i4]]+(((SAMPLE*LAMBDA)/(2*u))*diag(1,dim(XlXl[[i4]])[1])))%*%W
        j<-sort(unique(c(j,i4)))
      }
    }
  }
  if( (norm(vec(BETA)-vec(BETATempo),type="1") <= PRECISION) || (KK==1100)){
      CONDITION1<-FALSE
   }
   cat("\f")
   cat("\n","Convergence=", norm(vec(BETA)-vec(BETATempo),type="1"),"  KK=", KK," j=", j, "\n")
   #View(t(BETA))
  KK<-KK+1
}
########################## END: STEP4A #####

### COLLECT CURRENT BETA AND ACTIVE SET
liSTSLS<-list()
liSTSLS[[1]]<-j
liSTSLS[[2]]<-BETA
liSTSLS[[3]]<-KK
return(liSTSLS)
}
########################################################

SamSize<-1200

for(REPLICATE in minREP:maxREP){
### SIMULATION STUDY BEGIN ########################################
cat("REPLICATION=", REPLICATE,"\n")
timeSTART<-Sys.time()
DATACOLL<-DATA3A(SamSize,d=1,r1= -2.0, r2= 2.0, SEED=REPLICATE)
DATA3<-DATACOLL$DATA3
Y0<-DATACOLL$Y1
p<-2
Q<-1
delay<-DATACOLL$delay
TBETA<-matrix(vec(DATACOLL$tCOFF),ncol= (3), byrow=F) #ncol is the number of regime
TRUETHR<-c(-2.0,2.0)

#############################################################################################
DATAC<-cbind(1,DATA3[,c(-1:(-Q-1))])                       #IGNORE THE COLUMN OF THRESHOLD VARIABLE AND ADD A INTERCEPT TERM
ZERO<-as.matrix(replace(DATAC,DATAC!=is.na(DATAC),0))      #CREATE A ZERO MATRIX WITH SIMILAR DIMENSIONS WITH "DATAC"
SAMPLE<-dim(DATAC)[1]
Y01<-DATA3[,c(2:(Q+1))];   Y0<-vec(t(Y01));                #Y0= (y_{\pi(1)+d},...,y_{\pi(n*)+d}}: ARRANGED DEPENDENT VARIABLES 
######################################################
####  DATA DEFINITION: SECOND STEP ESTIMATE ##########
core<-3                                                    #NO OF PROCESSOR CORES THAT WILL BE USED
DATABEA<-DATACOLL$SSTEP
Xo<-as.matrix(cbind(1,DATABEA[,c(-1:(-Q-1))]))             #UNORDERED AR VARIABLES
Yo<-as.matrix((DATABEA[,c(2:(Q+1))]))                      #UNORDERED ENDOGENOUS VARIABLE
Z<-as.vector(DATABEA[,1])                                  #UNORDERED THRESHOLD VARIABLE
########################################################

BETAlist<-AcTiSETlist<-list()
BETAlist2<-AcTiSETlist2<-list()
SSElist1S<-IClist1S<-numeric()
KKlist<-numeric()
#### INITIALIZING ################################
BETAlist[[1]]<-matrix(0, nrow=(p+1), ncol=SAMPLE)
AcTiSETlist[[1]]<-1
##################################################

ij<-1
bb<-1
while(ij<=LAMBDAlength){ 
if(ij == 1){
listSLS<-SLSalgo(LAMBDAlist[ij],AcTiSETlist[[1]],BETAlist[[1]]) 
}else{
listSLS<-SLSalgo(LAMBDAlist[ij],AcTiSETlist[[bb]],BETAlist[[bb]])   
}
### listSLS[[1]] is the active set, listSLS[[2]] is Est. Parameters
  SSElist1S[[bb]]<-SSE1S<-SSESETAR1S(listSLS[[1]],listSLS[[2]])
  IClist1S[bb]<-BICSETAR1S(SAMPLE,SSE1S,listSLS[[1]],p,0.01)
  AcTiSETlist2[[bb]]<-AcTiSETlist[[bb+1]]<-listSLS[[1]]
  BETAlist2[[bb]]<-BETAlist[[bb+1]]<-listSLS[[2]]
  KKlist[bb]<-listSLS[[3]]
  bb<-bb+1
  
  ij<-ij+1  
}

if(bb==1){
FINALbetaS1<-BETAlist[[1]]
FINALlocaS1<-AcTiSETlist[[1]]
FINALKKK<-KKlist[1]
  }else{
FINALbetaS1<-BETAlist2[[which(min(IClist1S)==IClist1S)[1]]]
FINALlocaS1<-AcTiSETlist2[[which(min(IClist1S)==IClist1S)[1]]]
FINALKKK<-KKlist[which(min(IClist1S)==IClist1S)[1]]
}

if(length(FINALlocaS1)==1){
  thresholdA<-matrix(0,nrow=0,ncol=4)
}else{
  for (J in 2:length(FINALlocaS1)){
    if(J ==2 ){
      thresholdA<-matrix(0,nrow=1,ncol=4)
      thresholdA[1,1]<-FINALlocaS1[J]
      thresholdA[1,2]<-DATA3[FINALlocaS1[J],1]
      thresholdA[1,3]<-FINALlocaS1[J]-1
      thresholdA[1,4]<-DATA3[FINALlocaS1[J]-1,1]
    }else{
      thresholdB<-matrix(0,nrow=1,ncol=4)
      thresholdB[1,1]<-FINALlocaS1[J]
      thresholdB[1,2]<-DATA3[FINALlocaS1[J],1]
      thresholdB[1,3]<-FINALlocaS1[J]-1
      thresholdB[1,4]<-DATA3[FINALlocaS1[J]-1,1]
      thresholdA<-rbind(thresholdA,thresholdB)
    }
  }
}

#Final : ESTIMATED THRESHOLDS                     # (LINE 182)
colnames(thresholdA) <- c("Location","Threshold","NeighLocat","NeighThre")
FINALlocaS1<-(sort(FINALlocaS1)[-1])-1                                          #TAKE A LAG OF INDICES
THRESHOLD<-DATA3[c(sort(FINALlocaS1)),1]  

### COLLECT RESULT: BCD ALGO
BCDTHRES<-THRESHOLD
BCDNTHRES<-length(THRESHOLD)
BCDHausL<-OneSHAUS(THRESHOLD,TRUETHR)
BCDHausR<-OneSHAUS(TRUETHR,THRESHOLD)
BCDHaus<-max(c(BCDHausL,BCDHausR))
#######################################################################################################
timeEND<-Sys.time()
TIMESIMU<-print(difftime(timeEND, timeSTART, units = "mins")[[1]])

#######################################################################################################

LISTsimu<-list(
  list(BCDTHRES,BCDNTHRES,BCDHausL,BCDHausR,BCDHaus),
  list(TIMESIMU, FINALKKK),
  list(FINALlocaS1, THRESHOLD, SamSize),
  list(TBETA,TRUETHR)
    )
LISTA<-list(assign(paste("LISTBCDBEA",REPLICATE,sep=""),LISTsimu))
save(LISTA,file=paste("4. SLS vs aBCD/1.SLS/750/2.Second/SIMUBCDBEA.(",REPLICATE,").RData",sep = ""))
}
################################################################################################################################