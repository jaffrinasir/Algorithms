## THIS IS BACKWARD ELIMINATION ALGORITHM (BEA) FOR UNIVARIATE THRESHOLD AUTOREGRESSIVE MODEL (6.6.17)
## I WRITTEN CODES FOLLOWING STEP-BY-STEP IN MY RECENT DRAFTS OF ALGORITHM
## THE CODES ARE WRITTEN IN R AND C++ CODES

#### SET SYSTEM ENVIRONMENT (Cpp) ######################
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -std=c++11")
setwd("C:/Users/Jaffri/Dropbox/PhD Research/4. Third Year/3. Thesis Write up/Simulation Studies-Rcode/A. SETAR models/A. Coding development/GLARS")
########################################################

####  LOAD REQUIRED PACKAGES ###########################
library(Matrix)
library(fBasics)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(rootSolve)
########################################################

####  DATA DEFINITION  AND INITIALIZATIONS #############
core<-3                                                    #NO OF PROCESSOR CORES THAT WILL BE USED
DATABEA<-DATACOLL$SSTEP
Xo<-as.matrix(cbind(1,DATABEA[,c(-1:(-Q-1))]))             #UNORDERED AR VARIABLES
Yo<-as.matrix((DATABEA[,c(2:(Q+1))]))                      #UNORDERED ENDOGENOUS VARIABLE
Z<-as.vector(DATABEA[,1])                                  #UNORDERED THRESHOLD VARIABLE
########################################################

#### LOAD R FUNCTIONS  #################################
ICb<-function(jCard,No,SSEa,Cn){
    IC1b<-(No*log(SSEa/No))+(jCard*Cn*log(No))
  return(IC1b)
}
########################################################

#### LOAD C++ FUNCTIONS  ###############################
sourceCpp("./CppFiles/B1.COMPSSE.cpp")                    #FOR "COMPSSE" FUNCTIONS
sourceCpp("./CppFiles/B2.ARCOEF.cpp")                     #FOR "ARCOEF"  FUNCTION
########################################################

## COMPUTE INFORMATION CRITERION FOR INITIAL SET. ######
N<-dim(Yo)[1]                                              #NUMBER OF AVAILABLE OBSERVATIONS 
if(length(THRESHOLD)>0){
 THRESHOLD<-sort(unique(THRESHOLD))                        #REMOVE ANY SIMILAR VALUE OF THRESHOLD
}
## BEGIN : STEP 1-COMPUTE INFORMATION CRITERION FOR INITIAL SET
  o<-length(THRESHOLD)+1
  ICBEA<-c(rep(NA,o))
  THRESHOLDBEA<-sort(c(-Inf,Inf,THRESHOLD)) 
   for(i in 1: (length(THRESHOLDBEA) - 1)){
     if(i==1){
     LCSSE<-COMPSSE(Xo,Yo,Z,THRESHOLDBEA[i],THRESHOLDBEA[i+1],N)
     SSE<-(LCSSE$SSE)
     } else{
       LCSSE<-COMPSSE(Xo,Yo,Z,THRESHOLDBEA[i],THRESHOLDBEA[i+1],N)
       SSE<-SSE+(LCSSE$SSE)
      }
   }
  ICBEA[o]<-ICb(length(THRESHOLD),N,SSE,3)             #(m+1)-REGIME
## END : STEP 1 ######################################### 

#*** BEGIN - BACKWARD ELIMINATION ALGORITHM (BEA) *****#
if(length(THRESHOLD)!=0){
	COND<-TRUE
while(COND==TRUE){                                         #LOOP WILL CONTINUE UNTIL COND<-FALSE
print(o)
## STEP 2 : COMPUTE V_{m,i} #############################
V<-c(rep(NA,length(THRESHOLD)))	
for(j in 1:length(THRESHOLD))
{
THRESHOLDBEA<-sort(c(-Inf,Inf,THRESHOLD[-j])) 
for(i in 1: (length(THRESHOLDBEA) - 1)){
        if(i==1){
        LCSSE<-COMPSSE(Xo,Yo,Z,THRESHOLDBEA[i],THRESHOLDBEA[i+1],N)
        SSE<-LCSSE$SSE
        } else{
        LCSSE<-COMPSSE(Xo,Yo,Z,THRESHOLDBEA[i],THRESHOLDBEA[i+1],N)
        SSE<-SSE+LCSSE$SSE
     }
  }
  V[j]<-ICb((length(THRESHOLD)-1),N,SSE,3)                 #m-REGIME (REMOVE ONE THRESHOLD)
}
##SET V*_{m-1}=min_{i} V_{m,i}
ICBEA[o-1]<-min(V)
## END STEP 2 ###########################################

## STEP 3 - CHECK CONDITIONS ############################
   if(ICBEA[o-1]<ICBEA[o]){
    THRESHOLD<-THRESHOLD[-which(min(V)==V)[1]]                #CONDITION b.), REPEAT STEP 2-3
    o<-o-1
   }else if(ICBEA[o-1]>=ICBEA[o]){                         #CONDITION a.)
    COND<-FALSE
   }

   ##ADDITIONAL CONDITION FOR EMPTY THRESHOLD SET          #THIS SUPPORT CONDITION a.) b.), c.)
   if(length(THRESHOLD)==0){
     COND<-FALSE
   }
## END STEP 3 ###########################################
 }
}
#*** END   - BACKWARD ELIMINATION ALGORITHM (BEA) *****#

## COMPUTE AUTOREGRESSIVE COEFFICIENTS FOR EACH REGIME #
for(j in 1:length(THRESHOLD)+1)
{
  THRESHOLDCOEF<-sort(c(-Inf,Inf,THRESHOLD)) 
  for(i in 1: (length(THRESHOLDCOEF) - 1)){
    if(i==1){
      BETAa<-ARCOEF(Xo,Yo,Z,THRESHOLDCOEF[i],
        THRESHOLDCOEF[i+1],N)
      BETA<-BETAa$BETA
      COUNT<-BETAa$COUNT
      SSEC<-BETAa$SSE
    } else{
      BETAa<-ARCOEF(Xo,Yo,Z,THRESHOLDCOEF[i],
                    THRESHOLDCOEF[i+1],N)
      BETA<-cbind(BETA,BETAa$BETA)
      COUNT<-cbind(COUNT,BETAa$COUNT)
      SSEC<-cbind(SSEC,BETAa$SSE)
    }
  }
}
THRESHOLD
COUNT
BETA
SSEC