plot(TST3,type="l",main="Rearranged Series + GLARS",xlab="",ylab="")
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
########################################################

#### LOAD C++ FUNCTIONS  ###############################
sourceCpp("./CppFolder/A1.CORRFUNC.cpp")                  #FOR "CORRFUNC" FUNCTIONS
sourceCpp("./CppFolder/A2.SOLNORM2.cpp")                  #FOR "SOLNORM2"  FUNCTION
sourceCpp("./CppFolder/B1.COMPSSE.cpp")                    #FOR "COMPSSE" FUNCTIONS
sourceCpp("./CppFolder/B2.ARCOEF.cpp")                     #FOR "ARCOEF"  FUNCTION
########################################################

#### LOAD R FUNCTIONS  #################################
ICb<-function(jCard,No,SSEa,Cn){
    IC1b<-(No*log(SSEa/No))+(jCard*Cn*log(No))
  return(IC1b)
}
########################################################

####  DATA DEFINITION  AND INITIALIZATIONS (GLARS) #####
core<-2                                                      #NO OF PROCESSOR CORES THAT WILL BE USED
DATAC<-cbind(1,DATA3[,c(-1:(-Q-1))])                         #IGNORE THE COLUMN OF THRESHOLD VARIABLE AND ADD A INTERCEPT TERM
ZERO<-as.matrix(replace(DATAC,DATAC!=is.na(DATAC),0))        #CREATE A ZERO MATRIX WITH SIMILAR DIMENSIONS WITH "DATAC"
Y01<-DATA3[,c(2:(Q+1))];   Y0<-vec(t(Y01));                  #Y0 is : Y^{0}_{n}    
UPSI<-Y0                                                     #UPSI is : \upsilon^{[0]} = Y^{0}_{n}, THE RESPONSE VARIABLE
ACset<-c(1:dim(DATAC)[1])                                    #DECLARE all (n-h) - GROUPS (INTEGERS)  
MU<-as.matrix(replace(Y0,Y0!=is.na(Y0),0))                   #MU is : \mu^[0] = 0  
Kmax<-10                                                      #MAXIMUM ITERATIONS FOR GROUP LARS
DELTA<-3
########################################################

####  DATA DEFINITION  AND INITIALIZATIONS (BEA) ########
DATABEA<-DATACOLL$SSTEP
Xo<-as.matrix(cbind(1,DATABEA[,c(-1:(-Q-1))]))             #UNORDERED AR VARIABLES
Yo<-as.matrix((DATABEA[,c(2:(Q+1))]))                      #UNORDERED ENDOGENOUS VARIABLE
Z<-as.vector(DATABEA[,1])                                  #UNORDERED THRESHOLD VARIABLE
########################################################

#*** BEGIN - GLARS (NORM-2) ***************************#
j<-numeric()                                                 # j is : \mathcal(A)^{[0]} = \emptyset, AN ACTIVE SET
k<-1                                                         #SET INITIAL INDEX
#### STEP 2: COMPUTE THE "MOST CORRELATED SET" #########
#### CREATE FUNCTION TO COMPUTE CORRELATION    #########
CORR<-CORRFUNCcpp(DATAC,Y0,ACset,core)
COVMAX<-max(CORR)
j<-ACset[which(COVMAX==CORR)][1];  j<-sort(j)            #FIRST ELEMENT IN ACTIVE SET
jT<-sort(ACset[! ACset %in% seq(from=j-DELTA,
                                to=j+DELTA,by=1)])                     
#UPDATE FIRST INACTIVE SET                    
#### END: STEP 2   #####################################

if(Kmax==1){
ITER <- FALSE 	
}else{
ITER <- TRUE 	
}

while(ITER == TRUE){                                          #BEGIN LOOP
#### STEP 3 : COMPUTE DESCENT DIRECTION ################
for(l in 1:length(j)){	
   if(l == 1){
      Xa<-rbind(ZERO[c(0:(j[l]-1)),],DATAC[c((j[l]):dim(DATAC)[1]),])
    } else{
      Xa<-cbind(Xa,rbind(ZERO[c(0:(j[l]-1)),],DATAC[c((j[l]):dim(DATAC)[1]),]))
    }
} 

M1<-t(Xa)%*%Xa 
M2<-t(Xa)%*%UPSI
BETA<-ginv(M1)%*%(M2)
#### END: STEP 3   #####################################
C0<-Xa%*%BETA                                                  #COMPUTE C^{[k]}
#### STEP 4 : COMPUTE a_{j}, b_{j} and c_{j} ###########
LISTSOL<-SOLNORM2(DATAC, UPSI, C0, jT, core) 
aj<-as.vector(LISTSOL$aj)
bj<-as.vector(LISTSOL$bj)
cj<-as.vector(LISTSOL$cj)
d<-max(CORRFUNCcpp(DATAC,UPSI,ACset,core))^2
ALPHAj<-numeric()

#CHECK CONDITIONS
for(ij in 1:length(jT)){
	SS<-(bj[ij]-d)^2-(aj[ij]-d)*(cj[ij]-d)
	if(SS<0){
		ALPHAj[ij]<-1
		}else if((SS>=0) && ((cj[ij]-d)!=0)){
			aPlus<-((bj[ij]-d)+sqrt(SS))/(cj[ij]-d)
            aMinus<-((bj[ij]-d)-sqrt(SS))/(cj[ij]-d)
            if(((aMinus>= 0) && (aMinus<=1))==TRUE){
               ALPHAj[ij]<-aMinus
            }else{ALPHAj[ij]<-aPlus}
		}else if((SS>=0) && ((cj[ij]-d)==0)){
			   ALPHAj[ij]<-(aj[ij]-d)/(2*(bj[ij]-d))
		  }
}
ALPHA<-min(ALPHAj)
#### END: STEP 4     ####################################

#### STEP 5 : UPDATING INFOMATION #######################
j<-sort(c(j,jT[which(ALPHA==ALPHAj)][1]))                         #UPDATE ACTIVE SET
jT<-sort(jT[! jT %in% c(j,seq(from=jT[which(ALPHA==ALPHAj)][1]-DELTA,
    to=jT[which(ALPHA==ALPHAj)][1]+DELTA,by=1))])                     #UPDATE INACTIVE SET
MU<-MU+(ALPHA*C0)
UPSI<-Y0-MU
if(k==Kmax || ALPHA==1){
	ITER<- FALSE
}
k<- k+1
#### END: STEP 5     ####################################
}
#*** END   - GLARS (NORM-2) ***************************#
jFinal<-j
PrioTHR<-DATA3[c(jFinal),1]; print(PrioTHR)
jFinal<-(sort(jFinal)[-1])-1      #TAKE A LAG OF INDICES
abline(v = sort(jFinal))
THRESHOLD<-DATA3[c(sort(jFinal)),1]

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
  ICBEA[o]<-ICb(length(THRESHOLD),N,SSE,5)             #(m+1)-REGIME
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
  V[j]<-ICb((length(THRESHOLD)-1),N,SSE,5)                 #m-REGIME (REMOVE ONE THRESHOLD)
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
BETA
COUNT

### COMPUTE RESIDUALS
#R CODES FOR FITTING TAR MODEL
THRESHOLDFIT<-sort(c(-Inf,Inf,THRESHOLD))
onep<-rep(1,dim(Xo)[2])
Yh<-0
COUNTONE<-0
Z1<-as.matrix(Z)

for(i in 1:(length(THRESHOLD)+1)){
INDI2<-numeric()
for(l in 1:length(Z1[,1])){
  INDI2[l]<-ifelse((
    THRESHOLDFIT[i]<Z1[l,1] && Z1[l,1]<=THRESHOLDFIT[i+1] #PUT PROPER STATEMENT HERE
  ),1,0)  
}
INDI<-matrix(INDI2); COUNTONE<-COUNTONE+INDI
Xj<-Xo*kron(t(onep),INDI)
Yh<-Yh+(Xj%*%matrix(BETA[,i]))
if(i==1){MATONE<-INDI}else{MATONE<-cbind(MATONE,INDI)}
}

library(TSA)
ERRY0<-Y0-Yh
Box.test(ERRY0, 12)
McLeod.Li.test(y=ERRY0, gof.lag=12)