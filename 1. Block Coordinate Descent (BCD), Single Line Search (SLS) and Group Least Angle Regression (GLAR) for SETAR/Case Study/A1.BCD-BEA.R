par(mfrow=c(2,1), mai = c(0.5, 0.5, 0.5, 0.5))
TST3<-ts(DATA3[,2],frequency = 1,start=c(1,1))
plot(TST3,type="l",main="Rearranged Series + BCD",xlab="",ylab="")
#### SET SYSTEM ENVIRONMENT (Cpp) ######################
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -std=c++11")
setwd("C:/Users/Jaffri/Dropbox/PhD Research/4. Third Year/3. Thesis Write up/Simulation Studies-Rcode/A. SETAR models/A. Coding development")
#Set environment first
########################################################

####  LOAD REQUIRED PACKAGES ###########################
library(Matrix)
library(fBasics)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(rootSolve)
########################################################
sourceCpp("./CppFiles/1.YBySGRAM.cpp") 
# This will produce SGRAM() and YBy() C++ function
sourceCpp("./CppFiles/B1.COMPSSE.cpp")                    #FOR "COMPSSE" FUNCTIONS
sourceCpp("./CppFiles/B2.ARCOEF.cpp")                     #FOR "ARCOEF"  FUNCTION
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
Kmax<-10+1
DELTA<-10

LAMBDAlength<-20
LAMBDAlist<-seq(0.5,0.01,length.out=LAMBDAlength)  #GENERATE LIST OF LAMBDAS

BETAlist<-AcTiSETlist<-list()
SSElist1S<-IClist1S<-numeric()
#### INITIALIZING ################################
BETAlist[[1]]<-matrix(0, nrow=(p+1), ncol=SAMPLE)
AcTiSETlist[[1]]<-1
##################################################

#### LOAD R FUNCTIONS  #################################
ICb<-function(jCard,No,SSEa,Cn){
    IC1b<-(No*log(SSEa/No))+(jCard*Cn*log(No))
  return(IC1b)
}

BICSETAR1S<-function(SAMPLEa,SSEa,Aset,Pa,Cn){   #BIC FOR SETAR, FIRST STAGE ESTIMATE
BIC1a<-((SAMPLEa)*log(SSEa/SAMPLEa))+(length(Aset)*Cn*log(SAMPLEa))
return(BIC1a)
}

SECANT<- function(f, x0, x1, tol = 1e-9, n = 1000) {
  for (i in 1:n) {
    x2 <- x1 - f(x1) / ((f(x1) - f(x0)) / (x1 - x0)) # Calculate the new x value
    if (abs(x2 - x1) < tol) { # If the difference between the new value and the previous value is small enough, end iteration and output root.
      return(x2)
    }
    # If the root was not determined in the previous iteration, update the values and proceed to the next iteration.
    x0 <- x1
    x1 <- x2 
  }
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

BCDalgo<-function(LAMBDA,ActiVE,BETAa){
PRECISION<-1e-3                                            #PRESSION CONSTANT (smaller constant produce slower convergence)
BETA<-BETAa                                                #INITIALIZATION FOR BETA
j<-sort(ActiVE)
jA<-c(1:dim(DATAC)[1]); 
REMOVE<-numeric()
for(jj in 1:length(j)){REMOVE<-c(REMOVE,c((j[jj]-DELTA):(j[jj]+DELTA)))}
jA<-sort(jA[! jA %in% unique(REMOVE)])

Fj<-function(INDEX,j,B){
  if(length(j)<=1){AC<-0}else{
      jB<-j[! j %in% INDEX]
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

K<-1
CONDITION2<-TRUE

if(length(j)>0){
while(CONDITION2==TRUE){
CONDITION1<-TRUE
KK<-0
while(CONDITION1==TRUE){
      BETATempo<-BETA
  for(i4 in 1:length(j)){
    W<-Xlyl[[j[i4]]]-Fj(j[i4],j,BETA)
    V<-vl[[j[i4]]]%*%W
    D<-dl[[j[i4]]]
    if(j[i4]==1){
      BETA[,j[i4]]<-solve(XlXl[[j[i4]]])%*%W
      }else{
    if((2*norm(W,type="2")/SAMPLE)<LAMBDA){BETA[,j[i4]]<-0}else{
    #f2a<-function(u1){
    #         return(- 1 +
    #          sum((V)^2/((dl[[j[i4]]]*u1)+((2)^(-1)*SAMPLE*LAMBDA))^2))
    #     }
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
          #u<-uniroot(f2a, c(0.0001, 100))$root
          #u<-SECANT(f2a, 0.0001, 100)
          #print(u)
          #cat("u = ", u, "f(u) =", f2a(u), "\n")
          BETA[,j[i4]]<-solve(XlXl[[j[i4]]]+(((SAMPLE*LAMBDA)/(2*u))*diag(1,dim(XlXl[[j[i4]]])[1])))%*%W
    }
  }
}
  if( (norm(vec(BETA)-vec(BETATempo),type="1") <= PRECISION) || (KK==9999) 
  ){CONDITION1<-FALSE
  }
  KK<-KK+1
}
########################### STEP3A #######
j0<-numeric()
for(i5 in 1:length(j)){
  if(norm(BETA[,j[i5]],type="2")==0){
    j0<-c(j0,j[i5])
  }
}
j<-sort(j[! j %in% j0])
######################### END: STEP3A #####

########################### STEP4A ########
if(length(jA)>0){
CORR2<-numeric()
for(i6 in 1:length(jA)){
   CORR2[i6]<-2*(norm(Xlyl[[jA[i6]]]-Fj(jA[i6],j,BETA), type="2"))/SAMPLE
} 

if(max(CORR2) >= LAMBDA){
  j<-sort(c(j, jA[which(max(CORR2)==CORR2)]))
  jA<-sort(jA[! jA %in% c((jA[which(max(CORR2)==CORR2)]-DELTA):
                            (jA[which(max(CORR2)==CORR2)]+DELTA))
           ])
  K<-K+1
}else{CONDITION2<-FALSE}
}else{CONDITION2<-FALSE}
########################## END: STEP4A #####
if((length(j)>=Kmax) || (length(jA)==0)){CONDITION2<-FALSE}
 }
}
### COLLECT CURRENT BETA AND ACTIVE SET
liSTBCD<-list()
liSTBCD[[1]]<-j
liSTBCD[[2]]<-BETA
return(liSTBCD)
}
########################################################

ITERBCD<-TRUE
ij<-1
while( (ITERBCD==TRUE) && (ij<=LAMBDAlength)){
cat("\r", "LAMBDA iteration = ",  ij, "\n")  
if(ij == 1){
listBCD<-BCDalgo(LAMBDAlist[ij],AcTiSETlist[[1]],BETAlist[[1]])	
}else{
listBCD<-BCDalgo(LAMBDAlist[ij],AcTiSETlist[[ij-1]],BETAlist[[ij-1]])		
}
### listBCD[[1]] is the active set, listBCD[[2]] is Est. Parameters
if(length(listBCD[[1]])< Kmax){
  SSElist1S[[ij]]<-SSE1S<-SSESETAR1S(listBCD[[1]],listBCD[[2]])
  AcTiSETlist[[ij]]<-listBCD[[1]]
  BETAlist[[ij]]<-listBCD[[2]]
  IClist1S[ij]<-BICSETAR1S(SAMPLE,SSE1S,listBCD[[1]],p,0.001)
  ij<-ij+1
	}else{
	ITERBCD<-FALSE}
}

FINALbetaS1<-BETAlist[[which(min(IClist1S)==IClist1S)[1]]]
FINALlocaS1<-AcTiSETlist[[which(min(IClist1S)==IClist1S)[1]]]

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
thresholdA
PrioTHR<-DATA3[c(FINALlocaS1),1]; print(PrioTHR)
FINALlocaS1<-(sort(FINALlocaS1)[-1])-1                                          #TAKE A LAG OF INDICES
THRESHOLD<-DATA3[c(sort(FINALlocaS1)),1]  
abline(v = sort(FINALlocaS1))


## COMPUTE INFORMATION CRITERION FOR INITIAL SET. #####################################################
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
ERRY0<-Yo-Yh
Box.test(ERRY0, 12,type = c("Box-Pierce", "Ljung-Box"))
McLeod.Li.test(y=ERRY0, gof.lag=12)
plot(Yo,type="l",col="red",main="Original vs Predicted",ylab=expression('y'['t']), xlab="Time Index (t)")
points(Yh,type="l",col="blue",lty="dashed")
legend("topright", inset=.05,
       c("Original","Predicted"),lty=c(1,2), col=c("red","blue"),horiz=F)

par(mfrow=c(1,2))
acf(ERRY0, main=expression(paste('ACF for ', epsilon['t'])))
acf(ERRY0^2, main=expression(paste('ACF for ', epsilon['t']^2)))