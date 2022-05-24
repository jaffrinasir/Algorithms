## THIS IS GROUP LARS FOR UNIVARIATE THRESHOLD AUTOREGRESSIVE MODEL (5.6.17)
## I WRITTEN R CODES FOLLOWING STEP-BY-STEP IN MY RECENT DRAFTS OF ALGORITHM

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
sourceCpp("./CppFolder/A1.CORRFUNC.cpp")                     #FOR "CORRFUNC" FUNCTIONS
sourceCpp("./CppFolder/A2.SOLNORM2.cpp")                       #FOR "SOLNORM2"  FUNCTION
########################################################

####  DATA DEFINITION  AND INITIALIZATIONS #############
core<-3                                                      #NO OF PROCESSOR CORES THAT WILL BE USED
DATAC<-cbind(1,DATA3[,c(-1:(-Q-1))])                         #IGNORE THE COLUMN OF THRESHOLD VARIABLE AND ADD A INTERCEPT TERM
ZERO<-as.matrix(replace(DATAC,DATAC!=is.na(DATAC),0))        #CREATE A ZERO MATRIX WITH SIMILAR DIMENSIONS WITH "DATAC"
Y01<-DATA3[,c(2:(Q+1))];   Y0<-vec(t(Y01));                  #Y0 is : Y^{0}_{n}    
UPSI<-Y0                                                     #UPSI is : \upsilon^{[0]} = Y^{0}_{n}, THE RESPONSE VARIABLE
ACset<-c(1:dim(DATAC)[1])                                    #DECLARE all (n-h) - GROUPS (INTEGERS)  
MU<-as.matrix(replace(Y0,Y0!=is.na(Y0),0))                   #MU is : \mu^[0] = 0  
MAXITER<-60                                                   #MAXIMUM ITERATIONS FOR GROUP LARS
DELTA<-0
k<-1                                                         #SET INITIAL INDEX
j<-numeric()                                                 # j is : \mathcal(A)^{[0]} = \emptyset, AN ACTIVE SET
########################################################

#*** BEGIN - GLARS (NORM-2) ***************************#
#### STEP 2: COMPUTE THE "MOST CORRELATED SET" #########
#### CREATE FUNCTION TO COMPUTE CORRELATION    #########
CORR<-CORRFUNCcpp(DATAC,Y0,ACset,core)
COVMAX<-max(CORR)
j<-ACset[which(COVMAX==CORR)][1];  j<-sort(j)            #FIRST ELEMENT IN ACTIVE SET
jT<-sort(ACset[! ACset %in% seq(from=j-DELTA,
                                to=j+DELTA,by=1)])                     
#UPDATE FIRST INACTIVE SET                    
#### END: STEP 2   #####################################

if(MAXITER==1){
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
if(k==MAXITER || ALPHA==1){
	ITER<- FALSE
}
k<- k+1
#### END: STEP 5     ####################################
}
#*** END   - GLARS (NORM-2) ***************************#


jFinal<-j

if(length(jFinal)==1){
  thresholdA<-matrix(0,nrow=0,ncol=4)
}else{
  for (J in 2:length(jFinal)){
    if(J ==2 ){
      thresholdA<-matrix(0,nrow=1,ncol=4)
      thresholdA[1,1]<-jFinal[J]
      thresholdA[1,2]<-DATA3[jFinal[J],1]
      thresholdA[1,3]<-jFinal[J]-1
      thresholdA[1,4]<-DATA3[jFinal[J]-1,1]
    }else{
      thresholdB<-matrix(0,nrow=1,ncol=4)
      thresholdB[1,1]<-jFinal[J]
      thresholdB[1,2]<-DATA3[jFinal[J],1]
      thresholdB[1,3]<-jFinal[J]-1
      thresholdB[1,4]<-DATA3[jFinal[J]-1,1]
      thresholdA<-rbind(thresholdA,thresholdB)
    }
  }
}

#Final : ESTIMATED THRESHOLDS
colnames(thresholdA) <- c("Location","Threshold","NeighLocat","NeighThre")
thresholdA
jFinal<-(sort(jFinal)[-1])-1      #TAKE A LAG OF INDICES
THRESCOLL<-THRESHOLD<-DATA3[c(sort(jFinal)),1]
abline(v = sort(j)[-1])