#EXAMPLE 3. IN GOGA PAPER (2016)
#AUTOREGRESSIVE ORDER OF 1 FOR EACH REGIME WITH DELAY PARAMETER (d=1)
#LOAD ANY REQUIREMENT

library(mvtnorm)
library(MASS)
library(fBasics)

DATA3A<-function(n,d,r1,r2,SEED){
#A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
p<-2      #AR(2)
n<-n+200  #Below obs 201 will be deleted (See line 37)
K<- 1 #Number of Variate

#############################################################################
 if(SEED== 'NA'){
 set.seed(sample(1:100000, 1))
} else{
 set.seed(SEED)	
}
Z<-rmvnorm(n, mean=c(0), diag(1)) #Generate r.v
Z<-t(Z)
#############################################################################

d<- d   #d is delay parameter
r1<- r1  #r1-r2 is threshold variable
r2<- r2

paraC1<-matrix(c(0),nrow=K,byrow=TRUE)
paraC2<-matrix(c(0),nrow=K,byrow=TRUE)
paraC3<-matrix(c(0),nrow=K,byrow=TRUE)
#AR parameter for FIRST regime 
paraT11<-matrix(c(0.8),nrow=K,byrow=TRUE)
paraT12<-matrix(c(-0.2),nrow=K,byrow=TRUE)
#AR parameter for SECOND regime 
paraT21<-matrix(c(1.9),nrow=K,byrow=TRUE)
paraT22<-matrix(c(-0.81),nrow=K,byrow=TRUE)
#AR parameter for THIRD regime 
paraT31<-matrix(c(0.6),nrow=K,byrow=TRUE)
paraT32<-matrix(c(-1.0),nrow=K,byrow=TRUE)

#HERE, SIGMA1-SIGMA3 and SD1-SD3 matrices can have different diagonal values for each regime
SIGMA1<-matrix(c(1),ncol=K,byrow=TRUE); SD1<-chol(SIGMA1); #Var-Cov First regime
SIGMA2<-matrix(c(1),ncol=K,byrow=TRUE); SD2<-chol(SIGMA2); #Var-Cov Second regime
SIGMA3<-matrix(c(1),ncol=K,byrow=TRUE); SD3<-chol(SIGMA3); #Var-Cov Second regime
sigmaY<-matrix(c(1),ncol=K,byrow=TRUE) 

Y <- t(rmvnorm(n=p, mean=c(0), sigma=sigmaY))   #GENERATE INITAL SERIES OF Y, #MUST TRANSPOSE
ZERO<-matrix(0,nrow=K,ncol=(dim(Z)[2]-dim(Y)[2]))
Y<-cbind(Y,ZERO); 

for (i in (p+1):((dim(Z)[2]))) { 
	if(Y[,i-d]<=r1){                                       
	Y[,i]<-paraC1+(paraT11%*%(as.matrix(Y[,i-1])))+(paraT12%*%(as.matrix(Y[,i-2])))+(SD1%*%Z[,i])	
	}
	else if(Y[,i-d]> r1 && Y[,i-d]<=r2){
	Y[,i]<-paraC2+(paraT21%*%(as.matrix(Y[,i-1])))+(paraT22%*%(as.matrix(Y[,i-2])))+(SD2%*%Z[,i])	
	}else{
	Y[,i]<-paraC3+(paraT31%*%(as.matrix(Y[,i-1])))+(paraT32%*%(as.matrix(Y[,i-2])))+(SD3%*%Z[,i])	
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
n=1200
DATACOLL<-DATA3A(n,d=1,r1= -2, r2= 2,SEED="NA")
DATA3<-DATACOLL$DATA3
Y0<-DATACOLL$Y1
p<-2 
Q<-1
delay<-DATACOLL$delay
TBETA<-matrix(vec(DATACOLL$tCOFF),ncol= (3), byrow=F) #ncol is the number of regime

### PLOTING TIME SERIES #################################################################
par(mfrow=c(1,2))
plot(NA, xlim=c(1,n), ylim=c(-7,7), main="Original Series", ylab="", xlab="Time")
points(Y0[,1],col="red", type="l", pch=23,lty="dotted", cex=4)

plot(NA, xlim=c(1,n), ylim=c(-7,7), main="Rearranged observations", ylab="", xlab="Time")
points(DATA3[,2],col="red", type="l", pch=23,lty="dotted", cex=4)
#legend(n-(n/2 -70), 8.1, legend=c("First Series"),
#       col=c("red"), lty=c(3,1), cex=1)
#########################################################################################