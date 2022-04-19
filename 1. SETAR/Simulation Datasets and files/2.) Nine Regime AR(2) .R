library(mvtnorm)
library(MASS)
library(fBasics)

DATA3A<-function(n,d,r1,r2,r3,r4,r5,r6,r7,r8,SEED,C1,C2,C3,C4,C5,C6,C7,C8,C9){
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

paraC1<-matrix(c(C1),nrow=K,byrow=TRUE)
paraC2<-matrix(c(C2),nrow=K,byrow=TRUE)
paraC3<-matrix(c(C3),nrow=K,byrow=TRUE)
paraC4<-matrix(c(C4),nrow=K,byrow=TRUE)
paraC5<-matrix(c(C5),nrow=K,byrow=TRUE)
paraC6<-matrix(c(C6),nrow=K,byrow=TRUE)
paraC7<-matrix(c(C7),nrow=K,byrow=TRUE)
paraC8<-matrix(c(C8),nrow=K,byrow=TRUE)
paraC9<-matrix(c(C9),nrow=K,byrow=TRUE)
#AR(2) parameter for FIRST regime 
paraT11<-matrix(c(-0.6),nrow=K,byrow=TRUE)
paraT12<-matrix(c(0),nrow=K,byrow=TRUE)
#AR(2) parameter for SECOND regime 
paraT21<-matrix(c(0.3),nrow=K,byrow=TRUE)
paraT22<-matrix(c(0.9),nrow=K,byrow=TRUE)
#AR(2) parameter for THIRD regime 
paraT31<-matrix(c(-0.9),nrow=K,byrow=TRUE)
paraT32<-matrix(c(0),nrow=K,byrow=TRUE)
#AR(2) parameter for FORTH regime 
paraT41<-matrix(c(0.7),nrow=K,byrow=TRUE)
paraT42<-matrix(c(0.5),nrow=K,byrow=TRUE)
#AR(2) parameter for FIFTH regime 
paraT51<-matrix(c(0.1),nrow=K,byrow=TRUE)
paraT52<-matrix(c(0),nrow=K,byrow=TRUE)
#AR(2) parameter for SIXTH regime 
paraT61<-matrix(c(-0.9),nrow=K,byrow=TRUE)
paraT62<-matrix(c(0),nrow=K,byrow=TRUE)
#AR(2) parameter for SEVENTH regime 
paraT71<-matrix(c(0.9),nrow=K,byrow=TRUE)
paraT72<-matrix(c(0),nrow=K,byrow=TRUE)
#AR(2) parameter for EIGHTH regime 
paraT81<-matrix(c(-0.8),nrow=K,byrow=TRUE)
paraT82<-matrix(c(-0.2),nrow=K,byrow=TRUE)
#AR(2) parameter for EIGHTH regime 
paraT91<-matrix(c(-1.1),nrow=K,byrow=TRUE)
paraT92<-matrix(c(0),nrow=K,byrow=TRUE)

#HERE, SIGMA1-SIGMA3 and SD1-SD3 matrices can have different diagonal values for each regime
SIGMA1<-matrix(c(1),ncol=K,byrow=TRUE); SD1<-chol(SIGMA1); #Var-Cov First regime
SIGMA2<-matrix(c(1),ncol=K,byrow=TRUE); SD2<-chol(SIGMA2); #Var-Cov Second regime
SIGMA3<-matrix(c(1),ncol=K,byrow=TRUE); SD3<-chol(SIGMA3); #Var-Cov Second regime
SIGMA4<-matrix(c(1),ncol=K,byrow=TRUE); SD4<-chol(SIGMA4); #Var-Cov First regime
SIGMA5<-matrix(c(1),ncol=K,byrow=TRUE); SD5<-chol(SIGMA5); #Var-Cov Second regime
SIGMA6<-matrix(c(1),ncol=K,byrow=TRUE); SD6<-chol(SIGMA6); #Var-Cov Second regime
SIGMA7<-matrix(c(1),ncol=K,byrow=TRUE); SD7<-chol(SIGMA7); #Var-Cov First regime
SIGMA8<-matrix(c(1),ncol=K,byrow=TRUE); SD8<-chol(SIGMA8); #Var-Cov Second regime
SIGMA9<-matrix(c(1),ncol=K,byrow=TRUE); SD9<-chol(SIGMA9); #Var-Cov Second regime
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
  }  
  else if(Y[,i-d]> r2 && Y[,i-d]<=r3){
  Y[,i]<-paraC3+(paraT31%*%(as.matrix(Y[,i-1])))+(paraT32%*%(as.matrix(Y[,i-2])))+(SD3%*%Z[,i])  
  }
  else if(Y[,i-d]> r3 && Y[,i-d]<=r4){
  Y[,i]<-paraC4+(paraT41%*%(as.matrix(Y[,i-1])))+(paraT42%*%(as.matrix(Y[,i-2])))+(SD4%*%Z[,i])  
  }
  else if(Y[,i-d]> r4 && Y[,i-d]<=r5){
  Y[,i]<-paraC5+(paraT51%*%(as.matrix(Y[,i-1])))+(paraT52%*%(as.matrix(Y[,i-2])))+(SD5%*%Z[,i])  
  }
  else if(Y[,i-d]> r5 && Y[,i-d]<=r6){
  Y[,i]<-paraC6+(paraT61%*%(as.matrix(Y[,i-1])))+(paraT62%*%(as.matrix(Y[,i-2])))+(SD6%*%Z[,i])  
  }
  else if(Y[,i-d]> r6 && Y[,i-d]<=r7){
  Y[,i]<-paraC7+(paraT71%*%(as.matrix(Y[,i-1])))+(paraT72%*%(as.matrix(Y[,i-2])))+(SD7%*%Z[,i])  
  }
  else if(Y[,i-d]> r7 && Y[,i-d]<=r8){
  Y[,i]<-paraC8+(paraT81%*%(as.matrix(Y[,i-1])))+(paraT82%*%(as.matrix(Y[,i-2])))+(SD8%*%Z[,i])  
  }
  else{
  Y[,i]<-paraC9+(paraT91%*%(as.matrix(Y[,i-1])))+(paraT92%*%(as.matrix(Y[,i-2])))+(SD9%*%Z[,i])  
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

tCOFF<-cbind(paraC1,paraT11,paraT12,paraC2,paraT21,paraT22,paraC3,paraT31,paraT32,paraC4,paraT41,paraT42,
	paraC5,paraT51,paraT52,paraC6,paraT61,paraT62,paraC7,paraT71,paraT72,paraC8,paraT81,paraT82,paraC9,paraT91,paraT92)

ALLDATA<-list(DATA3,Y1,d,tCOFF,SSTEP)
names(ALLDATA)<-c("DATA3","Y1","delay","tCOFF","SSTEP")
return(ALLDATA)
}
n=5000
##REPLACE COEF WITH
# -4.5, 2.5, -2.0, 2.3, 1.0, 3.0, 1.6, -0.5, 1.5
# 2.0, 3.0, 4.0, 9.0, 8.0, 11.0, 9.0, 12.0, 9.0
# -0.6, 1.6, -0.6, 1.6, -0.6, 1.6, -0.6, 1.6, -0.6

DATACOLL<-DATA3A(n,d=1,r1= -3.5, r2= -2.5 ,r3= -1.5, r4= -0.5,r5= 0.5, r6=1.5, r7=2.5,r8=3.5,SEED="NA",
                 -4.5, 2.5, -2.0, 2.3, 1.0, 3.0, 1.6, -0.5, 1.5
	)
DATA3<-DATACOLL$DATA3
Y0<-DATACOLL$Y1
p<-2
Q<-1
delay<-DATACOLL$delay
TBETA<-matrix(vec(DATACOLL$tCOFF),ncol= (9), byrow=F) #ncol is the number of regime

### PLOTING TIME SERIES #################################################################
#par(mfrow=c(1,2))
#plot(NA, xlim=c(1,n), ylim=c(-8,10), main="First Series", ylab="", xlab="Time")
#points(Y0[,1],col="red", type="l", pch=23,lty="dotted", cex=4)

plot(NA, xlim=c(1,n), ylim=c(-8,10), main="Rearranged observations", ylab="", xlab="Time")
points(DATA3[,2],col="red", type="l", pch=23,lty="dotted", cex=4)
#legend(n-(n/2 -70), 8.1, legend=c("First Series"),
#       col=c("red"), lty=c(3,1), cex=1)
#########################################################################################
