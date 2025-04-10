#REAL US GNP

DATA3A<-function(d,p){
#A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
library(mvtnorm)
library(MASS)
library(fBasics)

Q<- 1 #Number of Variate

par(mfrow=c(1,2))
YA <- read.csv("C:/Users/Jaffri/Dropbox/PhD Research/4. Third Year/3. Thesis Write up/Simulation Studies-Rcode/A. SETAR models/C. Real case Study/GNP2018.csv")
TST1<-ts(YA[,2],frequency = 4,start=c(1947,1))
Y1<-as.matrix(100*diff(log(YA[,2]), differences =1))
TST2<-ts(Y1,frequency = 4,start=c(1947,2))
#plot.ts(TST1,,main="Original US GNP")
#plot.ts(TST2,,main="Growth rate of US GNP")
plot(TST1,type="l",main="Original US GNP",xlab="",ylab="")
plot(TST2,type="l",main="Growth rate of US GNP",xlab="",ylab="")

library(cowplot)
library(readxl)
library(ggplot2)
library(forecast)
plotA<-ggAcf(TST2,lag.max = 60,lag=60, main="ACF of log return U.S GNP")
plotB<-ggPacf(TST2,lag.max = 60,lag=60, main="PACF of log return U.S GNP")
plot_grid(plotA, plotB,  labels = c('', '' ,'') , ncol = 2)

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

ALLDATA<-list(DATA3,Y1,d,p,SSTEP)
names(ALLDATA)<-c("DATA3","Y1","delay","p","SSTEP")
return(ALLDATA)
}
DATACOLL<-DATA3A(d=6,p=11)
DATA3<-DATACOLL$DATA3
Y0<-DATACOLL$Y1 
Q<-1
p<-DATACOLL$p
delay<-DATACOLL$delay
TBETA<-NA

library(TSA)
pvaluem=NULL
for (d in 1:11){
     res=tlrt(Y0,p=11,d=d,a=0.05,b=0.95)
     pvaluem= cbind( pvaluem, round(c(d,signif(c(res$test.statistic, 
                                                 res$p.value))),3))
 }
 rownames(pvaluem)=c('d','test statistic','p-value')
pvaluem

