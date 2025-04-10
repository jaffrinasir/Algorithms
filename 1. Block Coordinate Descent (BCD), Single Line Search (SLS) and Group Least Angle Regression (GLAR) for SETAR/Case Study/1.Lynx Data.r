#LYNX DATA

DATA3A<-function(d,p){
#A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
library(mvtnorm)
library(MASS)
library(fBasics)

Q<- 1 #Number of Variate

Y1<-as.matrix(log10(lynx))
#Y1<-as.matrix(lynx)
par(mfrow=c(1,2))
TST1<-ts(lynx,frequency = 1,start=c(1821,1))
TST2<-ts(Y1,frequency = 1,start=c(1821,1))
plot(TST1,type="l",main="Original Lynx data",xlab="",ylab="")
plot(TST2,type="l",main="Transformed Lynx data",xlab="",ylab="")

library(cowplot)
library(readxl)
library(ggplot2)
library(forecast)
plotA<-ggAcf(TST2,lag.max = 60,lag=60, main="ACF of log transformed Canadian Lynx")
plotB<-ggPacf(TST2,lag.max = 60,lag=60, main="PACF of log transformed Canadian Lynx")
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
DATACOLL<-DATA3A(d=3,p=8)
DATA3<-DATACOLL$DATA3
Y0<-DATACOLL$Y1 
Q<-1
p<-DATACOLL$p
delay<-DATACOLL$delay
TBETA<-NA

library(TSA)
pvaluem=NULL
for (d in 1:5){
     res=tlrt(Y0,p=2,d=d,a=0.05,b=0.95)
     pvaluem= cbind( pvaluem, round(c(d,signif(c(res$test.statistic, 
                                                 res$p.value))),3))
 }
 rownames(pvaluem)=c('d','test statistic','p-value')
pvaluem

