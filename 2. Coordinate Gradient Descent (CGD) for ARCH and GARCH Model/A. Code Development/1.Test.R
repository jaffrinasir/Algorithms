library(fGarch)

###################### CASE 1 ########################################################################################################
spec1 = garchSpec(model = list(omega = 1, alpha = c(0.1), beta = c(0.1)))  ##
for(i in 1:1000){
GFIT<-0
ERRM1<-garchSim(spec1, n = 5000,n.start=500)
GFIT1<-garchFit(data = ERRM1$garch, cond.dist = "QMLE", include.mean = FALSE,formula =~garch(1,1) , description = NULL,trace = F)@fit$par
if(i==1){
COLLE1<-GFIT1
GBIAS1<-GFIT1-c(1,0.1,0.1)
GVAR1<-(GFIT1-c(1,0.1,0.1))^2
}else{
  COLLE1<-cbind(COLLE1,GFIT1)
  GBIAS1<-cbind(GBIAS1,(GFIT1-c(1,0.1,0.1)))
  GVAR1<-cbind(GVAR1,(GFIT1-c(1,0.1,0.1))^2)
  }
}
par(mfrow=c(3,1))
plot(COLLE[1,],type="l",xlab="Iterations",ylab=expression(omega))
plot(COLLE[2,],type="l",xlab="Iterations",ylab=expression(alpha))
plot(COLLE[3,],type="l",xlab="Iterations",ylab=expression(beta))
MEAN1<-rowSums(COLLE1)/1000
BIAS1<-rowSums(GBIAS1)/1000
ESD1<-sqrt(rowSums(GVAR1)/1000)

##################### CASE 2 #########################################################################################################
spec2 = garchSpec(model = list(omega = 0.01, alpha = c(0.2), beta = c(0.75)))  ##TAKEN FROM CHAN & McAleer (2002)
for(i in 1:1000){
  GFIT2<-0
  ERRM2<-garchSim(spec2, n = 5000,n.start=500)
  GFIT2<-garchFit(data = ERRM2$garch, cond.dist = "QMLE", include.mean = FALSE,formula =~garch(1,1) , description = NULL,trace = F)@fit$par
  if(i==1){
    COLLE2<-GFIT2
    GBIAS2<-GFIT2-c(0.01,0.2,0.75)
    GVAR2<-(GFIT2-c(0.01,0.2,0.75))^2
  }else{
    COLLE2<-cbind(COLLE2,GFIT2)
    GBIAS2<-cbind(GBIAS2,(GFIT2-c(0.01,0.2,0.75)))
    GVAR2<-cbind(GVAR2,(GFIT2-c(0.01,0.2,0.75))^2)
  }
}
par(mfrow=c(3,1))
plot(COLLE2[1,],type="l",xlab="Iterations",ylab=expression(omega))
plot(COLLE2[2,],type="l",xlab="Iterations",ylab=expression(alpha))
plot(COLLE2[3,],type="l",xlab="Iterations",ylab=expression(beta))
MEAN2<-rowSums(COLLE2)/1000
BIAS2<-rowSums(GBIAS2)/1000
ESD2<-sqrt(rowSums(GVAR2)/1000)