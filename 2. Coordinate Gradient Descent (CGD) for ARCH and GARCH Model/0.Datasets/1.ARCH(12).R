### LOAD DATASET ARCH(12) ####################################################
ARCH<-function(n){
  #A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
  library(mvtnorm)
  library(MASS)
  library(fBasics)
  library(matrixcalc)
  rCh<- 12
  n<-n+200  #Obs. below 201 will be deleted 
  Q<- 1 #Number of Variate
  
  #######  Generate r.v of white noise ########################################
  eta<-rmvnorm(n, mean=c(0), diag(1)) 
  eta<-t(eta)
  #############################################################################
  
  ##INTERCEPT
  paraGA0<-matrix(c(0.01),nrow=Q,byrow=TRUE)
  #### MA(p) representative
  paraGAa1<-matrix(c(0.15),nrow=Q,byrow=TRUE)
  paraGAa2<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa3<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa4<-matrix(c(0.3),nrow=Q,byrow=TRUE)
  paraGAa5<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa6<-matrix(c(0.2),nrow=Q,byrow=TRUE)
  paraGAa7<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa8<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa9<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa10<-matrix(c(0.15),nrow=Q,byrow=TRUE)
  paraGAa11<-matrix(c(0),nrow=Q,byrow=TRUE)
  paraGAa12<-matrix(c(0.19),nrow=Q,byrow=TRUE)

  H<-max(rCh)
  ht<-abs(t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1)))))  #HETEROS. OBS
  err<-t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1)))) #ERR. TERM
  Y <-t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1))))
  ZERO<-matrix(0,nrow=Q,ncol=(dim(eta)[2]-dim(Y)[2]))
  ht<-cbind(ht,ZERO) 
  err<-cbind(err,ZERO)
  Y<-cbind(Y,ZERO) 

ARCHr<-sapply( paste('paraGAa', 1:rCh, sep=''), get ,envir=sys.frame(sys.parent(0)))

    for (i in (H+1):((dim(eta)[2]))) {
      sumRg<-0
       for(j in 1:rCh){
        sumRg<-sumRg + (ARCHr[[j]]%*%(err[,(i-j)]^{2}))
       }
    ht[,i]<-paraGA0+(sumRg)
    err[,i]<-eta[,i]%*%sqrt(ht[,i])
    Y[,i]<-err[,i]
  }
Y<-Y[,c(-1:-200)]; #THIS TO REMOVE FIRST 200 OBSERVATIONS TO AVOID SET.SEED EFFECT
Y1<-as.matrix(Y)

coefARCH<-cbind(paraGA0, paraGAa1, paraGAa2, paraGAa3,paraGAa4,
                         paraGAa5, paraGAa6, paraGAa7,paraGAa8,
                         paraGAa9, paraGAa10, paraGAa11,paraGAa12)

LAGNONZERO<-c(1,2,5,7,11,13)   #LOCATION OF NONZERO PARAMETERS
LAGNONZERO2<-c(1,1,0,0,1,0,1,0,0,0,1,0,1)  #1: non-zero 0: zero
ALLDATA<-list(Y1,err,coefARCH,rCh,LAGNONZERO,LAGNONZERO2,sum(LAGNONZERO2))
names(ALLDATA)<-c("Y1","err","coefARCH","Gr","LAGNONZERO","LAGNONZERO2","COUNTNZ")
return(ALLDATA)
}
n=1200
ARCH1<-ARCH(n)
ERR<-as.vector(ARCH1$Y1)
#plot(ERR, type="l")
par(mfrow=c(1,2))
plot(NA, xlim=c(1,n), ylim=c(-2,2), main="Original Series", ylab="", xlab="Time")
points(ERR,col="red", type="l", pch=23,lty="dotted", cex=4)
plot(NA, xlim=c(1,n), ylim=c(0,1), main="Squared Series", ylab="", xlab="Time")
points(ERR,col="red", type="l", pch=23,lty="dotted", cex=4)
##################################################################################################

#all(c(1,2,3)%in%c(1,2,3,4))
#isTRUE(all.equal(c(1,2),c(1,3)))