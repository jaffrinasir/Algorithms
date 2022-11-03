##FUNCTION pure GARCH(1,3)
GARCH<-function(n){
  #A.] GENERATE SIMULATED DATA WITH PREDEFINED PARAMETERS
  library(mvtnorm)
  library(MASS)
  library(fBasics)
  library(matrixcalc)
  rCh<- 1
  sCh<- 3 
  n<-n+200  #Obs. below 201 will be deleted 
  Q<- 1 #Number of Variate
  
  ####### Generate r.v of white noise #########################################
  eta<-rmvnorm(n, mean=c(0), diag(1)) 
  eta<-t(eta)
  #############################################################################
  
  ##INTERCEPT
  paraGA0<-matrix(c(0.5),nrow=Q,byrow=TRUE)
  #### MA(p) representative
  paraGAa1<-matrix(c(0.1),nrow=Q,byrow=TRUE)
  #### AR(q) representative
  paraGAb1<-matrix(c(0.1),nrow=Q,byrow=TRUE)
  paraGAb2<-matrix(c(0.0),nrow=Q,byrow=TRUE)
  paraGAb3<-matrix(c(0.55),nrow=Q,byrow=TRUE)

  H<-max(rCh,sCh)
  ht<-abs(t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1)))))  #HETEROS. OBS
  err<-t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1)))) #ERR. TERM
  Y <-t(rmvnorm(n=H, mean=c(0), sigma=matrix(c(1))))
  ZERO<-matrix(0,nrow=Q,ncol=(dim(eta)[2]-dim(Y)[2]))
  ht<-cbind(ht,ZERO) 
  err<-cbind(err,ZERO)
  Y<-cbind(Y,ZERO) 

GARCHr<-sapply( paste('paraGAa', 1:rCh, sep=''), get ,envir=sys.frame(sys.parent(0)))
GARCHs<-sapply( paste('paraGAb', 1:sCh, sep=''), get ,envir=sys.frame(sys.parent(0)))

    for (i in (H+1):((dim(eta)[2]))) {
      sumRg<-sumSg<-0
       for(j in 1:rCh){
        sumRg<-sumRg + (GARCHr[[j]]%*%((err[,(i-j)])^{2}))
       }
       for(k in 1:sCh){
        sumSg<-sumSg + (GARCHs[[k]]%*%ht[,(i-k)])
       }
     
    ht[,i]<-paraGA0+(sumRg)+(sumSg)
    err[,i]<-eta[,i]%*%sqrt(ht[,i])
    Y[,i]<-err[,i]
  }
Y<-Y[,c(-1:-200)]; #THIS TO REMOVE FIRST 200 OBSERVATIONS TO AVOID SET.SEED EFFECT
Y1<-as.matrix(Y)

coefGARCH<-cbind(paraGA0, paraGAa1, 
  paraGAb1,paraGAb2,paraGAb3)

LAGNONZERO<-c(1,2,3,5)
LAGNONZERO2<-c(1,1,1,0,1)  #1: non-zero 0: zero
ALLDATA<-list(Y1,err,coefGARCH,rCh,sCh,LAGNONZERO,LAGNONZERO2,sum(LAGNONZERO2),"LAGNONZERO","LAGNONZERO2","COUNTNZ")
names(ALLDATA)<-c("Y1","err","coefGARCH","Gr","Gs")
return(ALLDATA)
}
n=1500
GARCH1<-GARCH(n)
ERR<-as.vector(GARCH1$Y1)

plot(NA, xlim=c(1,n), ylim=c(-2.1,2.1), main="Original Series", ylab="", xlab="Time")
points(ERR,col="red", type="l", pch=23,lty="dotted", cex=4)
##################################################################################################
#all(c(1,2,3)%in%c(1,2,3,4))==TRUE
#isTRUE(all.equal(c(1,2),c(1,3)))