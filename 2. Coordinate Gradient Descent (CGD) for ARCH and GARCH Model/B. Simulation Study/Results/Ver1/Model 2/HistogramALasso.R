par(mfrow=c(2,3))


hist(colPARAmeanadapLASSO[1,], breaks=20, 
	xlab=expression(paste(alpha,"[", 0, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[1],col="red", lwd=3, lty=2)

hist(colPARAmeanadapLASSO[2,], breaks=20, 
	xlab=expression(paste(alpha,"[", 1, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[2],col="red", lwd=3, lty=2)

hist(colPARAmeanadapLASSO[3,], breaks=20, 
	xlab=expression(paste(beta,"[", 1, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[3],col="red", lwd=3, lty=2)

hist(colPARAmeanadapLASSO[4,], breaks=20, 
	xlab=expression(paste(beta,"[", 2, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[4],col="red", lwd=3, lty=2)

hist(colPARAmeanadapLASSO[5,], breaks=20, 
	xlab=expression(paste(beta,"[", 3, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[5],col="red", lwd=3, lty=2)