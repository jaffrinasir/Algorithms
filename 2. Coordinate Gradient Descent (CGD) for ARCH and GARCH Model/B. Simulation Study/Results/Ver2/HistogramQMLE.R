par(mfrow=c(5,3))

for(ii in 1:13){
hist(colPARAmeanQMLE[ii,], breaks=20, xlab=expression(paste("alpha[", ii-1, "]", sep="")), main="",col="lightblue"); 
#xlab=expression(paste("alpha[", ii-1, "]")
abline(v=TRUEPARA[ii],col="red", lwd=3, lty=2)
}

par(mfrow=c(5,3))


hist(colPARAmeanQMLE[1,], breaks=20, 
	xlab=expression(paste(alpha,"[", 0, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[1],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[2,], breaks=20, 
	xlab=expression(paste(alpha,"[", 1, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[2],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[3,], breaks=20, 
	xlab=expression(paste(alpha,"[", 2, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[3],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[4,], breaks=20, 
	xlab=expression(paste(alpha,"[", 3, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[4],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[5,], breaks=20, 
	xlab=expression(paste(alpha,"[", 4, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[5],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[6,], breaks=20, 
	xlab=expression(paste(alpha,"[", 5, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[6],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[7,], breaks=20, 
	xlab=expression(paste(alpha,"[", 6, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[7],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[8,], breaks=20, 
	xlab=expression(paste(alpha,"[", 7, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[8],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[9,], breaks=20, 
	xlab=expression(paste(alpha,"[", 8, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[9],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[10,], breaks=20, 
	xlab=expression(paste(alpha,"[", 9, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[10],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[11,], breaks=20, 
	xlab=expression(paste(alpha,"[", 10, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[11],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[12,], breaks=20, 
	xlab=expression(paste(alpha,"[", 11, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[12],col="red", lwd=3, lty=2)

hist(colPARAmeanQMLE[13,], breaks=20, 
	xlab=expression(paste(alpha,"[", 12, "]", sep="")), main="",col="lightblue"); 
abline(v=TRUEPARA[13],col="red", lwd=3, lty=2)