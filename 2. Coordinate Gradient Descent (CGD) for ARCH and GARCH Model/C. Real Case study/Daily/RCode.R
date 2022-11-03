library(cowplot)
library(readxl)
library(ggplot2)
library(forecast)
ASX_Ordinary <- read_excel("C:/Users/Jaffri/Dropbox/PhD Research/4. Third Year/3. Thesis Write up/Simulation Studies-Rcode/B. Pure GARCH models/C. Real Case study/Daily/ASX Ordinary.xlsx", col_types = c("date", "numeric", "numeric", "numeric", "numeric"))
library("ggplot2", lib.loc="~/R/win-library/3.5")
plot1<-ggplot(data = ASX_Ordinary, aes(x = Date, y = Close))+geom_line(color = "black", size = 1) + labs(title = "ASX Ordinary", y = "Daily Closing Price", x="Time")
Return<-as.matrix(100*diff(log(ASX_Ordinary$Close), differences =1))
Date1<-ASX_Ordinary$Date
ASX_Return<-data.frame(Date1[2:length(Date1)],Return); 
colnames(ASX_Return)[1]<-"Date"
plot2<-ggplot(data = ASX_Return, aes(x = Date, y = Return))+geom_line(color = "red", size = 1) + labs(title = "ASX Ordinary Log Return", y="Log Return, Daily Closing Price (%)", x="Time")
plot3<-ggplot(data = ASX_Return, aes(x = Date, y = Return^2))+geom_line(color = "blue", size = 1) + labs(title = "Squared Log Return", y="", x="Time")
ABOVEP<-plot_grid(plot1, plot2, labels = c('A', 'B'))
plot_grid(ABOVEP, plot3, labels = c('', 'C') , ncol = 1)

ggAcf(Return^2,lag.max = 30, main="ACF of Squared ASX Ordinary Log Return")
ggPacf(Return^2,lag.max = 30, main="PACF of Squared ASX Ordinary Log Return")

#par(mfrow=c(2,1))
#plot((Return^2), type="l", main="Squared Return", ylab="")
#pacf(Return^2,lag=60, main="PACF of Squared ASX Ordinary Log Return")

ERR2<-Return^2
V<-10
MWin<-numeric()
for(II in (V+1):length(ERR2)){
MWin[II-V]<-sum(ERR2[(II-V):(II)])/V
}

MOv_Win<-data.frame(Date1[(V+2):length(Date1)],MWin); 
colnames(MOv_Win)[1]<-"Date"; colnames(MOv_Win)[2]<-"MW"
SigARCH <- read.csv("C:/Users/Jaffri/Dropbox/PhD Research/4. Third Year/3. Thesis Write up/Simulation Studies-Rcode/B. Pure GARCH models/C. Real Case study/Daily/SigARCH.csv")
SigARCH<-SigARCH[,2]; 
Fit_ARCH<-data.frame(Date1[(V+2):length(Date1)],SigARCH[(V+1):length(SigARCH)]);
colnames(Fit_ARCH)[1]<-"Date"; colnames(Fit_ARCH)[2]<-"Fit"
SigGARCH <- read.csv("C:/Users/Jaffri/Dropbox/PhD Research/4. Third Year/3. Thesis Write up/Simulation Studies-Rcode/B. Pure GARCH models/C. Real Case study//Daily/SigGARCH.csv")
SigGARCH<-SigGARCH[,2]; 
Fit_GARCH<-data.frame(Date1[(V+2):length(Date1)],SigGARCH[(V+1):length(SigGARCH)]);
colnames(Fit_GARCH)[1]<-"Date"; colnames(Fit_GARCH)[2]<-"Fit"

plotA1<-ggplot(data = MOv_Win, aes(x = Date, y = MW))+geom_line(color = "black", size = 1) + labs(title = "Moving Windows of Squared Return", y="", x="Time")
plotA2<-ggplot(data = Fit_ARCH, aes(x = Date, y = Fit))+geom_line(color = "black", size = 1) + labs(title = "Fitted ARCH Model", y="", x="Time")
plotA3<-ggplot(data = Fit_GARCH, aes(x = Date, y = Fit))+geom_line(color = "black", size = 1) + labs(title = "Fitted GARCH Model", y="", x="Time")
plot_grid(plotA1, plotA2, plotA3, labels = c('A', 'B' ,'C') , ncol = 1)