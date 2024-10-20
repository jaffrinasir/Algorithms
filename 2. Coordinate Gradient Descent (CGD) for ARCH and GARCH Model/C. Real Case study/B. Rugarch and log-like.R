#### LOAD R LIBRARY  ########################################################
library(mvtnorm)
library(MASS)
library(fBasics)
library(matrixcalc)
library(rugarch)
library(lmtest)
#############################################################################
library(cowplot)
library(readxl)
library(ggplot2)
library(forecast)
ASX_Ordinary <- read_excel("P:/1. Draft Research Articles/2. Paper (GARCH)/(R code) Pure GARCH models/C. Real Case study/Daily/ASX Ordinary.xlsx", col_types = c("date", "numeric", "numeric", "numeric", "numeric"))
library("ggplot2", lib.loc="~/R/win-library/3.5")
plot1<-ggplot(data = ASX_Ordinary, aes(x = Date, y = Close))+geom_line(color = "black", size = 1) + labs(title = "ASX Ordinary", y = "Daily Closing Price", x="Time")
Return<-as.matrix(100*diff(log(ASX_Ordinary$Close), differences =1))
TRAIN<-90 #IN PERCENTAGE
ERR<-as.vector(Return)[1:(TRAIN/100*length(Return))]

garchspec1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), submodel = NULL, external.regressors = NULL, 
	variance.targeting = FALSE), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE, archm = FALSE))
garchspec2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 2), submodel = NULL, external.regressors = NULL, 
	variance.targeting = FALSE), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE, archm = FALSE))
garchspec3<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 1), submodel = NULL, external.regressors = NULL, 
	variance.targeting = FALSE), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE, archm = FALSE))
garchspec4<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2), submodel = NULL, external.regressors = NULL, 
	variance.targeting = FALSE), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE, archm = FALSE))

GARCH1n1<-ugarchfit(garchspec1, ERR)
GARCH1n2<-ugarchfit(garchspec2, ERR)
GARCH2n1<-ugarchfit(garchspec3, ERR)
GARCH2n2<-ugarchfit(garchspec4, ERR)

pvalA<- pchisq(-2*(GARCH1n1@fit$LLH-GARCH1n2@fit$LLH),df = 1, lower.tail = F)
pvalB<- pchisq(-2*(GARCH1n1@fit$LLH-GARCH2n1@fit$LLH),df = 1, lower.tail = F)
pvalC<- pchisq(-2*(GARCH1n1@fit$LLH-GARCH2n2@fit$LLH),df = 2, lower.tail = F)

