setwd("C:/Users/Jaffri/Dropbox/Journal Papers drafts/2020/Simulation Studies-Rcode/B. Pure GARCH models/B. Simulation Study/Special")
#setwd("C:/Users/mjaf8/Dropbox/Journal Papers drafts/2020/Simulation Studies-Rcode/B. Pure GARCH models/B. Simulation Study/Special")

### LOAD R PACKAGES #######################################
library(rpgm)
library(robustbase)
library(mvtnorm)
library(MASS)
library(fBasics)
library(matrixcalc)
###########################################################

### LOAD R FUNCTIONS #########################################################
Mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux)); 
  if(length(x)==length(ux[tab == max(tab)])){print("No mode")}else{ux[tab == max(tab)]}
}
##############################################################################

#### COUNT ESTIMATED NONZERO VARIABLES
cNonZerofull<-cNonZeroadaplasso<-integer()
cBIASNonZerofull<-cBIASNonZeroadaplasso<-integer()
cEVARNonZerofull<-cEVARNonZeroadaplasso<-integer()
#### CHECK TRUE OR FALSE MODEL
pickAllSigfull<-pickAllSigadaplasso<-0
isTRUEinfull<-isTRUEinadaplasso<-0
#### VARIABLES FOR ESTIMATED PARAMETERS
colESTfull<-colESTadaplasso<-list()
colBIASfull<-colBIASadaplasso<-list()
colEVARfull<-colEVARadaplasso<-list()
cpropoESTfull<-cpropoESTadaplasso<-list()

for(ii in 1:100){
  load(paste("./LIST(1)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 101:200){
  load(paste("./LIST(2)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 201:300){
  load(paste("./LIST(3)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 301:400){
  load(paste("./LIST(4)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 401:500){
  load(paste("./LIST(5)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 501:600){
  load(paste("./LIST(6)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 601:700){
  load(paste("./LIST(7)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 701:800){
  load(paste("./LIST(8)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 801:900){
  load(paste("./LIST(9)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(ii in 901:1000){
  load(paste("./LIST(10)/SIMUGARCH.(", ii ,").RData", sep=""))
  assign(paste("LISTGARCH",ii,sep=""),LISTA)
  rm(LISTA)
}

for(REPLICATE in 1:1000){
LISTB<-lapply( paste('LISTGARCH', REPLICATE, sep=''), get ,envir=sys.frame(sys.parent(0)))
### BEGIN : FULL CONSTRAINT QMLE ############################################
colESTfull[[REPLICATE]]<-LISTB[[1]][[1]][[1]][[1]]
colBIASfull[[REPLICATE]]<-LISTB[[1]][[1]][[1]][[2]]
colEVARfull[[REPLICATE]]<-LISTB[[1]][[1]][[1]][[3]]
pickAllSigfull[REPLICATE]<-LISTB[[1]][[1]][[1]][[4]]
isTRUEinfull[REPLICATE]<-LISTB[[1]][[1]][[1]][[5]]
cNonZerofull[REPLICATE]<-LISTB[[1]][[1]][[1]][[6]]
cBIASNonZerofull[REPLICATE]<- LISTB[[1]][[1]][[1]][[7]]
cEVARNonZerofull[REPLICATE]<- LISTB[[1]][[1]][[1]][[8]]
cpropoESTfull[[REPLICATE]]<-LISTB[[1]][[1]][[1]][[9]]
### END : FULL CONSTRAINT QMLE ##############################################
#### COLLECT RESULT: ADAPTIVE LASSO
colESTadaplasso[[REPLICATE]]<-LISTB[[1]][[1]][[2]][[1]]
colBIASadaplasso[[REPLICATE]]<-LISTB[[1]][[1]][[2]][[2]]
colEVARadaplasso[[REPLICATE]]<-LISTB[[1]][[1]][[2]][[3]]
pickAllSigadaplasso[REPLICATE]<-LISTB[[1]][[1]][[2]][[4]]
isTRUEinadaplasso[REPLICATE]<-LISTB[[1]][[1]][[2]][[5]]
cNonZeroadaplasso[REPLICATE]<-LISTB[[1]][[1]][[2]][[6]]
cBIASNonZeroadaplasso[REPLICATE]<- LISTB[[1]][[1]][[2]][[7]]
cEVARNonZeroadaplasso[REPLICATE]<- LISTB[[1]][[1]][[2]][[8]]
cpropoESTadaplasso[[REPLICATE]]<-LISTB[[1]][[1]][[2]][[9]]
### END : ADAPTIVE LASSO  ##############################################
rm(LISTB)
rm(list=paste("LISTGARCH",REPLICATE,sep=""))
}

#### COLLECT RESULTS ########################################################################
for(REPLICATE in 1:1000){
  if(REPLICATE==1){
    ### RESULT: FULL QMLE
    colPARAmeanQMLE<-cbind(colESTfull[[REPLICATE]])
    colPARAbiasQMLE<-cbind(t(colBIASfull[[REPLICATE]]))
    colPARAVARQMLE<-cbind(t(colEVARfull[[REPLICATE]]))
    colPropofull<-cbind(cpropoESTfull[[REPLICATE]])
    ### RESULT: ADAPTIVE LASSO QMLE
    colPARAmeanadapLASSO<-cbind(colESTadaplasso[[REPLICATE]])
    colPARAbiasadapLASSO<-cbind(t(colBIASadaplasso[[REPLICATE]]))
    colPARAVARadapLASSO<-cbind(t(colEVARadaplasso[[REPLICATE]]))
    colPropoadapLASSO<-cbind(cpropoESTadaplasso[[REPLICATE]])
  }else{
    ### RESULT: FULL QMLE
    colPARAmeanQMLE<-cbind(colPARAmeanQMLE,as.matrix(colESTfull[[REPLICATE]]))
    colPARAbiasQMLE<-cbind(colPARAbiasQMLE,t(colBIASfull[[REPLICATE]]))
    colPARAVARQMLE<-cbind(colPARAVARQMLE,t(colEVARfull[[REPLICATE]]))
    colPropofull<-cbind(colPropofull,as.matrix(cpropoESTfull[[REPLICATE]]))
    ### RESULT: ADAPTIVE LASSO QMLE
    colPARAmeanadapLASSO<-cbind(colPARAmeanadapLASSO,as.matrix(colESTadaplasso[[REPLICATE]]))
    colPARAbiasadapLASSO<-cbind(colPARAbiasadapLASSO,t(colBIASadaplasso[[REPLICATE]]))
    colPARAVARadapLASSO<-cbind(colPARAVARadapLASSO,t(colEVARadaplasso[[REPLICATE]]))
    colPropoadapLASSO<-cbind(colPropoadapLASSO,as.matrix(cpropoESTadaplasso[[REPLICATE]]))
  }
}

### RESULT: FULL QMLE ##################################################
## COUNT NON-ZERO FOR FULL QMLE
min(cNonZerofull)                #MINIMUM
max(cNonZerofull)                #MAXIMUM
round(mean(cNonZerofull))        #MEAN (rounded)
Mode(cNonZerofull)               #MODE
sum(cBIASNonZerofull)/1000  #BIAS
sqrt(sum(cEVARNonZerofull)/1000)  #ESD

table(pickAllSigfull)
table(isTRUEinfull)

PARAmeanQMLE<-rowSums(colPARAmeanQMLE)/1000; 
rowMins(colPARAmeanQMLE)
rowMedians(colPARAmeanQMLE)
rowMaxs(colPARAmeanQMLE)
PARAbiasQMLE<-rowSums(colPARAbiasQMLE)/1000; 
PARAesdQMLE<-sqrt(rowSums(colPARAVARQMLE)/1000); 
PARApropoQMLE<-rowSums(colPropofull)/1000; 
########################################################################

### RESULT: ADAPTIVE LASSO #############################################
## COUNT NON-ZERO FOR ADAPTIVE LASSO
min(cNonZeroadaplasso)                #MINIMUM
max(cNonZeroadaplasso)                #MAXIMUM
round(mean(cNonZeroadaplasso))        #MEAN (rounded)
Mode(cNonZeroadaplasso)               #MODE
sum(cBIASNonZeroadaplasso)/1000  #BIAS
sqrt(sum(cEVARNonZeroadaplasso)/1000)  #ESD

table(pickAllSigadaplasso)
table(isTRUEinadaplasso)

PARAmeanadapLASSO<-rowSums(colPARAmeanadapLASSO)/1000; 
rowMins(colPARAmeanadapLASSO)
rowMedians(colPARAmeanadapLASSO)
rowMaxs(colPARAmeanadapLASSO)
PARAbiasadapLASSO<-rowSums(colPARAbiasadapLASSO)/1000; 
PARAesdadapLASSO<-sqrt(rowSums(colPARAVARadapLASSO)/1000); 
PARApropadapLASSO<-rowSums(colPropoadapLASSO)/1000;
##########################################################################

TEXT01<-"#### RESULT: FULL QMLE ####" 
TEXT02<-"## COUNT NON-ZERO FOR FULL QMLE ##"
TEXT03<-"## PARAMETER ESTIMATES: FULL QMLE ##"
TEXT04<-"#### RESULT: ADAPTIVE LASSO ####"  
TEXT05<-"## COUNT NON-ZERO FOR ADAPTIVE LASSO ##"
TEXT06<-"## PARAMETER ESTIMATES: ADAPTIVE LASSO ##"

LISTFULL0<-list(TEXT01, TEXT02, min(cNonZerofull), max(cNonZerofull), round(mean(cNonZerofull)), Mode(cNonZerofull),
                sum(cBIASNonZerofull)/1000, sqrt(sum(cEVARNonZerofull)/1000),
                table(pickAllSigfull), table(isTRUEinfull), 
                TEXT03, PARAmeanQMLE, rowMins(colPARAmeanQMLE), rowMedians(colPARAmeanQMLE),
                rowMaxs(colPARAmeanQMLE), PARAbiasQMLE, PARAesdQMLE, PARApropoQMLE,
                TEXT04, TEXT05, min(cNonZeroadaplasso), max(cNonZeroadaplasso), round(mean(cNonZeroadaplasso)), Mode(cNonZeroadaplasso),
                sum(cBIASNonZeroadaplasso)/1000, sqrt(sum(cEVARNonZeroadaplasso)/1000),
                table(pickAllSigadaplasso), table(isTRUEinadaplasso),
                TEXT06, PARAmeanadapLASSO, rowMins(colPARAmeanadapLASSO), rowMedians(colPARAmeanadapLASSO), 
                rowMaxs(colPARAmeanadapLASSO), PARAbiasadapLASSO,  PARAesdadapLASSO, PARApropadapLASSO
)

names(LISTFULL0)<-c( "", "", "Minimum", "Maximum", "Average", "Mode",
                    "Bias", "ESD",
                    "Pick All Relevant", "Is Correct model??", 
                    "", "Est. Avg. Full", "Est. Min. Full", "Est. Med. Full",
                    "Est. Max. Full", "Est. Bias Full", "Est. ESD Full", "Est. propo Full",
                    "", "", "Minimum", "Maximum", "Average", "Mode",
                    "Bias", "ESD",
                    "Pick All Relevant", "Is Correct model??",
                    "", "Est. Avg. A.LASSO", "Est. Min. A.LASSO", "Est. Med. A.LASSO", 
                    "Est. Max. A.LASSO", "Est. Bias A.LASSO",  "Est. ESD A.LASSO", "Est. propo A.LASSO"
)

setwd("C:/Users/Jaffri/Dropbox/Journal Papers drafts/2020/Simulation Studies-Rcode/B. Pure GARCH models/B. Simulation Study")
OUTA<-capture.output(list(LISTFULL0)); 
cat("n=",1200,"Replication=",1000,OUTA,file="./Results/Ver1/Model 1/Model1n1200q50s3.txt",sep="\n", append = FALSE)