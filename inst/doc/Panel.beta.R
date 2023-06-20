## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(saeHB.panel.beta)
data("dataPanelbeta")

## -----------------------------------------------------------------------------
dataPanelbeta <- dataPanelbeta[1:25,] #for the example only use part of the dataset
area <- max(dataPanelbeta[,2])
period <- max(dataPanelbeta[,3])
result<-Panel.beta(ydi~xdi1+xdi2,area=area, period=period ,iter.mcmc = 10000,thin=5,burn.in = 1000,data=dataPanelbeta)

## -----------------------------------------------------------------------------
result$Est

## -----------------------------------------------------------------------------
result$coefficient

## -----------------------------------------------------------------------------
result$refvar

## -----------------------------------------------------------------------------
MSE_HB<-result$Est$SD^2
summary(MSE_HB)

## -----------------------------------------------------------------------------
RSE_HB<-sqrt(MSE_HB)/result$Est$MEAN*100
summary(RSE_HB)

## -----------------------------------------------------------------------------
y_dir<-dataPanelbeta[,1]
y_HB<-result$Est$MEAN
y<-as.data.frame(cbind(y_dir,y_HB))
summary(y)
MSE_dir<-dataPanelbeta[,4]
MSE<-as.data.frame(cbind(MSE_dir, MSE_HB))
summary(MSE)
RSE_dir<-sqrt(MSE_dir)/y_dir*100
RSE<-as.data.frame(cbind(RSE_dir, RSE_HB))
summary(RSE)

