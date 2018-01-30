
setwd("C:/Users/freek_000/Documents/PhD_Otto/AileneProject/PopGenDensityDependence/")
library(Rcpp)
library(plot3D)
library(animation)

sourceCpp("RcppNumericIntegration.cpp")

input <-  matrix(0, nrow = 200, ncol = 200)
input[10,10] = 1.0

z <- input
for(i in 0:20){
  name <- paste(i,".png",sep="")
  png(name)
  z = ItteratorNP(z, 0.5 ,0.4 ,1)
  image2D(z=z)
  dev.off()
}


