source('optimal_withinHostTrans_triple_Median_Symptom.R')
source('withinHostTrans.R')
source('withinHostTrans_RCpp240.R')
library("deSolve")
library("adaptMCMC")
library(parallel)

lb = c( -1,  0,    1,  -1,  -1,   -1,  -1)-0.5
ub = c(  1,  5,     6,  1,   1,    1,   1)+0.5
x0 = c(0.35,  3.5,  3.5,  0.5, -0.9, 0.5,  0.35)+runif(1)*0.1
# set.seed(runif(1))

dimension <- 7
global.min <- 0
tol <- 1e-13
lower <- lb
upper <- ub
out_SA = GenSA(par = x0, lower = lb, upper = ub, fn = optimal_withinHostTrans_triple_NoMedian_minimize,
               control=list(threshold.stop=global.min+tol, nb.stop.improvement=10^4, max.call=10^4,# max.call=10^4,
                            verbose=TRUE, max.time=3600), 
               PopDist_Baloxavir=data$PopDist.Baloxavir, BaloxavirScaleHour=data$BaloxavirScaleHour,
               HourScaling=1, PopDist_Placebo=data$PopDist.Placebo, PlaceboScaleHour=data$PlaceboScaleHour,
               PopDist_Oseltamivir=data$PopDist.Oseltamivir, OseltamivirScaleHour=data$OseltamivirScaleHour,
               lb=lb, ub=ub, returnLog=100, eiSingle=1, parasNum=length(x0), otherParas=c()) 
