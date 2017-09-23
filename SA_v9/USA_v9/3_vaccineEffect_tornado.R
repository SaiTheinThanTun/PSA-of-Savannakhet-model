
#setwd("D:\\OneDrive\\MORU\\Projects\\PSA of Savannakhet model\\SA_v9\\USA_v9")
setwd("/Users/sai/OneDrive/MORU/Projects/PSA of Savannakhet model/SA_v9/USA_v9")
sourceCpp('modGMS.cpp')
source('1_shiny2ode.R')
source('runGMSFunction.R')
library(deSolve)
library(shiny)
library(TSA)
library(Rcpp)


####intervention switches####
#scenario 1
EDATon <- TRUE
ITNon <- TRUE
IRSon <- FALSE
MDAon <- TRUE
VACon <- FALSE
v_same <- TRUE
MSATon <- TRUE
primon <- FALSE

# define the number of weeks to run the model
dt<-1/12
startyear<-2007
stopyear<-2023
maxt<-stopyear-startyear
times <- seq(0, maxt, by = dt)
tsteps<-length(times)



###common input to default and min/max####
scenario_iR<-(c(EDATon = EDATon,
                ITNon = ITNon,
                IRSon = IRSon,
                MDAon = MDAon,
                primon = primon,
                MSATon = MSATon,
                VACon = as.numeric(VACon),
                v_same = v_same))

###default####
#change 2
API <- 10
bh_max <- 20
eta <- 30
covEDAT0 <- 25
covITN0 <- 70
effITN <- 30
covIRS0 <- 0
effIRS <- 15
muC <- 1
muA <- 1
muU <- 1
percfail2018 <- 5
percfail2019 <- 15
percfail2020 <- 30
EDATscale <- 1
covEDATi <- 70
ITNscale <- 1
covITNi <- 90
IRSscale <- 1
covIRSi <- 90
lossd <- 30
dm <- 6
ka <- 2
delta <- 5
cmda_1 <- 50
cmda_2 <- 50
cmda_3 <- 50
kf <- 45
ks <- 634
tm_1 <- 9
tm_2 <- 10
tm_3 <- 11
effv_3 <- 92
effv_4 <- 92
MSATscale <- 1
covMSATi <- 90
MSATsensC <- 99
MSATsensA <- 87
MSATsensU <- 4

#####default run####
default.result <- matrix(NA, 1, 2) 

#change 3
parametersR <- (c(
  bh_max = bh_max,                 # bites per human per night
  eta = eta,
  covEDAT0 = covEDAT0,
  covITN0 = covITN0,
  effITN = effITN,
  covIRS0 = covIRS0,
  effIRS = effIRS,
  muC = muC,
  muA = muA,
  muU = muU,
  percfail2018 = percfail2018,
  percfail2019 = percfail2019,
  percfail2020 = percfail2020,
  
  EDATscale = EDATscale,
  covEDATi = covEDATi,
  ITNscale = ITNscale,
  covITNi = covITNi,
  IRSscale = IRSscale,
  covIRSi = covIRSi,
  cmda_1 = cmda_1,
  cmda_2 = cmda_2,
  cmda_3 = cmda_3,
  tm_1 = tm_1,          # timing of 1st round [2018 to 2021 - 1 month steps]
  tm_2 = tm_2,          # timing of 2nd round [2018+(1/12) to 2021 - 1 month steps]
  tm_3 = tm_3,          # timing of 3rd round [2018+(2/12) to 2021 - 1 month steps]
  dm = dm,
  lossd = lossd,
  
  MSATscale = MSATscale,
  covMSATi = covMSATi,
  MSATsensC = MSATsensC,
  MSATsensA = MSATsensA,
  MSATsensU = MSATsensU,
  
  #effv_1 = effv_1,
  #effv_2 = effv_2,
  effv_3 = effv_3,
  effv_4 = effv_4,
  ka = ka,
  delta = delta,
  kf = kf,
  ks = ks
  #vh = vh
))

# initial prevalence
initprevR <- (0.001*API)

GMSout <- runGMS(initprevR, scenario_iR,parametersR)

#GMSout[nrow(GMSout),c(3,4)] #output total incidence and prevelance

default.result[1,1] <- (GMSout[GMSout[,1]==2018,2]-GMSout[nrow(GMSout),2])/GMSout[GMSout[,1]==2018,2] #GMSout[nrow(GMSout),2]
default.result[1,2] <- (GMSout[GMSout[,1]==2018,4]-GMSout[nrow(GMSout),4])/GMSout[GMSout[,1]==2018,4] #GMSout[nrow(GMSout),4]




###Scenario 2####
#source('1_shiny2ode.R')
EDATon <- TRUE
ITNon <- TRUE
IRSon <- FALSE
MDAon <- TRUE
VACon <- TRUE
v_same <- TRUE
MSATon <- TRUE
primon <- FALSE

# define the number of weeks to run the model
dt<-1/12
startyear<-2007
stopyear<-2023
maxt<-stopyear-startyear
times <- seq(0, maxt, by = dt)
tsteps<-length(times)


###common input to default and min/max####
scenario_iR<-(c(EDATon = EDATon,
                ITNon = ITNon,
                IRSon = IRSon,
                MDAon = MDAon,
                primon = primon,
                MSATon = MSATon,
                VACon = as.numeric(VACon),
                v_same = v_same))


####for####
result <- matrix(NA, nrow(valueRange)*no.s, 2) #valueTable and valueRange are from 'shiny2ode.R'

for(i in 1:nrow(simValueTable)){
  
  #change 4
  API <- simValueTable[i,1]
  bh_max <- simValueTable[i,2]
  eta <- simValueTable[i,3]
  covEDAT0 <- simValueTable[i,4]
  covITN0 <- simValueTable[i,5]
  effITN <- simValueTable[i,6]
  covIRS0 <- simValueTable[i,7]
  effIRS <- simValueTable[i,8]
  muC <- simValueTable[i,9]
  muA <- simValueTable[i,10]
  muU <- simValueTable[i,11]
  percfail2018 <- simValueTable[i,12]
  percfail2019 <- simValueTable[i,13]
  percfail2020 <- simValueTable[i,14]
  EDATscale <- simValueTable[i,15]
  covEDATi <- simValueTable[i,16]
  ITNscale <- simValueTable[i,17]
  covITNi <- simValueTable[i,18]
  IRSscale <- simValueTable[i,19]
  covIRSi <- simValueTable[i,20]
  lossd <- simValueTable[i,21]
  dm <- simValueTable[i,22]
  ka <- simValueTable[i,23]
  delta <- simValueTable[i,24]
  cmda_1 <- simValueTable[i,25]
  cmda_2 <- simValueTable[i,26]
  cmda_3 <- simValueTable[i,27]
  kf <- simValueTable[i,28]
  ks <- simValueTable[i,29]
  tm_1 <- simValueTable[i,30]
  tm_2 <- simValueTable[i,31]
  tm_3 <- simValueTable[i,32]
  effv_3 <- simValueTable[i,33]
  effv_4 <- simValueTable[i,34]
  MSATscale <- simValueTable[i,35]
  covMSATi <- simValueTable[i,36]
  MSATsensC <- simValueTable[i,37]
  MSATsensA <- simValueTable[i,38]
  MSATsensU <- simValueTable[i,39]
  
  #change 3
  parametersR <- (c(
    bh_max = bh_max,                 # bites per human per night
    eta = eta,
    covEDAT0 = covEDAT0,
    covITN0 = covITN0,
    effITN = effITN,
    covIRS0 = covIRS0,
    effIRS = effIRS,
    muC = muC,
    muA = muA,
    muU = muU,
    percfail2018 = percfail2018,
    percfail2019 = percfail2019,
    percfail2020 = percfail2020,
    
    EDATscale = EDATscale,
    covEDATi = covEDATi,
    ITNscale = ITNscale,
    covITNi = covITNi,
    IRSscale = IRSscale,
    covIRSi = covIRSi,
    cmda_1 = cmda_1,
    cmda_2 = cmda_2,
    cmda_3 = cmda_3,
    tm_1 = tm_1,          # timing of 1st round [2018 to 2021 - 1 month steps]
    tm_2 = tm_2,          # timing of 2nd round [2018+(1/12) to 2021 - 1 month steps]
    tm_3 = tm_3,          # timing of 3rd round [2018+(2/12) to 2021 - 1 month steps]
    dm = dm,
    lossd = lossd,
    
    MSATscale = MSATscale,
    covMSATi = covMSATi,
    MSATsensC = MSATsensC,
    MSATsensA = MSATsensA,
    MSATsensU = MSATsensU,
    
    #effv_1 = effv_1,
    #effv_2 = effv_2,
    effv_3 = effv_3,
    effv_4 = effv_4,
    ka = ka,
    delta = delta,
    kf = kf,
    ks = ks
    #vh = vh
  ))
  
  # initial prevalence
  initprevR <- (0.001*API)
  
  GMSout <- runGMS(initprevR, scenario_iR,parametersR)
  
  #GMSout[nrow(GMSout),c(3,4)] #output total incidence and prevelance
  
  result[i,1] <- (GMSout[GMSout[,1]==2018,2]-GMSout[nrow(GMSout),2])/GMSout[GMSout[,1]==2018,2] #GMSout[nrow(GMSout),2]
  result[i,2] <- (GMSout[GMSout[,1]==2018,4]-GMSout[nrow(GMSout),4])/GMSout[GMSout[,1]==2018,4] #GMSout[nrow(GMSout),4]
}

####subsetting vaccine parameters####
par_vaccine <- c(22,23,24,28,29,33,34) #dm, ka, delta, kf, ks, effv_3, effv_4
par_vaccine_double <- sort(c(par_vaccine*2,(par_vaccine*2+1)))
result <- result[par_vaccine_double,]
simValueTable <- simValueTable[par_vaccine_double,]

####arranging results with vaccine on####
incidenceT <- matrix(result[,1], nrow(result)/no.s,no.s, byrow = T)
prevalenceT <- matrix(result[,2], nrow(result)/no.s,no.s, byrow = T)


###write scenario comparison results####

incidenceDiff <- incidenceT-default.result[,1]
incidenceLowValue <- incidenceDiff[,1]
incidenceHiValue <- incidenceDiff[,2]

incidenceRange <- cbind(as.data.frame(cbind(incidenceLowValue,incidenceHiValue)), valueTable[par_vaccine,1])
rangeSize <- abs(incidenceRange[,2]-incidenceRange[,1])
incidenceRange <- incidenceRange[order(rangeSize),]


prevalenceDiff <- prevalenceT-default.result[,2]
prevalenceLowValue <- prevalenceDiff[,1]
prevalenceHiValue <- prevalenceDiff[,2]

prevalenceRange <- cbind(as.data.frame(cbind(prevalenceLowValue,prevalenceHiValue)), valueTable[par_vaccine,1])
rangeSize2 <- abs(prevalenceRange[,2]-prevalenceRange[,1])
prevalenceRange <- prevalenceRange[order(rangeSize2),]

write.csv(incidenceRange,paste('result/incidenceRange',gsub(':','_',Sys.time()),'.csv', sep = ''))
write.csv(prevalenceRange,paste('result/prevalenceRange',gsub(':','_',Sys.time()),'.csv', sep = '')) 

####plots####
png(file=paste('result/compare_incidenceTornado_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(mar=c(5, 7, 4, 2), cex=1.25)
barplot(incidenceRange[,1], names.arg = incidenceRange[,3], las=1, horiz=T, xlim=c(-.05,0.05), beside = T, axes=F, col='light blue', main=paste("Sensitivity on incidence, Savannakhet model \n EDAT+ITN+MDA vs EDAT+ITN+MDA+VAC"), xlab='% reduction in incidence')
barplot(incidenceRange[,2], horiz=T, axes=F, beside=T, add=T, col=adjustcolor( "gold", alpha.f = 0.7)) #col='gold', alpha=.8)
axis(1, at=seq(-.05,0.05,by=.01), labels = 100*seq(-.05,0.05,by=.01)+round(default.result[,1]*100,0), ylab='% reduction in incidence')
legend(x=-.045,y=3,legend=c('lower value of parameter','higher value of parameter'), fill=c('light blue','gold'), cex=1)
dev.off()




png(file=paste('result/compare_prevalenceTornado_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(mar=c(5, 7, 4, 2), cex=1.25)
barplot(prevalenceRange[,1], names.arg = prevalenceRange[,3], axes=F,las=1, horiz=T, xlim=c(-.05,0.05), beside = T, col='light blue', main=paste("Sensitivity on prevalence, Savannakhet model \n EDAT+ITN+MDA vs EDAT+ITN+MDA+VAC"), xlab='% reduction in prevalence')
barplot(prevalenceRange[,2], horiz=T, axes=F, beside=T, add=T, col=adjustcolor( "gold", alpha.f = 0.7))
axis(1, at=seq(-.05,0.05,by=.01), labels = 100*seq(-.05,0.05,by=.01)+round(default.result[,2]*100,0), ylab='% reduction in prevalence')
legend(x=-.045,y=3,legend=c('lower value of parameter','higher value of parameter'), fill=c('light blue','gold'), cex=1)
dev.off()