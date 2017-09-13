############################################################
####Scenario analysis: Vaccine effect######
#Developed by Sai Thein Than Tun, sai@tropmedres.ac
#First development date: 2017 June 20
#Updated: 2017 September 12
#run the model for default value and then min and max values varying each value at a time

# library(deSolve)
# library(shiny)
# library(TSA)
# library(Rcpp)

#setwd("D:\\OneDrive\\MORU\\Projects\\PSA of Savannakhet model\\SA_v9\\USA_v9")
setwd("/Users/sai/OneDrive/MORU/Projects/PSA of Savannakhet model/SA_v9/USA_v9")
#sourceCpp('modGMS.cpp')
source('1_shiny2ode.R')
source('runGMSFunction.R')

####"min values"####
####intervention switches####
####scenario 1: min values####
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
API <- 1
bh_max <- 15
eta <- 5
covEDAT0 <- 10
covITN0 <- 50
effITN <- 20
covIRS0 <- 0
effIRS <- 5
muC <- 0.5
muA <- 0.5
muU <- 0.5
percfail2018 <- 2
percfail2019 <- 7
percfail2020 <- 15
EDATscale <- 0.25
covEDATi <- 60
ITNscale <- 0.25
covITNi <- 70
IRSscale <- 0.25
covIRSi <- 70
lossd <- 15
dm <- 1
ka <- 1
delta <- 1
cmda_1 <- 30
cmda_2 <- 30
cmda_3 <- 30
kf <- 15
ks <- 500
tm_1 <- 9
tm_2 <- 10
tm_3 <- 11
effv_3 <- 70
effv_4 <- 70
MSATscale <- 0.25
covMSATi <- 70
MSATsensC <- 90
MSATsensA <- 70
MSATsensU <- 1

#####default run####
min_default.result0 <- matrix(NA, 1, 2) 

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

min_default.result0[1,1] <- (GMSout[GMSout[,1]==2018,2]-GMSout[nrow(GMSout),2])/GMSout[GMSout[,1]==2018,2] #GMSout[nrow(GMSout),2]
min_default.result0[1,2] <- (GMSout[GMSout[,1]==2018,4]-GMSout[nrow(GMSout),4])/GMSout[GMSout[,1]==2018,4] #GMSout[nrow(GMSout),4]

###Scenario 2: min values####
EDATon <- TRUE
ITNon <- TRUE
IRSon <- FALSE
MDAon <- TRUE
VACon <- TRUE
v_same <- TRUE
MSATon <- TRUE
primon <- FALSE

###common input to default and min/max####
scenario_iR<-(c(EDATon = EDATon,
                ITNon = ITNon,
                IRSon = IRSon,
                MDAon = MDAon,
                primon = primon,
                MSATon = MSATon,
                VACon = as.numeric(VACon),
                v_same = v_same))

#####default run####
min_default.result1 <- matrix(NA, 1, 2) 

GMSout <- runGMS(initprevR, scenario_iR,parametersR)

#GMSout[nrow(GMSout),c(3,4)] #output total incidence and prevelance

min_default.result1[1,1] <- (GMSout[GMSout[,1]==2018,2]-GMSout[nrow(GMSout),2])/GMSout[GMSout[,1]==2018,2] #GMSout[nrow(GMSout),2]
min_default.result1[1,2] <- (GMSout[GMSout[,1]==2018,4]-GMSout[nrow(GMSout),4])/GMSout[GMSout[,1]==2018,4] #GMSout[nrow(GMSout),4]




####"max values"####
####intervention switches####
####scenario 1: max values####
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
API <- 30
bh_max <- 25
eta <- 60
covEDAT0 <- 50
covITN0 <- 90
effITN <- 50
covIRS0 <- 90
effIRS <- 25
muC <- 10
muA <- 20
muU <- 20
percfail2018 <- 10
percfail2019 <- 20
percfail2020 <- 45
EDATscale <- 3
covEDATi <- 95
ITNscale <- 3
covITNi <- 90
IRSscale <- 3
covIRSi <- 90
lossd <- 30
dm <- 12
ka <- 8
delta <- 12
cmda_1 <- 70
cmda_2 <- 70
cmda_3 <- 70
kf <- 60
ks <- 800
tm_1 <- 9
tm_2 <- 10
tm_3 <- 11
effv_3 <- 99
effv_4 <- 99
MSATscale <- 3
covMSATi <- 95
MSATsensC <- 100
MSATsensA <- 95
MSATsensU <- 30

#####default run####
max_default.result0 <- matrix(NA, 1, 2) 

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

max_default.result0[1,1] <- (GMSout[GMSout[,1]==2018,2]-GMSout[nrow(GMSout),2])/GMSout[GMSout[,1]==2018,2] #GMSout[nrow(GMSout),2]
max_default.result0[1,2] <- (GMSout[GMSout[,1]==2018,4]-GMSout[nrow(GMSout),4])/GMSout[GMSout[,1]==2018,4] #GMSout[nrow(GMSout),4]

###Scenario 2: min values####
EDATon <- TRUE
ITNon <- TRUE
IRSon <- FALSE
MDAon <- TRUE
VACon <- TRUE
v_same <- TRUE
MSATon <- TRUE
primon <- FALSE

###common input to default and min/max####
scenario_iR<-(c(EDATon = EDATon,
                ITNon = ITNon,
                IRSon = IRSon,
                MDAon = MDAon,
                primon = primon,
                MSATon = MSATon,
                VACon = as.numeric(VACon),
                v_same = v_same))

#####default run####
max_default.result1 <- matrix(NA, 1, 2) 

GMSout <- runGMS(initprevR, scenario_iR,parametersR)

#GMSout[nrow(GMSout),c(3,4)] #output total incidence and prevelance

max_default.result1[1,1] <- (GMSout[GMSout[,1]==2018,2]-GMSout[nrow(GMSout),2])/GMSout[GMSout[,1]==2018,2] #GMSout[nrow(GMSout),2]
max_default.result1[1,2] <- (GMSout[GMSout[,1]==2018,4]-GMSout[nrow(GMSout),4])/GMSout[GMSout[,1]==2018,4] #GMSout[nrow(GMSout),4]


#% reduction with vaccination
#low value parameter set
lv <- (min_default.result1-min_default.result0)*100/min_default.result0
hv <- (max_default.result1-max_default.result0)*100/max_default.result0
val.name <- c('incidence', 'incidence', 'prevalence', 'prevalence')
plotVector <- c(lv[1],hv[1],lv[2],hv[2])
names(plotVector) <- val.name

png(file=paste('result/vaccineEffect_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.25)
barplot(plotVector, main='% reduction in incidence & prevalence \n due to vaccination (EDAT+ITN+MDA  VS EDAT+ITN+MDA+RTSS)', col=c('light blue','gold','light blue','gold'),ylab='% reduction')
legend(x=3,y=3,legend=c('lower value of parameter','higher value of parameter'), fill=c('light blue','gold'), cex=1.5)
dev.off()