############################################################
####Univariate Sensitivity Analysis: Scenario analysis######
#Developed by Sai Thein Than Tun, sai@tropmedres.ac
#First development date: 2017 June 20
#Updated: 2017 September 08
#run the model for default value and then min and max values varying each value at a time

setwd("D:\\OneDrive\\MORU\\Projects\\PSA of Savannakhet model\\SA_v9\\USA_v9")
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

####arranging results with vaccine off####
incidenceT0 <- matrix(result[,1], nrow(result)/no.s,no.s, byrow = T)
prevalenceT0 <- matrix(result[,2], nrow(result)/no.s,no.s, byrow = T)




###Scenario 2####
source('1_shiny2ode.R')
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

#non-reactive function runGMS is now outside of the server function
runGMS<-function(initprev, scenario, param) 
{
  #MODEL PARAMETERS
  parameters <- c(scenario,
                  timei = 2018,
                  nuTr = 14,                   # days of infectiosness after treatment ACT [N]
                  nuTrp = 7,                   # days of infectiosness after treatment ACT+primaquine [N]
                  alpha = 0.7,                   # relative amplitude seasonality [N]
                  phi = 0.0,                   # phase angle seasonality [N]
                  epsilonh=0.23,                 # per bite probability of an infectious mosquito infecting a human
                  epsilonm=0.5,                  # per bite probability of an infectious human infecting a mosquito
                  b=365/3,                       # per mosquito rate of biting
                  deltam=365/14,                 #
                  gammam=365/10,#Rate of becoming infectious from the latent phase for mosquitos, Kai Matuschewski: Getting infectious
                  cm_1=80,
                  cm_2=95,
                  cm_3=95,
                  covMSAT0=0,
                  omega = 2,                   # average duration of immunity (years) [N]
                  nuC = 3,                     # days of symptoms in the absence of treatment [N], #change 9 -> 3
                  nuA = 60,                    # days of asymptomatic microscopically detectable carriage [N]
                  nuU = 100,                    # days of asymptomatic microscopically undetectable carriage [N], #change 60 -> 100, Mean duration of a malaria untreated infection: 160 days, 
                  rhoa = 55,                   # relative infectivity of asymptomatic microscopically detectable carriers compared with clinical infections (%) [N]
                  rhou = 17,                   # relative infectivity of asymptomatic microscopically undetectable carriers compared with clinical infections (%) [N]
                  ps = 90,                     # % of all non-immune new infections that are clinical [N]
                  pr = 20,                     # % of all immune new infections that are clinical [N]
                  mu = 50,                      # life expectancy (years) [N]
                  param)
  
  
  
  # MODEL INITIAL CONDITIONS
  # population size
  initP<-10000 
  
  initS_0<-0.5*(1-initprev)*initP
  initIC_0<-0
  initIA_0<-initprev*initP
  initIU_0<-0
  initR_0<-0.5*(1-initprev)*initP
  initTr_0<-0
  
  state <- c(Y = 0, Cinc_det = 0, Cinc_tot = 0, 
             S_0 = initS_0, IC_0 = initIC_0, IA_0 = initIA_0, IU_0 = initIU_0, R_0 = initR_0, Tr_0 = initTr_0, Sm_0 = 0, Rm_0 = 0,
             S_1 = 0, IC_1 = 0, IA_1 = 0, IU_1 = 0, R_1 = 0, Tr_1 = 0, Sm_1 = 0, Rm_1 = 0,
             S_2 = 0, IC_2 = 0, IA_2 = 0, IU_2 = 0, R_2 = 0, Tr_2 = 0, Sm_2 = 0, Rm_2 = 0,
             S_3 = 0, IC_3 = 0, IA_3 = 0, IU_3 = 0, R_3 = 0, Tr_3 = 0, Sm_3 = 0, Rm_3 = 0,
             S_4 = 0, IC_4 = 0, IA_4 = 0, IU_4 = 0, R_4 = 0, Tr_4 = 0, Sm_4 = 0, Rm_4 = 0,
             y01_1=as.vector(param["effv_3"])*.008, y02_1 = 0, yf_1=0, ys_1=0,
             y01_2=as.vector(param["effv_3"])*.009, y02_2 = 0, yf_2=0, ys_2=0,
             y01_3=as.vector(param["effv_3"])/100, y02_3 = 0, yf_3=0, ys_3=0,
             y01_4=as.vector(param["effv_4"])/100, y02_4 = 0, yf_4=0, ys_4=0
  )
  
  
  #out <- ode(y = state, times = times, func = modGMS, parms = parameters)
  WmodGMSrcpp<-function(t,state,parameters){
    tmp<-modGMSrcpp(t,state,parameters)
    return(list(tmp))
  }
  #out <- ode(y = state, times = times, func = WmodGMSrcpp, parms = parameters)
  out <- ode(y = state, times = times, func = WmodGMSrcpp, parms = parameters, method="vode")
  
  # MODEL OUTPUTS
  ipop <- 5:44
  iinc_det <- 3
  iinc_tot <- 4
  iprev <- c(6,  7,  8, 10, 14, 15, 16, 18, 22, 23, 24, 26, 30, 31, 32, 34, 38, 39, 40, 42)
  
  # population
  times<-out[,1]+startyear
  pop<-rowSums(out[,ipop])
  
  
  # clinical incidence detected per 1000 per month
  tci_det <- out[,iinc_det]
  clinmonth_det <- tci_det
  clinmonth_det[1] <- 0
  clinmonth_det[2:length(times)] <- 1000*(tci_det[2:length(times)] - tci_det[1:(length(times)-1)])/pop[2:length(times)]
  
  # clinical incidence total per 1000 per month
  tci_tot <- out[,iinc_tot]
  clinmonth_tot <- tci_tot
  clinmonth_tot[1] <- 0
  clinmonth_tot[2:length(times)] <- 1000*(tci_tot[2:length(times)] - tci_tot[1:(length(times)-1)])/pop[2:length(times)]
  
  
  # % prevalence
  prevalence <- 100*rowSums(out[,iprev])/pop # Additional file: Equation no.13
  GMSout<-matrix(NA,nrow=length(times),ncol=4)
  GMSout[,1]<-times
  GMSout[,2]<-clinmonth_det
  GMSout[,3]<-clinmonth_tot
  GMSout[,4]<-prevalence
  
  return(GMSout)
}

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
incidenceT1 <- matrix(result[,1], nrow(result)/no.s,no.s, byrow = T)
prevalenceT1 <- matrix(result[,2], nrow(result)/no.s,no.s, byrow = T)


###write scenario comparison results####
incidenceT <- rbind(incidenceT0, incidenceT1)
prevalenceT <- rbind(prevalenceT0, prevalenceT1)

incidenceDiff <- incidenceT-default.result[,1]
incidenceLowValue <- incidenceDiff[,1]
incidenceHiValue <- incidenceDiff[,2]

incidenceRange <- cbind(as.data.frame(cbind(incidenceLowValue,incidenceHiValue)), c(paste(valueTable[par_vaccine,1],'V_0',sep='_'),paste(valueTable[par_vaccine,1],'V_1',sep='_')))
rangeSize <- abs(incidenceRange[,2]-incidenceRange[,1])
incidenceRange <- incidenceRange[order(rangeSize),]


prevalenceDiff <- prevalenceT-default.result[,2]
prevalenceLowValue <- prevalenceDiff[,1]
prevalenceHiValue <- prevalenceDiff[,2]

prevalenceRange <- cbind(as.data.frame(cbind(prevalenceLowValue,prevalenceHiValue)), c(paste(valueTable[par_vaccine,1],'V_0',sep='_'),paste(valueTable[par_vaccine,1],'V_1',sep='_')))
rangeSize2 <- abs(prevalenceRange[,2]-prevalenceRange[,1])
prevalenceRange <- prevalenceRange[order(rangeSize2),]
# result_combined <- cbind(result,simValueTable)
# colnames(result_combined) <- c('detected_incidence','prevalence',valueTable[,1])
# write.csv(result_combined,paste("result/usa_run_",gsub(':','_',Sys.time()),".csv", sep = ''))
write.csv(incidenceRange,'result/incidenceRange.csv')
write.csv(prevalenceRange,'result/prevalenceRange.csv')
###plot####

png(file=paste('result/compare_incidenceTornado_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 1800)
par(mar=c(5, 8, 4, 2), cex=1.25)
barplot(incidenceRange[,1], names.arg = incidenceRange[,3], las=1, horiz=T, xlim=c(-.85,0.1), beside = T, axes=F, col='light blue', main=paste("Sensitivity on incidence, Savannakhet model"), xlab='% reduction in incidence')
barplot(incidenceRange[,2], horiz=T, axes=F, beside=T, add=T, col=adjustcolor( "gold", alpha.f = 0.7)) #col='gold', alpha=.8)
axis(1, at=seq(-.85,0.1,by=.1), labels = seq(-.85,0.1,by=.1)+round(default.result[,1],2), ylab='% reduction in incidence')
legend(x=-.6,y=8,legend=c('lower value of parameter','higher value of parameter'), fill=c('light blue','gold'), cex=1.5)
dev.off()




png(file=paste('result/compare_prevalenceTornado_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 1800)
par(mar=c(5, 8, 4, 2), cex=1.25)
barplot(prevalenceRange[,1], names.arg = prevalenceRange[,3], axes=F,las=1, horiz=T, xlim=c(-.2,.05), beside = T, col='light blue', main=paste("Sensitivity on prevalence, Savannakhet model"), xlab='% reduction in prevalence')
barplot(prevalenceRange[,2], horiz=T, axes=F, beside=T, add=T, col=adjustcolor( "gold", alpha.f = 0.7))
axis(1, at=seq(-.2,.05,by=.01), labels = seq(-.2,.05,by=.01)+round(default.result[,2],2), ylab='% reduction in prevalence')
legend(x=-.175,y=8,legend=c('lower value of parameter','higher value of parameter'), fill=c('light blue','gold'), cex=1.5)
dev.off()