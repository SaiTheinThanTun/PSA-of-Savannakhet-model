setwd("~/OneDrive/MORU/Projects/PSA of Savannakhet model/")
library(deSolve)
library(shiny)
library(TSA)
library(Rcpp)
sourceCpp("~/OneDrive/MORU/Projects/PSA of Savannakhet model/modGMS.cpp")
source("~/OneDrive/MORU/Projects/PSA of Savannakhet model/modified copy of shiny2ode.R")

#code to construct the 'for' loop for PSA will be generated from sourcing 'modified copy of shiny2ode.R'
#valueTable was also in there

EDATon <- TRUE
ITNon <- TRUE
IRSon <- TRUE
MDAon <- TRUE
VACon <- TRUE
MSATon <- TRUE
primon <- FALSE

#no.s <- 4 #no. of simulations/ no. of values between the range
default.result <- matrix(NA, 1, 2) #valueTable and valueRange are from 'modified copy of shiny2ode.R'

simValueTable <- genSimValue(valueRange,no.s)

  #input (later will serve as getting data from for loop)
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
cmda_1 <- 50
cmda_2 <- 50
cmda_3 <- 50
tm_1 <- 9
tm_2 <- 10
tm_3 <- 11
effv_1 <- 75
effv_2 <- 80
effv_3 <- 92
vh <- 90
MSATscale <- 1
covMSATi <- 90
MSATsensC <- 99
MSATsensA <- 87
MSATsensU <- 4  
  #non-reactive parameters
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
                    # phi = 0.0,                   # phase angle seasonality [N]
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
               S_4 = 0, IC_4 = 0, IA_4 = 0, IU_4 = 0, R_4 = 0, Tr_4 = 0, Sm_4 = 0, Rm_4 = 0
    )
    
    
    #out <- ode(y = state, times = times, func = modGMS, parms = parameters)
    WmodGMSrcpp<-function(t,state,parameters){
      tmp<-modGMSrcpp(t,state,parameters)
      return(list(tmp))
    }
    out <- ode(y = state, times = times, func = WmodGMSrcpp, parms = parameters)
    
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
  
  
  scenario_0<-c(EDATon = 0,
                ITNon = 0,
                IRSon = 0,
                MDAon = 0,
                primon = 0,
                MSATon = 0,
                VACon = 0)
  
  scenario_iR<-c(EDATon = EDATon,
                 ITNon = ITNon,
                 IRSon = IRSon,
                 MDAon = MDAon,
                 primon = primon,
                 MSATon = MSATon,
                 VACon = as.numeric(VACon))
  
  parametersR <- c(
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
    
    effv_1 = effv_1,
    effv_2 = effv_2,
    effv_3 = effv_3,
    vh = vh
  )
  
  
  # initial prevalence
  initprevR <- (0.001*API)
  
  GMSout <- runGMS(initprevR, scenario_iR,parametersR)
  
  #GMSout[nrow(GMSout),c(3,4)] #output total incidence and prevelance
  
  default.result[1,1] <- GMSout[nrow(GMSout),3]
  default.result[1,2] <- GMSout[nrow(GMSout),4]
