#runGMSFunction
library(deSolve)
library(shiny)
library(TSA)
library(Rcpp)
setwd("D:\\OneDrive\\MORU\\Projects\\PSA of Savannakhet model\\SA_v9\\USA_v9")
sourceCpp('modGMS.cpp')


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