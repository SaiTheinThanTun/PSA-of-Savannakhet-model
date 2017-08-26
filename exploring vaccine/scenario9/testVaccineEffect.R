#initialize/read in final parameter set 
setwd('/Users/sai/OneDrive/MORU/Projects/PSA of Savannakhet model/exploring vaccine/scenario9/')

#code to construct the 'for' loop for PSA will be generated from sourcing 'modified copy of shiny2ode.R'
#valueTable was also in there
library(ggplot2)
library(stringr)
library(plyr)
library(deSolve)
library(Rcpp)
sourceCpp('modGMS.cpp')


resultList <- list()


a <- getwd()
foldername <- str_split(a,'/')[[1]][length(str_split(a,'/')[[1]])]
#simValueTable <- read.csv('plausibleParamFinal_0MDAcoverage.csv', header=TRUE)
#simValueTable <- read.csv('result/plausibleParamFinal_0MDAcoverage.csv')
var37 <- read.csv('var37.csv', stringsAsFactors = FALSE)[,2]
test <- read.csv('MalMod_2017-08-02 203845.csv', stringsAsFactors = FALSE)
test[8,4] <- 'API'
test <- test[-(1:7),]

test$VarName <- factor(test$VarName, levels=var37)
#levels(test$VarName) <- var37
test <- test[order(test$VarName),]

simValueTable <- as.data.frame(rbind(test[,2]))
#simValueTable <- as.data.frame(rbind(test[,2],test[,2]))
colnames(simValueTable) <- test[,4]
  
# for(j in 1:7){
#   EDATon <- j>1
#   ITNon <- j>2
#   IRSon <- j>3
#   MDAon <- j>4
#   VACon <- j>5
#   MSATon <- j>6
#   primon <- FALSE
for(j in 1:2){
  EDATon <- TRUE
  ITNon <- TRUE
  IRSon <- TRUE
  MDAon <- TRUE
  VACon <- j>1
  MSATon <- FALSE
  primon <- FALSE  
  #no.s <- 2 #no. of simulations/ no. of values between the range
  #result <- matrix(NA, nrow(valueRange)*no.s, 2) #valueTable and valueRange are from 'modified copy of shiny2ode.R'
  
  #simValueTable <- genSimValue(valueRange,no.s) 
  #valueRange was for the univariate sensitivity analysis
  #we will use here psaTableBaseline that is generated in 'modified copy of shiny2ode.R'
  #here we have to use 'plausibleParamFinal151.csv' which is initialized out of j loop
  #simValueTable <- psaTableBaseline
  result <- matrix(NA, nrow(simValueTable),6)
  for(i in 1:nrow(simValueTable)){
    #input (later will serve as getting data from for loop)
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
    cmda_1 <- simValueTable[i,23]
    cmda_2 <- simValueTable[i,24]
    cmda_3 <- simValueTable[i,25]
    tm_1 <- simValueTable[i,26]
    tm_2 <- simValueTable[i,27]
    tm_3 <- simValueTable[i,28]
    effv_1 <- simValueTable[i,29]
    effv_2 <- simValueTable[i,30]
    effv_3 <- simValueTable[i,31]
    vh <- simValueTable[i,32]
    MSATscale <- simValueTable[i,33]
    covMSATi <- simValueTable[i,34]
    MSATsensC <- simValueTable[i,35]
    MSATsensA <- simValueTable[i,36]
    MSATsensU <- simValueTable[i,37]
    
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
    
    result[i,1] <- GMSout[nrow(GMSout),2] #incidence at year 2023
    result[i,2] <- GMSout[nrow(GMSout),4] #prevalence at year 2023
    result[i,3] <- min(GMSout[GMSout[,1]>2018,2]) #minimal incidence
    result[i,4] <- min(GMSout[GMSout[,1]>2018,4]) #minimal prevalence
    result[i,5] <- GMSout[GMSout[,1]>2018,1][which.min(GMSout[GMSout[,1]>2018,2])] #time at min incidence
    result[i,6] <- GMSout[GMSout[,1]>2018,1][which.min(GMSout[GMSout[,1]>2018,4])] #time at min prevalence
  }
  
  #result_combined <- cbind(result,simValueTable)
  #colnames(result_combined) <- c('detected_incidence','prevalence',colnames(simValueTable)) #valueTable[,1])
  #write.csv(result_combined,paste(foldername,"_Intervention_",j,"_",gsub(':','_',Sys.time()),".csv", sep = ''))
  
  resultList[[j]] <- result
}


######################################
######for VACCINATION off and on######
######################################
resultTable <- data.frame(cbind(resultList[[1]],resultList[[2]]))
xlabels <- c("Off","On")
names(resultTable) <- paste(rep(xlabels, each=6),c('inc','prev','min.inc','min.prev','min.inc.t','min.prev.t'))

resultTable <- mutate(resultTable, incidence.red=resultTable[,1]-resultTable[,7])
resultTable <- mutate(resultTable, prevalence.red=resultTable[,2]-resultTable[,8])
write.csv(resultTable,paste(foldername,"_IncPrev_result",gsub(':','_',Sys.time()),".csv", sep = ''), row.names = FALSE)


#PLOT
png(file=paste(foldername,'_IncidenceByVACCINATION_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.5)
boxplot(resultTable[,c(1,7,3,9)], main=paste('detectedIncidence/1000/month\n MinValue@\n',resultTable[,5],' vs ', resultTable[,11],sep=''))
dev.off()

png(file=paste(foldername,'_PrevalenceByVACCINATION_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.5)
boxplot(resultTable[,c(2,8,4,10)], main=paste('TruePrevalence\n MinValue@\n',resultTable[,6],' vs ',resultTable[,12],sep=''))
dev.off()
