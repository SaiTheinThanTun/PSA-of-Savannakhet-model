#no longer need to run 3 if running this
###########################################
###########################################
#setwd("~/OneDrive/MORU/Projects/PSA of Savannakhet model/truePSA/output_Inc_Prev_seas/") #mac
setwd("D:/OneDrive/MORU/Projects/PSA of Savannakhet model/truePSA/output_Inc_Prev_seas/") #windows
library(deSolve)
library(shiny)
library(TSA)
library(Rcpp)
library(stringr)
library(ggplot2)
library(reshape)

#sourceCpp("~/OneDrive/MORU/Projects/PSA of Savannakhet model/modGMS.cpp")
#source("~/OneDrive/MORU/Projects/PSA of Savannakhet model/modified copy of shiny2ode.R")
sourceCpp("0_modGMS.cpp")
source("1_modified copy of shiny2ode.R")

result <- read.csv('result/psa_run_2017-07-25 18_26_18.csv')[,-1]
# 
library(manipulate)
manipulate(
  #a <- sum((result[,2]>=prevmin&result[,2]<prevmax)&(result[,1]>=incmin&result[,1]<incmax))
  plot(sum((result[,2]>=prevmin&result[,2]<prevmax)&(result[,1]>=incmin&result[,1]<incmax)), main=paste(sum((result[,2]>=prevmin&result[,2]<prevmax)&(result[,1]>=incmin&result[,1]<incmax)))),
  prevmin=slider(2.1,4),
  prevmax=slider(4.1,5),
  incmin=slider(1,2.5),
  incmax=slider(2.5,3)
)

#chosen values
#prev 3.8 to 4.5
#inc 2.3 to 2.8
params_generated <- read.csv('result/psa_run_2017-07-25 18_26_18.csv')[,-1]
params_plausible_base <- params_generated[(params_generated[,2]>=3.48&params_generated[,2]<5)&(params_generated[,1]>=2.08&params_generated[,1]<2.814),]
write.csv(params_plausible_base,paste('result/plausibleParametersBase_',nrow(params_plausible_base),'.csv', sep = ''), row.names = FALSE)


#COMBINING RANDOM Intervention parameters

#1. getting default intervention parameters (THIS IS NOT USED IN THIS VERSION!!! JUST FOR THE NAMES)
defaultParams <- c(
  API = 10,
  bh_max = 20,
  eta = 30,
  covEDAT0 = 25,
  covITN0 = 70,
  effITN = 30,
  covIRS0 = 0,
  effIRS = 15,
  muC = 1,
  muA = 1,
  muU = 1,
  percfail2018 = 5,
  percfail2019 = 15,
  percfail2020 = 30,
  EDATscale = 1,
  covEDATi = 70,
  ITNscale = 1,
  covITNi = 90,
  IRSscale = 1,
  covIRSi = 90,
  lossd = 30,
  dm = 6,
  cmda_1 = 50,
  cmda_2 = 50,
  cmda_3 = 50,
  tm_1 = 9,
  tm_2 = 10,
  tm_3 = 11,
  effv_1 = 75,
  effv_2 = 80,
  effv_3 = 92,
  vh = 90,
  MSATscale = 1,
  covMSATi = 90,
  MSATsensC = 99,
  MSATsensA = 87,
  MSATsensU = 4 
)

length(defaultParams) #37
dim(params_plausible_base) #151 39 #additional 2 columns are incidence & prevalence values

#2. Generating random intervention parameters
#since baseline has 102 sets of parameters, we'll try to have 10200 random intervention parameters
#into which 102 sets can be repeated 100 times

#initialize for the PSA values for intervention 20170725
set.seed <- 50
noBaseline <- 10200 #will depend on the number of sets you want and the number of basline sets with plausible values (that fall into Savannakhet incidence/prevalence)
psaTableIntv <- matrix(NA, noBaseline, nrow(MinMax))
for(j in 1:noBaseline){
  for(i in 1:nrow(MinMax)){
    psaTableIntv[j,i] <- runif(1,MinMax[i,1],MinMax[i,2])
  }
}

colnames(psaTableIntv) <- names(defaultParams)
psaTableIntv_backup <- psaTableIntv

#then remodify the sequentially rising values such as effv_1, effv_2, effv_3 and tm_x parameters
####this could be very important for the vaccine effect
tm_x <- t(sapply(psaTableIntv[,28], function(y){rev(seq(y,length=3,by=-1))})) #decrease in reverse by 1 for tm_x variable [26:28]
effv_x <- t(sapply(psaTableIntv[,31], function(y){y*sort(c(1,runif(2)))})) #decrease in reverse by a proportion of effv_x variable [29:31]


part1 <- psaTableIntv[,1:25]
part2 <- psaTableIntv[,32:37]

psaTableIntv <- cbind(part1,tm_x,effv_x,part2)
colnames(psaTableIntv) <- names(defaultParams)

#3. combining random intervention parameters with generated plausible parameters for baseline
#that is Baseline+Random Intervention parameters
#baseline 1:14, intervention 15:37

baseline <- params_plausible_base[,-c(1,2,17:39)]
baseline <- baseline[rep(seq_len(nrow(baseline)), each=100),]
params_plausible_rInt <- cbind(baseline,psaTableIntv[,15:37])
# intervention <- data.frame(matrix(rep(defaultParams[15:37],nrow(baseline)),nrow(baseline), byrow = TRUE))
# names(intervention) <- names(defaultParams)[15:37]
# 
# params_plausible_rInt <- rbind(cbind(baseline,intervention)[-1,],defaultParams) #adding the default parameters at 102nd row
dim(params_plausible_rInt) #10200 37
write.csv(params_plausible_rInt, 'result/plausibleParamFinal.csv', row.names = FALSE)



#code to construct the 'for' loop for PSA will be generated from sourcing 'modified copy of shiny2ode.R'
#valueTable was also in there
resultList <- list()

#initialize/read in final parameter set 
simValueTable <- read.csv('result/plausibleParamFinal.csv')

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
  VACon <- FALSE
  MSATon <- j>1
  primon <- FALSE  
  #no.s <- 2 #no. of simulations/ no. of values between the range
  #result <- matrix(NA, nrow(valueRange)*no.s, 2) #valueTable and valueRange are from 'modified copy of shiny2ode.R'
  
  #simValueTable <- genSimValue(valueRange,no.s) 
  #valueRange was for the univariate sensitivity analysis
  #we will use here psaTableBaseline that is generated in 'modified copy of shiny2ode.R'
  #here we have to use 'plausibleParamFinal151.csv' which is initialized out of j loop
  #simValueTable <- psaTableBaseline
  result <- matrix(NA, nrow(simValueTable),2)
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
    
    result[i,1] <- GMSout[nrow(GMSout),2]
    result[i,2] <- GMSout[nrow(GMSout),4]
  }
  
  result_combined <- cbind(result,simValueTable)
  colnames(result_combined) <- c('detected_incidence','prevalence',valueTable[,1])
  write.csv(result_combined,paste("result/Intervention_",j,"_",gsub(':','_',Sys.time()),".csv", sep = ''))
  
  resultList[[j]] <- result
}

#writing the result of incidence and prevalence values after each intervention
######################################
#######for MULTIPLE interventions#####
######################################

# resultTable <- data.frame(cbind(resultList[[1]],resultList[[2]],resultList[[3]],resultList[[4]],resultList[[5]],resultList[[6]],resultList[[7]]))
# xlabels <- c("Base","EDAT","ITN","IRS","MDA","VAC","MSAT")
# names(resultTable) <- rep(xlabels, each=2)
# write.csv(resultTable,paste("result/IncPrev_result",gsub(':','_',Sys.time()),".csv", sep = ''), row.names = FALSE)
# 
# 
# #PLOT
# png(file=paste('result/IncidenceByIntervention_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
# par(cex=1.5)
# boxplot(resultTable[,seq(1,14,by=2)], main='detectedIncidence/1000/month')
# dev.off()
# 
# png(file=paste('result/PrevalenceByIntervention_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
# par(cex=1.5)
# boxplot(resultTable[,seq(2,14,by=2)], main='TruePrevalence')
# dev.off()
# 
# 
# ###########
# #INCIDENCE#
# ###########
# #reshaping data for line plots Incidence
# incId <- cbind(1:nrow(resultTable),resultTable[,seq(1,14,by=2)])
# names(incId)[1] <- 'id'
# sum(incId$MDA<incId$VAC) #17 why 17 cases in which VAC bing up the incidence
# 
# incIdMolten <- melt(incId, id='id')
# 
# #ggplot2
# #ggplot(data=data, aes(x=measurement, y=weight, col=group, group=group:rat)) + 
# #geom_point() + geom_line() + facet_wrap(~group, ncol=2)
# png(file=paste('result/IncidenceLine_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
# par(cex=1.5)
# g <- ggplot(data=incIdMolten, aes(x=variable, y=value,group=id))
# g <- g+geom_point(shape=1, fill=NA)+geom_line(colour=rainbow(nrow(incIdMolten)))+ggtitle('Incidence/1000/Month')
# g <- g+geom_point(data = incIdMolten[seq(102, by=102, length=7),], aes(shape=factor(id)))
# g
# dev.off()
# 
# #reshaping data for line plots with facets
# resultTable2 <- resultTable #[-length(resultTable),]
# incId <- cbind(1:nrow(resultTable2),rep(1:6, each=17),resultTable2[,seq(1,14,by=2)])
# names(incId)[1:2] <- c('id','group')
# sum(incId$MDA<incId$VAC) #17 why 17 cases in which VAC bing up the incidence
# 
# incIdMolten <- melt(incId, id=c('id', 'group'))
# 
# #ggplot2
# #ggplot(data=data, aes(x=measurement, y=weight, col=group, group=group:rat)) + 
# #geom_point() + geom_line() + facet_wrap(~group, ncol=2)
# png(file=paste('result/IncidenceLineFacet_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
# par(cex=1.5)
# g <- ggplot(data=incIdMolten, aes(x=variable, y=value,group=id, col=group))
# g+geom_point()+geom_line()+ggtitle('Incidence/1000/Month')+ facet_wrap(~group, ncol=2)
# #colour=rainbow(1050)
# dev.off()
# 
# ############
# #PREVALENCE#
# ############
# #reshaping data for line plots Prevalence
# prevId <- cbind(1:nrow(resultTable2),resultTable2[,seq(2,14,by=2)])
# names(prevId)[1] <- 'id'
# sum(prevId$MDA<prevId$VAC) #0 such cases
# 
# prevIdMolten <- melt(prevId, id='id')
# 
# #ggplot2
# #ggplot(data=data, aes(x=measurement, y=weight, col=group, group=group:rat)) + 
# #geom_point() + geom_line() + facet_wrap(~group, ncol=2)
# png(file=paste('result/PrevalenceLine_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
# par(cex=1.5)
# g <- ggplot(data=prevIdMolten, aes(x=variable, y=value,group=id))
# g <- g+geom_point(shape=1, fill=NA)+geom_line(colour=rainbow(nrow(prevIdMolten)))+ggtitle('Prevalence')
# g <- g+geom_point(data = prevIdMolten[seq(102, by=102, length=7),], aes(shape=factor(id)))
# g
# dev.off()
# 
# #reshaping data for line plots with facets
# resultTable2 <- resultTable #[-length(resultTable),]
# prevId <- cbind(1:nrow(resultTable2),rep(1:6, each=17),resultTable2[,seq(1,14,by=2)])
# names(prevId)[1:2] <- c('id','group')
# 
# prevIdMolten <- melt(prevId, id=c('id', 'group'))
# 
# #ggplot2
# #ggplot(data=data, aes(x=measurement, y=weight, col=group, group=group:rat)) + 
# #geom_point() + geom_line() + facet_wrap(~group, ncol=2)
# png(file=paste('result/PrevalenceLineFacet_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
# par(cex=1.5)
# g <- ggplot(data=prevIdMolten, aes(x=variable, y=value,group=id, col=group))
# g+geom_point()+geom_line()+ggtitle('Prevalence')+ facet_wrap(~group, ncol=2)
# #colour=rainbow(1050)
# dev.off()


######################################
######for MSAT off and on######
######################################
resultTable <- data.frame(cbind(resultList[[1]],resultList[[2]]))
xlabels <- c("Off","On")
names(resultTable) <- rep(xlabels, each=2)
write.csv(resultTable,paste("result/IncPrev_result",gsub(':','_',Sys.time()),".csv", sep = ''), row.names = FALSE)


#PLOT
png(file=paste('result/IncidenceByMSAT_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.5)
boxplot(resultTable[,seq(1,4,by=2)], main='detectedIncidence/1000/month')
dev.off()

png(file=paste('result/PrevalenceByMSAT_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.5)
boxplot(resultTable[,seq(2,4,by=2)], main='TruePrevalence')
dev.off()


###########
#INCIDENCE#
###########
#reshaping data for line plots Incidence
incId <- cbind(1:nrow(resultTable),resultTable[,seq(1,4,by=2)])
names(incId)[1] <- 'id'
sum(incId$MDA<incId$VAC) #17 why 17 cases in which VAC bing up the incidence

incIdMolten <- melt(incId, id='id')

#ggplot2
#ggplot(data=data, aes(x=measurement, y=weight, col=group, group=group:rat)) + 
#geom_point() + geom_line() + facet_wrap(~group, ncol=2)
png(file=paste('result/IncidenceLine_MSAT_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.5)
g <- ggplot(data=incIdMolten, aes(x=variable, y=value,group=id))
g <- g+geom_point(shape=1, fill=NA)+geom_line(colour=rainbow(nrow(incIdMolten)))+ggtitle('Incidence/1000/Month')
#g <- g+geom_point(data = incIdMolten[seq(102, by=102, length=7),], aes(shape=factor(id)))
g
dev.off()

#reshaping data for line plots with facets
resultTable2 <- resultTable #[-length(resultTable),]
incId <- cbind(1:nrow(resultTable2),rep(1:2, each=5100),resultTable2[,seq(1,4,by=2)])
names(incId)[1:2] <- c('id','group')
sum(incId$MDA<incId$VAC) #17 why 17 cases in which VAC bing up the incidence

incIdMolten <- melt(incId, id=c('id', 'group'))

#ggplot2
#ggplot(data=data, aes(x=measurement, y=weight, col=group, group=group:rat)) + 
#geom_point() + geom_line() + facet_wrap(~group, ncol=2)
png(file=paste('result/IncidenceLineFacet_MSAT_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.5)
g <- ggplot(data=incIdMolten, aes(x=variable, y=value,group=id, col=group))
g+geom_point()+geom_line()+ggtitle('Incidence/1000/Month')+ facet_wrap(~group, ncol=2)
#colour=rainbow(1050)
dev.off()

############
#PREVALENCE#
############
#reshaping data for line plots Prevalence
prevId <- cbind(1:nrow(resultTable2),resultTable2[,seq(2,4,by=2)])
names(prevId)[1] <- 'id'
sum(prevId$MDA<prevId$VAC) #0 such cases

prevIdMolten <- melt(prevId, id='id')

#ggplot2
#ggplot(data=data, aes(x=measurement, y=weight, col=group, group=group:rat)) + 
#geom_point() + geom_line() + facet_wrap(~group, ncol=2)
png(file=paste('result/PrevalenceLine_MSAT_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.5)
g <- ggplot(data=prevIdMolten, aes(x=variable, y=value,group=id))
g <- g+geom_point(shape=1, fill=NA)+geom_line(colour=rainbow(nrow(prevIdMolten)))+ggtitle('Prevalence')
#g <- g+geom_point(data = prevIdMolten[seq(102, by=102, length=7),], aes(shape=factor(id)))
g
dev.off()

#reshaping data for line plots with facets
resultTable2 <- resultTable #[-length(resultTable),]
prevId <- cbind(1:nrow(resultTable2),rep(1:2, each=5100),resultTable2[,seq(1,4,by=2)])
names(prevId)[1:2] <- c('id','group')

prevIdMolten <- melt(prevId, id=c('id', 'group'))

#ggplot2
#ggplot(data=data, aes(x=measurement, y=weight, col=group, group=group:rat)) + 
#geom_point() + geom_line() + facet_wrap(~group, ncol=2)
png(file=paste('result/PrevalenceLineFacet_MSAT_',gsub(':','_',Sys.time()),'.png',sep = ''), width = 1280, height = 800)
par(cex=1.5)
g <- ggplot(data=prevIdMolten, aes(x=variable, y=value,group=id, col=group))
g+geom_point()+geom_line()+ggtitle('Prevalence')+ facet_wrap(~group, ncol=2)
#colour=rainbow(1050)
dev.off()