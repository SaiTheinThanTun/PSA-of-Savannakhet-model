#To explore the following 3 scenarios####
#1. EDAT+ITN+MDA
#2. EDAT+ITN+MDA+RTS same coverage
#3. EDAT+ITN+MDA+RTS full coverage

#write.csv(x,'xSeries.csv',row.names = FALSE) #to initiate the plot
setwd('/Users/sai/OneDrive/MORU/Projects/PSA of Savannakhet model/SA_v9/USA_v9')
source('runGMSFunction.R')
xSeries <- read.csv('xSeries.csv')[,1]

#loop values####
bh_max_i <- c(18,23,25)
VACon_j <- c(0,1,1)
v_same_j <- c(0,1,0)
scen_labels <- c('EDAT+ITN+MDA','EDAT+ITN+MDA+RTS same coverage','EDAT+ITN+MDA+RTS full coverage')

#fixed parameters####
EDATon <- TRUE
ITNon <- TRUE
IRSon <- FALSE
MDAon <- TRUE
MSATon <- TRUE
primon <- FALSE



API <- 10
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



#plot setting####
#choose incidence or prevalence####
#png('3Scenarios_incidence.png', width = 1800, height=1600)
png('3Scenarios_prevalence.png', width = 1800, height=1600)
par(mfrow=c(3,1), cex=1.5)

for(i in 1:3){
  #choose incidence or prevalence####
  #plot(xSeries,rep(0,length(xSeries)),ylim=c(0,20),type = 'n',xlab='Time', ylab = 'Incidence')
  plot(xSeries,rep(0,length(xSeries)),ylim=c(0,38),type = 'n',xlab='Time', ylab = 'Prevalence')
  legend(2021,20,legend=scen_labels,lty=1:3)
  for(j in 1:3){
    #print(bh_max_i[i])
    #print(VACon_j[j])
    # print(v_same_j[j])
    
    ####intervention switches####
    
    VACon <- VACon_j[j] #TRUE
    v_same <- v_same_j[j] #TRUE
    
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
    
    bh_max <- bh_max_i[i] #20
    
    
    #####default run####
    
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
    
    #plot here
    #lty=j for plots and legends
    
    times<-GMSout[,1]
    clinmonth_det<-GMSout[,2] #cbind(GMSout0[,2],GMSouti[,2])
    clinmonth_tot<-GMSout[,3] #cbind(GMSout0[,3],GMSouti[,3])
    prevalence<-GMSout[,4] #cbind(GMSout0[,4],GMSouti[,4])
    
    runin<-(2016-startyear)/dt
    
    finclin<-max(clinmonth_tot[(runin:length(clinmonth_det))])
    finprev<-max(prevalence[(runin:length(prevalence))])
    
    
    # PLOTTING
    #par(mfrow=c(1,2), cex=1.5)
    
    maxy<-max(finclin,API/12)
    x<-times[(runin:length(clinmonth_det))]
    
    #for incidence, turn on the 2 lines here####
    # y1<-clinmonth_det[runin:length(clinmonth_det)]
    # lines(x,y1,type = 'l', lty=j)
    
    #for prevalence, turn on the 2 lines here####
    y2<-prevalence[runin:length(prevalence)]
    lines(x,y2,type = 'l', lty=j)
    
  }
}
dev.off()
par(mfrow=c(1,1), cex=1)