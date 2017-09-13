###################
####shiny2ode######
#Developed by Sai Thein Than Tun, sai@tropmedres.ac
#First development date: 2017 June 20
#Updated: 2017 September 05
#Turn Shiny UI Scripts into a matrix with
#Variable name, Default value, Min, and Max (and Steps if applicable)

setwd("D:\\OneDrive\\MORU\\Projects\\PSA of Savannakhet model\\SA_v9\\USA_v9")
library(stringr)
library(TSA)

#############################################################################################
####Parameters from shiny app################################################################
#############################################################################################
#change 1
z <- '
sliderInput(inputId="API", label = "baseline API", value = 10, min=1, max=30,step=0.5),
                    sliderInput(inputId="bh_max", label = "number of mosquito bites per human per night (peak season)", value = 20, min=15, max=25,step=1), 
sliderInput(inputId="eta", label = "% of all infections that are caught outside the village (forest)", value = 30, min=5, max=60,step=10),
sliderInput(inputId="covEDAT0", label = "baseline % of all clinical cases treated", value = 25, min=10, max=50)

       sliderInput(inputId="covITN0", label = "baseline coverage of ITN (%) ", value = 70, min=50, max=90,step=.5),
       sliderInput(inputId="effITN", label = "% of infections averted due to ownership of ITN ", value = 30, min=20, max=50), 
       sliderInput(inputId="covIRS0", label = "baseline coverage of IRS (%) ", value = 0, min=0, max=90,step=10),
       sliderInput(inputId="effIRS", label = "% reduction in biting rate due to IRS ", value = 15, min=5, max=25,step=5)

       sliderInput(inputId="muC", label = "imported clinical cases per 1000 population per year ", value = 1, min=0.5, max=10,step=1),
       sliderInput(inputId="muA", label = "imported asymptomatic microscopically detectable carriers per 1000 population per year ", value = 1, min=0.5, max=20,step=1),
       sliderInput(inputId="muU", label = "imported asymptomatic microscopically undetectable carriers per 1000 population per year ", value = 1, min=0.5, max=20,step=1)

       sliderInput(inputId="percfail2018", label = "% of cases failing treatment in 2018 and before ", value = 5, min=2, max=10,step=5),
       sliderInput(inputId="percfail2019", label = "% of cases failing treatment in 2019  ", value = 15, min=7, max=20,step=5),
       sliderInput(inputId="percfail2020", label = "% of cases failing treatment in 2020 and after  ", value = 30, min=15, max=45,step=5)

                  sliderInput(inputId="EDATscale", label = "years to scale up EDAT ", value = 1, min=.25, max=3, step=.25),
                  sliderInput(inputId="covEDATi", label = "new % of all clinical cases treated", value = 70, min=60, max=95,step=5)

           sliderInput(inputId="ITNscale", label = "years to universal access to LLIN", value = 1, min=.25, max=3, step=.25),
           sliderInput(inputId="covITNi", label = "new bed-net use of LLIN (%)", value = 90, min=70, max=90,step=5)

           sliderInput(inputId="IRSscale", label = "years to scale up IRS ", value = 1, min=.25, max=3, step=.25),
           sliderInput(inputId="covIRSi", label = "new coverage of IRS (%) ", value = 90, min=70, max=90,step=5)

                sliderInput(inputId="lossd", label = "days prophylaxis provided by the ACT", value = 30, min=15, max=30,step=1),
                sliderInput(inputId="dm", label = "months to complete each round ", value = 6, min=1, max=12,step=0.5),
                sliderInput(inputId = "ka", label = "weeks to full protective effect of vaccine after dose", value = 2, min=1, max=8),
                sliderInput(inputId = "delta", label = "months to transition between fast and slow component of vaccine", value = 5, min=1, max=12)

                sliderInput(inputId="cmda_1", label = "effective population coverage of focal MDA in round 1 ", value = 50, min=30, max=70,step=10),
                sliderInput(inputId="cmda_2", label = "effective population coverage of focal MDA in round 2 ", value = 50, min=30, max=70,step=10),
                sliderInput(inputId="cmda_3", label = "effective population coverage of focal MDA in round 3 ", value = 50, min=30, max=70,step=10),
                sliderInput(inputId = "kf", label = "days with fast component of vaccine", value = 45, min=15, max=60),
                sliderInput(inputId = "ks", label = "days with slow component of vaccine", value = 634, min=500, max=800)

                sliderInput(inputId="tm_1", label = "timing of 1st round [2018+ no. of month]", value = 9, min=9, max=9,step=1),
                sliderInput(inputId="tm_2", label = "timing of 2nd round [2018+ no. of month]", value = 10, min=10, max=10,step=1),
                sliderInput(inputId="tm_3", label = "timing of 3rd round [2018+ no. of month]", value = 11, min=11, max=11,step=1)

                sliderInput(inputId="effv_3", label = "% protective efficacy of RTS,S with 3rd dose", value = 92, min=70, max=99),
                sliderInput(inputId="effv_4", label = "% protective efficacy of RTS,S with Booster dose", value = 92, min=70, max=99)

                sliderInput(inputId="MSATscale", label = "years to scale up MSAT ", value = 1, min=.25, max=3, step=.25), 
                sliderInput(inputId="covMSATi", label = "new coverage of MSAT (%)", value = 90, min=70, max=95,step=10)

                sliderInput(inputId="MSATsensC", label = "sensitivity HS RDT (clinical) ", value = 99, min=90, max=100,step=5),
                sliderInput(inputId="MSATsensA", label = "sensitivity HS RDT (micro detectable, asym)", value = 87, min=70, max=95,step=5),
                sliderInput(inputId="MSATsensU", label = "sensitivity HS RDT (micro undetectable, asym)", value = 4, min=1, max=30,step=5)
'

#############################################################################################
####1. Grab and extract relevent numbers from "slider inputs" from shiny code (in text)######
####also tests the number of values during the process#######################################
#############################################################################################

odeExtractor <- function(shinyString){
  split.vector <- str_split(shinyString, '=|,|\\)')
  step1 <- lapply(split.vector,str_trim)
  varName <- str_sub(step1[[1]][2],2,-2)
  step2 <- lapply(step1, as.numeric)
  Values <- t(sapply(step2, na.remove))
  result <- c(varName,Values)
  
  if(length(result)==5){
    return(result)
  } else if(length(result)<5){
    print(paste(varName,"is of length",length(result), "probably missing the steps in the slider"))
    while(length(result)<5){
      result <- c(result,NA)
    }
    # for(i in 1:(5-length(result))){
    #   result <- c(result,NA)
    # }
    # if(length(result)<5){
    #   result <- c(result,NA)
    # }
    return(result)
  } else if(length(result)>5){
    print(paste(varName,"is of length",length(result)))
    return(rep(9999,5))
  }
}

#############################################################################################
####2. Grab output from odeExtractor in 1, construct default, min, max table#################
####also turns the numbers into ODE input formats############################################
#############################################################################################

shiny2ode <- function(shinyStringS){
  shinyStringS <- str_split(shinyStringS, '\n')
  trimmed <- lapply(shinyStringS, str_trim)
  y <- lapply(trimmed[[1]][str_length(trimmed[[1]])>10], odeExtractor)
  
  matrix(unlist(y), nrow=length(y), byrow=T)
}

#############################################################################################
####Running the 2 functions above############################################################
#############################################################################################

valueTable <- shiny2ode(z)

#output for reference
valueTableWrite <- valueTable[,1:4]
colnames(valueTableWrite) <- c('varName','default','min','max')
write.csv(valueTableWrite,paste("Data file ",gsub(':','_',Sys.time()),".csv", sep=""))


valueRange <- matrix(as.numeric(valueTable[,c(2,3,4)]),nrow(valueTable),3)

MinMax <- valueRange[,-1]

#transform the valueTable into 
#1. ODE code
cat(paste(valueTable[,1],"<-",valueTable[,2]), sep='\n')


#2. to use in for loop
cat(paste(valueTable[,1]," <- simValueTable[i,",1:nrow(valueTable),"]", sep=""), sep='\n')


#############################################################################################
####3. Construct a stack of parameters table from min, max value table#######################
#############################################################################################

genSimValue <- function(vrl, no.s){ #vlr Value Range Input, no.s # of simulation/ #of values between the range
  simValueTable <- matrix(NA,nrow(vrl)*no.s,nrow(vrl)) #has row as scenario and column as variables
  for(i in 1:nrow(vrl)){
    for(j in 1:(no.s)){
      simValueTable[j+(no.s*(i-1)),i] <- seq(vrl[i,2],vrl[i,3],length=no.s)[j]
    }
    simValueTable[is.na(simValueTable[,i]),i] <- vrl[i,1]
  }
  simValueTable
}

no.s <- 2 #no. of simulations/ no. of values between the range

simValueTable <- genSimValue(valueRange,no.s)



#############################################################################################
####X. initialize for the PSA values for baseline############################################
#############################################################################################
#
# set.seed <- 501
# noBaseline <- 10000
# psaTableBaseline <- matrix(NA, noBaseline, nrow(MinMax))
# for(j in 1:noBaseline){
#   for(i in 1:nrow(MinMax)){
#     psaTableBaseline[j,i] <- runif(1,MinMax[i,1],MinMax[i,2])
#   }
# }








