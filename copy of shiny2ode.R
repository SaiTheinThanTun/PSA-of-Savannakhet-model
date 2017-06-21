###################
####shiny2ode######
#Developed by Sai Thein Than Tun, sai@tropmedres.ac
#2017 June 20
#Turn Shiny UI Scripts into a matrix with
#Variable name, Default value, Min, and Max (and Steps if applicable)

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
    print(paste(varName,"is of length",length(result)))
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

shiny2ode <- function(shinyStringS){
  shinyStringS <- str_split(shinyStringS, '\n')
  trimmed <- lapply(shinyStringS, str_trim)
  y <- lapply(trimmed[[1]][str_length(trimmed[[1]])>10], odeExtractor)
  
  matrix(unlist(y), nrow=length(y), byrow=T)
}

#example  
z <- 'sliderInput(inputId="API", label = "baseline API", value = 10, min=1, max=100,step=0.5),
sliderInput(inputId="bh_max", label = "number of mosquito bites per human per night (peak season)", value = 20, min=0, max=80,step=1), 
sliderInput(inputId="eta", label = "% of all infections that are caught outside the village (forest)", value = 30, min=0, max=100,step=10),
sliderInput(inputId="covEDAT0", label = "baseline % of all clinical cases treated", value = 25, min=0, max=100)
),
column(3,
sliderInput(inputId="covITN0", label = "baseline coverage of ITN (%) ", value = 70, min=0, max=90,step=.5),
sliderInput(inputId="effITN", label = "% of infections averted due to ownership of ITN ", value = 30, min=0, max=50), 
sliderInput(inputId="covIRS0", label = "baseline coverage of IRS (%) ", value = 0, min=0, max=90,step=10),
sliderInput(inputId="effIRS", label = "% reduction in biting rate due to IRS ", value = 15, min=0, max=25,step=5)
),
column(3,
sliderInput(inputId="muC", label = "imported clinical cases per 1000 population per year ", value = 1, min=0, max=10,step=1),
sliderInput(inputId="muA", label = "imported asymptomatic microscopically detectable carriers per 1000 population per year ", value = 1, min=0, max=100,step=1),
sliderInput(inputId="muU", label = "imported asymptomatic microscopically undetectable carriers per 1000 population per year ", value = 1, min=0, max=100,step=1)
),
column(3,
sliderInput(inputId="percfail2018", label = "% of cases failing treatment in 2018 and before ", value = 5, min=0, max=100,step=5),
sliderInput(inputId="percfail2019", label = "% of cases failing treatment in 2019  ", value = 15, min=0, max=100,step=5),
sliderInput(inputId="percfail2020", label = "% of cases failing treatment in 2020 and after  ", value = 30, min=0, max=100,step=5)
)
sliderInput(inputId="EDATscale", label = "years to scale up EDAT ", value = 1, min=.25, max=3, step=.25),
sliderInput(inputId="covEDATi", label = "new % of all clinical cases treated", value = 70, min=0, max=100,step=5)
)),

sliderInput(inputId="ITNscale", label = "years to universal access to LLIN", value = 1, min=.25, max=3, step=.25),
sliderInput(inputId="covITNi", label = "new bed-net use of LLIN (%)", value = 90, min=0, max=90,step=5)
)),
sliderInput(inputId="IRSscale", label = "years to scale up IRS ", value = 1, min=.25, max=3, step=.25),
sliderInput(inputId="covIRSi", label = "new coverage of IRS (%) ", value = 90, min=0, max=90,step=5)
))
sliderInput(inputId="lossd", label = "days prophylaxis provided by the ACT", value = 30, min=15, max=30,step=1),
sliderInput(inputId="dm", label = "months to complete each round ", value = 6, min=1, max=24,step=0.5)
sliderInput(inputId="cmda_1", label = "effective population coverage of focal MDA in round 1 ", value = 50, min=0, max=100,step=10),
sliderInput(inputId="cmda_2", label = "effective population coverage of focal MDA in round 2 ", value = 50, min=0, max=100,step=10),
sliderInput(inputId="cmda_3", label = "effective population coverage of focal MDA in round 3 ", value = 50, min=0, max=100,step=10)
),

column(3,
sliderInput(inputId="tm_1", label = "timing of 1st round ", value = 9, min=1, max=36,step=1),
sliderInput(inputId="tm_2", label = "timing of 2nd round ", value = 10, min=2, max=36,step=1),
sliderInput(inputId="tm_3", label = "timing of 3rd round ", value = 11, min=3, max=36,step=1)
),
column(3,
sliderInput(inputId="effv_1", label = "% protective efficacy of RTS,S with 1st dose", value = 75, min=0, max=100),
sliderInput(inputId="effv_2", label = "% protective efficacy of RTS,S with 2nd dose", value = 80, min=0, max=100),
sliderInput(inputId="effv_3", label = "% protective efficacy of RTS,S with 3rd dose", value = 92, min=0, max=100),
sliderInput(inputId="vh", label = "half-life of vaccine protection (days)", value = 90, min=10, max=500,step=10)
)
sliderInput(inputId="MSATscale", label = "years to scale up MSAT ", value = 1, min=.25, max=3, step=.25), 
sliderInput(inputId="covMSATi", label = "new coverage of MSAT (%)", value = 90, min=0, max=100,step=10)
),
column(3,
sliderInput(inputId="MSATsensC", label = "sensitivity HS RDT (clinical) ", value = 99, min=0, max=100,step=5),
sliderInput(inputId="MSATsensA", label = "sensitivity HS RDT (micro detectable, asym)", value = 87, min=0, max=100,step=5),
sliderInput(inputId="MSATsensU", label = "sensitivity HS RDT (micro undetectable, asym)", value = 4, min=0, max=100,step=5)
)'


valueTable <- shiny2ode(z)



#draft
# sec.part <- str_split('sliderInput(inputId="API", label = "baseline API", value = 10, min=1, max=100,step=0.5),', 'inputId=')[[1]][2]
# 
# varName <- str_sub(str_split(sec.part,",.*label")[[1]][1],2,-2)
# 
# Values0 <- str_split(sec.part,",.*label")[[1]][2]
# 
# defValue.part <- str_split(Values0,",?[ ]min")[[1]][1]
# 
# Values1 <- str_split(Values0,",?[ ]min")[[1]][2]
# 
# defValue <- str_trim(str_split(defValue.part, ",.*value.*=")[[1]][2])
# 
# #another way
# 
# split.vector <- str_split('sliderInput(inputId="API", label = "baseline API", value = 10, min=1, max=100,step=0.5),', '=|,|\\)')
# step1 <- lapply(split.vector,str_trim)
# varName <- str_sub(step1[[1]][2],2,-2)
# step2 <- lapply(step1, as.numeric)
# Values <- t(sapply(step2, na.remove))
# c(varName,Values)
# x <- lapply(z, str_trim)
# y <- lapply(x[[1]][str_length(x[[1]])>10], odeExtractor)
# 
# matrix(unlist(y), nrow=4, byrow=T)
