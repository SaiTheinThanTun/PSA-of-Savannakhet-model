library(deSolve)
library(shiny)
library(TSA)
library(Rcpp)
sourceCpp("modGMS.cpp")


ui <- fluidPage(
  tags$head(includeScript("google-analytics.js")),
  tabsetPanel(
    id="panels",
    tabPanel(title = strong("Baseline"),
             column(3,
                    sliderInput(inputId="API", label = "baseline API", value = 10, min=1, max=100,step=0.5),
                    sliderInput(inputId="bh_max", label = "number of mosquito bites per human per night (peak season)", value = 20, min=0, max=80,step=1), #change range 0-80, Dan data
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
    ),
    
    tabPanel(title = strong("Interventions currently available"),
             column(4,
                    wellPanel(
                      h3("Early Diagnosis and Treatment"),
                      checkboxInput(inputId="EDATon", label = "switch on scale up of EDAT ", value = FALSE),
                      checkboxInput(inputId="primon", label = "ACT+primaquine for EDAT and MDA ", value = FALSE), #under EDAT checkbox
                      sliderInput(inputId="EDATscale", label = "years to scale up EDAT ", value = 1, min=.25, max=3, step=.25),
                      sliderInput(inputId="covEDATi", label = "new % of all clinical cases treated", value = 70, min=0, max=100,step=5)
                    )), 
             column(4,wellPanel(
                      h3("Insecticide Treated Net (LLIN)"),
                      checkboxInput(inputId="ITNon", label = "switch on scale up of LLIN", value = FALSE),
                      sliderInput(inputId="ITNscale", label = "years to universal access to LLIN", value = 1, min=.25, max=3, step=.25),
                      sliderInput(inputId="covITNi", label = "new bed-net use of LLIN (%)", value = 90, min=0, max=90,step=5)
                    )),
             column(4,wellPanel(
               h3("Indoor Residual Spray"),
               checkboxInput(inputId="IRSon", label = "switch on scale up of IRS ", value = FALSE),
               sliderInput(inputId="IRSscale", label = "years to scale up IRS ", value = 1, min=.25, max=3, step=.25),
               sliderInput(inputId="covIRSi", label = "new coverage of IRS (%) ", value = 90, min=0, max=90,step=5)
             ))
  ),
  tabPanel(title = strong("Interventions under trial: Focal MVDA (hotspot)"),
                     column(3,
                            checkboxInput(inputId="MDAon", label = "switch on MDA", value = FALSE), #6
                            sliderInput(inputId="lossd", label = "days prophylaxis provided by the ACT", value = 30, min=15, max=30,step=1),
                            sliderInput(inputId="dm", label = "months to complete each round ", value = 6, min=1, max=24,step=0.5),
                            sliderInput(inputId = "ka", label = "weeks to full protective effect of vaccine after dose", value = 2, min=1, max=8),
                            sliderInput(inputId = "delta", label = "months to transition between fast and slow component of vaccine", value = 5, min=1, max=12)
                     ),
                     column(3,
                            sliderInput(inputId="cmda_1", label = "effective population coverage of focal MDA in round 1 ", value = 50, min=0, max=100,step=10),
                            sliderInput(inputId="cmda_2", label = "effective population coverage of focal MDA in round 2 ", value = 50, min=0, max=100,step=10),
                            sliderInput(inputId="cmda_3", label = "effective population coverage of focal MDA in round 3 ", value = 50, min=0, max=100,step=10),
                            sliderInput(inputId = "kf", label = "days with fast component of vaccine", value = 45, min=15, max=60),
                            sliderInput(inputId = "ks", label = "days with slow component of vaccine", value = 634, min=500, max=800)
                     ),
                     
                     column(3,
                            sliderInput(inputId="tm_1", label = "timing of 1st round [2018+ no. of month, 1 means Jan'2018, 13 means Jan'2019]", value = 9, min=1, max=36,step=1),
                            sliderInput(inputId="tm_2", label = "timing of 2nd round [2018+ no. of month]", value = 10, min=2, max=36,step=1),
                            sliderInput(inputId="tm_3", label = "timing of 3rd round [2018+ no. of month]", value = 11, min=3, max=36,step=1)
                     ),
                     column(3,
                            radioButtons(inputId="VACon", label = "With vaccination: ", choices = c("No"=0, "Yes"=1), selected = 0, inline=TRUE),
                            checkboxInput(inputId="v_same", label = "Vaccinate at same coverage as MDA", value = TRUE), #1/Vaccine at same coverage as MDA
                            sliderInput(inputId="effv_3", label = "% protective efficacy of RTS,S with 3rd dose", value = 92, min=0, max=100),
                            sliderInput(inputId="effv_4", label = "% protective efficacy of RTS,S with Booster dose", value = 92, min=0, max=100)
                     )
            ),
            tabPanel(title = strong("Interventions under trial: Focal MSAT (mobile)"),
                     column(3,
                            checkboxInput(inputId="MSATon", label = "switch on MSAT for imported cases", value = FALSE),
                            sliderInput(inputId="MSATscale", label = "years to scale up MSAT ", value = 1, min=.25, max=3, step=.25), 
                            sliderInput(inputId="covMSATi", label = "new coverage of MSAT (%)", value = 90, min=0, max=100,step=10)
                     ),
                     column(3,
                            sliderInput(inputId="MSATsensC", label = "sensitivity HS RDT (clinical) ", value = 99, min=0, max=100,step=5),
                            sliderInput(inputId="MSATsensA", label = "sensitivity HS RDT (micro detectable, asym)", value = 87, min=0, max=100,step=5),
                            sliderInput(inputId="MSATsensU", label = "sensitivity HS RDT (micro undetectable, asym)", value = 4, min=0, max=100,step=5)
                     )
            ),
            tabPanel(title= strong("Download"),
                     br(),
                     downloadButton("downloadTable", "Download current values of parameters"),
                     downloadButton("downloadplot","Download high resolution figure")),
            tabPanel(title= strong("Restore your parameters"),
                     wellPanel(
                       fileInput(inputId = "file", label ="Your input file:", accept = c(".csv"))
                     )
            ),
            tabPanel(title=strong("User Manual & Help"),
                     br(),
                     tags$ul(tags$li(strong(a(href="https://www.dropbox.com/s/d5q4ldkxtm2az6m/RAI_strategydesigntool_usermanual_03032017.pdf?dl=0", "Download User Manual")))),
                     strong("Contact the developers for any questions and feedback"),
                     tags$ul(
                       tags$li(a(href="http://www.tropmedres.ac/sai-thein-than-tun","Sai Thein Than Tun, "), a(href="mailto:sai@tropmedres.ac","sai@tropmedres.ac")),
                       tags$li(a(href="http://www.tropmedres.ac/researchers/researcher/sompob-saralamba","Sompob Saralamba, "),a(href="mailto:sompob@tropmedres.ac","sompob@tropmedres.ac")),
                       tags$li("Shwe Sin Kyaw"),
                       tags$li("Phetsavanh Chanthavilay"),
                       tags$li("Olivier Celhay, ", a(href="mailto:olivier.celhay@gmail.com","olivier.celhay@gmail.com")),
                       tags$li("Trần Đăng Nguyên"),
                       tags$li("Trần Nguyễn Anh Thư"),
                       tags$li("Daniel M Parker"),
                       tags$li("Professor Maciej F Boni"),
                       tags$li("Professor Arjen M Dondorp"),
                       tags$li(a(href="http://www.tropmedres.ac/researchers/researcher/lisa-white","Professor Lisa White, "), a(href="mailto:lisa@tropmedres.ac","lisa@tropmedres.ac"))
                     ))
  ),
  fluidRow(plotOutput(outputId = "MODEL")),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  hr(),
  fluidRow(h4("          Legend")),
  fluidRow(h4("          Grey solid line: baseline scenario. Blue solid line: elimination strategy scenario.")), 
  fluidRow(h4("          Dark blue solid line: target baseline API. Grey dashed lines: start and end of elimination activities.")),
  fluidRow(h4("          Red dashed line: pre-elimination threshold (API = 1 per 1000 per year)"))
  
)

#non-reactive parameters
# define the number of weeks to run the model
dt<-1/12
startyear<-2007
stopyear<-2023
maxt<-stopyear-startyear
times <- seq(0, maxt, by = dt)
tsteps<-length(times)

ParLabel <- read.table('ParLabel.csv', sep=",", as.is=TRUE)

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


server <- function(input, output, session) {
  scenario_0<-c(EDATon = 0,
                ITNon = 0,
                IRSon = 0,
                MDAon = 0,
                primon = 0,
                MSATon = 0,
                VACon = 0,
                v_same = 0)
  
  scenario_iR<-reactive(c(EDATon = input$EDATon,
                          ITNon = input$ITNon,
                          IRSon = input$IRSon,
                          MDAon = input$MDAon,
                          primon = input$primon,
                          MSATon = input$MSATon,
                          VACon = as.numeric(input$VACon),
                          v_same = input$v_same))
  
  parametersR <- reactive(c(
    bh_max = input$bh_max,                 # bites per human per night
    eta = input$eta,
    covEDAT0 = input$covEDAT0,
    covITN0 = input$covITN0,
    effITN = input$effITN,
    covIRS0 = input$covIRS0,
    effIRS = input$effIRS,
    muC = input$muC,
    muA = input$muA,
    muU = input$muU,
    percfail2018 = input$percfail2018,
    percfail2019 = input$percfail2019,
    percfail2020 = input$percfail2020,
    
    EDATscale = input$EDATscale,
    covEDATi = input$covEDATi,
    ITNscale = input$ITNscale,
    covITNi = input$covITNi,
    IRSscale = input$IRSscale,
    covIRSi = input$covIRSi,
    cmda_1 = input$cmda_1,
    cmda_2 = input$cmda_2,
    cmda_3 = input$cmda_3,
    tm_1 = input$tm_1,          # timing of 1st round [2018 to 2021 - 1 month steps]
    tm_2 = input$tm_2,          # timing of 2nd round [2018+(1/12) to 2021 - 1 month steps]
    tm_3 = input$tm_3,          # timing of 3rd round [2018+(2/12) to 2021 - 1 month steps]
    dm = input$dm,
    lossd = input$lossd,
    
    MSATscale = input$MSATscale,
    covMSATi = input$covMSATi,
    MSATsensC = input$MSATsensC,
    MSATsensA = input$MSATsensA,
    MSATsensU = input$MSATsensU,
    
    #effv_1 = input$effv_1,
    #effv_2 = input$effv_2,
    effv_3 = input$effv_3,
    effv_4 = input$effv_4,
    ka = input$ka,
    delta = input$delta,
    kf = input$kf,
    ks = input$ks
    #vh = input$vh
  ))
  
  #getting back previous parameters
  data <- reactive({read.csv(input$file$datapath)})
  datavalue <- reactive(data()[,2])
  observeEvent(input$file,{
    updateCheckboxInput(session, "EDATon", value = datavalue()[1])
    updateCheckboxInput(session, "ITNon", value = datavalue()[2])
    updateCheckboxInput(session, "IRSon", value = datavalue()[3])
    updateCheckboxInput(session, "MDAon", value = datavalue()[4])
    updateCheckboxInput(session, "primon", value = datavalue()[5])
    updateCheckboxInput(session, "MSATon", value = datavalue()[6])
    updateSliderInput(session, "VACon", value = datavalue()[7])
    
    updateSliderInput(session, "API", value = datavalue()[8])
    
    updateSliderInput(session, "bh_max", value = datavalue()[9])
    updateSliderInput(session, "eta", value = datavalue()[10])
    updateSliderInput(session, "covEDAT0", value = datavalue()[11])
    updateSliderInput(session, "covITN0", value = datavalue()[12])
    updateSliderInput(session, "effITN", value = datavalue()[13])
    updateSliderInput(session, "covIRS0", value = datavalue()[14])
    updateSliderInput(session, "effIRS", value = datavalue()[15])
    updateSliderInput(session, "muC", value = datavalue()[16])
    updateSliderInput(session, "muA", value = datavalue()[17])
    updateSliderInput(session, "muU", value = datavalue()[18])
    updateSliderInput(session, "percfail2018", value = datavalue()[19])
    updateSliderInput(session, "percfail2019", value = datavalue()[20])
    updateSliderInput(session, "percfail2020", value = datavalue()[21])
    updateSliderInput(session, "EDATscale", value = datavalue()[22])
    updateSliderInput(session, "covEDATi", value = datavalue()[23])
    updateSliderInput(session, "ITNscale", value = datavalue()[24])
    updateSliderInput(session, "covITNi", value = datavalue()[25])
    updateSliderInput(session, "IRSscale", value = datavalue()[26])
    updateSliderInput(session, "covIRSi", value = datavalue()[27])
    updateSliderInput(session, "cmda_1", value = datavalue()[28])
    updateSliderInput(session, "cmda_2", value = datavalue()[29])
    updateSliderInput(session, "cmda_3", value = datavalue()[30])
    updateSliderInput(session, "tm_1", value = datavalue()[31])
    updateSliderInput(session, "tm_2", value = datavalue()[32])
    updateSliderInput(session, "tm_3", value = datavalue()[33])
    updateSliderInput(session, "dm", value = datavalue()[34])
    updateSliderInput(session, "lossd", value = datavalue()[35])
    updateSliderInput(session, "MSATscale", value = datavalue()[36])
    updateSliderInput(session, "covMSATi", value = datavalue()[37])
    updateSliderInput(session, "MSATsensC", value = datavalue()[38])
    updateSliderInput(session, "MSATsensA", value = datavalue()[39])
    updateSliderInput(session, "MSATsensU", value = datavalue()[40])
    #updateSliderInput(session, "effv_1", value = datavalue()[41])
    #updateSliderInput(session, "effv_2", value = datavalue()[42])
    updateSliderInput(session, "effv_3", value = datavalue()[43])
    #updateSliderInput(session, "vh", value = datavalue()[44])
    
  })
  
  # initial prevalence
  initprevR <- reactive(0.001*input$API)
  
  GMSout0R <- reactive(runGMS(initprevR(), scenario_0,parametersR()))
  
  GMSoutiR <- reactive(runGMS(initprevR(), scenario_iR(),parametersR()))
  
  plotR <- function()
  {
    GMSout0<-GMSout0R()
    
    GMSouti<-GMSoutiR()
    
    times<-GMSout0[,1]
    clinmonth_det<-cbind(GMSout0[,2],GMSouti[,2])
    clinmonth_tot<-cbind(GMSout0[,3],GMSouti[,3])
    prevalence<-cbind(GMSout0[,4],GMSouti[,4])
    
    runin<-(2016-startyear)/dt
    
    finclin<-max(clinmonth_tot[(runin:length(clinmonth_det[,1])),])
    finprev<-max(prevalence[(runin:length(prevalence[,1])),])
    
    
    # PLOTTING
    par(mfrow=c(1,2), cex=1.5)
    
    maxy<-max(finclin,input$API/12)
    x<-times[(runin:length(clinmonth_det[,1]))]
    y1<-clinmonth_det[runin:length(clinmonth_det[,1]),1]
    y2<-clinmonth_tot[runin:length(clinmonth_tot[,1]),1]
    
    plot(x,y1, type='l',lty=1,col=rgb(0,0,0,alpha=0.1),xlab = "Time",ylab="incidence per 1000 per month",main="Monthly cases per 1000 population",ylim=c(0,maxy),lwd=2)
    lines(x,y2, type='l',lty=1,col=rgb(0,0,0,alpha=0.1),lwd=2)
    
    polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(0,0,0,alpha=0.1),border=NA)
    
    y1<-clinmonth_det[runin:length(clinmonth_det[,1]),2]
    y2<-clinmonth_tot[runin:length(clinmonth_tot[,1]),2]
    lines(x,y1, type='l',lty=1,col=rgb(0,0,1,alpha=0.4),lwd=2)
    lines(x,y2, type='l',lty=1,col=rgb(0,0,1,alpha=0.4),lwd=2)
    
    polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(0,0,1,alpha=0.4),border=NA)
    
    lines(c(2018,2018),c(-maxy,2*maxy),col="dark grey",lty=3,lwd=2)
    
    abline(h=input$API/12,col="dark blue",lty=1,lwd=1)
    abline(h=1/12,col="red",lty=3,lwd=3)
    maxy<-finprev
    plot(times[(runin:length(prevalence[,1]))],prevalence[(runin:length(prevalence[,1])),1], type='l',lty=1,col=rgb(0,0,0,alpha=0.25),xlab = "Time",ylab="% prevalence",main="Predicted true prevalence",ylim=c(0,maxy),lwd=6)
    lines(times[(runin:length(prevalence[,1]))],prevalence[(runin:length(prevalence[,1])),2], type='l',lty=1,col=rgb(0,0,1,alpha=0.6),xlab = "Time",ylab="% prevalence",main="Predicted true prevalence",ylim=c(0,maxy),lwd=6)
    lines(c(2018,2018),c(-maxy,2*maxy),col="dark grey",lty=3,lwd=2)
    
  }
  
  output$MODEL <- renderPlot({
    plotR()
  })
  
  output$downloadplot <- downloadHandler(
    filename = function(){paste('MalMod_',gsub("\\:","",Sys.time()),'.png',sep='')},
    content = function(file) {
      png(filename=file, height= 1600, width=4800, units= "px", res=300) #if(...=="png"){png(file)} else if(...=="pdf"){pdf(file)}
      plotR()
      dev.off()
    })
  
  tableContentR <- reactive({
    tmp <- c(scenario_iR(), input$API, parametersR())
    tmp2 <- cbind(ParLabel[,1], tmp, ParLabel[,2], names(tmp))
    colnames(tmp2) <- c("Name","Value","Unit","VarName")
    tmp2
  })
  
  output$downloadTable <- downloadHandler(
    filename = function(){paste('MalMod_',gsub("\\:","",Sys.time()),'.csv',sep='')},
    content = function(file) {
      write.csv(tableContentR(), file, row.names = FALSE)
    })
}

shinyApp(ui = ui, server = server)