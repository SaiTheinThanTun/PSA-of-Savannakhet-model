#MDA and VACCINE effect explored through a compartmental model
library(deSolve)

times <- Y <- seq(0,(2023-2007), by=1/12)
startyear <- 2018

#tsteps<-length(times)


parameters <- c(
  MDAon = 1,
  tm_1 = 2018+(9/12),
  tm_2 = 2018+(10/12),
  tm_3 = 2018+(11/12),
  tm_4 = 2018+(9/12)+1,
  dm = 6/12, #dur to complete each round
  cm_1 = .80,
  cm_2 = .95,
  cm_3 = .95,
  cmda_1 =.5,
  cmda_2 =.5,
  cmda_3 =.5,
  lossd = 365/15,
  VACon = 1,
  effv_1 = .75,
  effv_2 = .8,
  effv_3 = .92,
  vh = 90/365
)


state <- c(Y=0, X0 = 10000, X1 = 0, Xe = 0, V0 = 10000, V1 = 0, V2 = 0, Ve = 0)

MVDAsub <- function(t, state, parameters){
  with(as.list(c(state, parameters)),
       {
         
         # define parameters
         m_1 <- MDAon*(Y>(tm_1-startyear))*(Y<=(tm_1+dm-startyear))*(-log((1-cm_1))/dm)
         m_2 <- MDAon*(Y>(tm_2-startyear))*(Y<=(tm_2+dm-startyear))*(-log((1-cm_2))/dm)
         m_3 <- MDAon*(Y>(tm_3-startyear))*(Y<=(tm_3+dm-startyear))*(-log((1-cm_3))/dm) 
         m_4 <- MDAon*(Y>(tm_4-startyear))*(Y<=(tm_4+dm-startyear))*(-log((1-cm_3))/dm) #// coverage of m4 is the same as m3
         m_5 <- 0
         
         # rate of change for MDA
         dY <- 1
         dX0 <- -X0*m_1*cmda_1
         dX1 <- X0*m_1*cmda_1-X1*lossd
         dXe <- X1*lossd
         
         # rate of change for vaccination
         dV0 <- 
         
         # return the rate of change
         list(c(dY, dX0, dX1, dXe, dV0, dV1, dV2, dVe))
       }
  ) 
  
}
out <- ode(y = state, times = times, func = MVDAsub, parms = parameters)

percDrug <- out[,4]/rowSums(out[,3:5])
plot(percDrug[1:30], type='l')
