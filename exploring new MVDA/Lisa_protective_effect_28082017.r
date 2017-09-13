# the formula for protective effect of MVDA

library(deSolve)

times <- seq(0, 3, by = (1/52))

#MODEL PARAMETERS
parameters <- c(cm = 365, #12/1, #365,          # rate of mass intervetion deployment
                ka = 52/2,       # 1/(time to full protective effect of vaccine after dose)
                kf = 365*log(2)/45,       # fast rate of loss of vaccine protective effect
                ks = 365*log(2)/634,        # slow rate of loss of vaccine protective effect
                delta = 12/5,    # 1/(time to transition between fast and slow rate for loss of vaccine protection)
                pv0 = 0.9        # maximum starting vaccine protection level
                )

inity <- 1

# MODEL INITIAL CONDITIONS
state <- c(y01=inity, y02 = 0, yf=0, ys=0)

# set up a function to solve the equations
prot<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         cmon<- (t>0.5)*(t<=0.6)*cm
         
         dy01 <- -cmon*y01
         dy02 <- cmon*y01-ka*y02
         dyf <-  ka*y02-(kf+delta)*yf
         dys <-  delta*yf - ks*ys
         
         # return the rate of change
         list(c(dy01, dy02, dyf, dys))
       }
  ) 
  
}


# run the model
out <- ode(y = state, times = times, func = prot, parms = parameters)

pv <- out[,"yf"]+out[,"ys"]

plot(times,pv,type="l",col='blue',ylim=c(0,1))



