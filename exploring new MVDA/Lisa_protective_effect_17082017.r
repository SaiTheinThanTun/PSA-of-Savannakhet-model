# the formula for protective effect of MVDA

library(deSolve)

times <- seq(0, 3, by = (1/52))

#MODEL PARAMETERS
parameters <- c(cm = 12,          # rate of mass intervetion deployment
                ld = 52/2,       # rate of loss of drug effect
                ka = 365,#52/2,       # 1/(time to full protective effect of vaccine after dose)
                kf = 12/1,       # fast rate of loss of vaccine protective effect
                ks = 1/3,        # slow rate of loss of vaccine protective effect
                delta = 12/1,    # 1/(time to transition between fast and slow rate for loss of vaccine protection)
                pv0 = 0.9        # maximum starting vaccine protection level
                )

inity <- 1

# MODEL INITIAL CONDITIONS
state <- c(x0 = 1, x1=0, y01=1, y02 = 0, yf=0, ys=0)

# set up a function to solve the equations
prot<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         dx0 <- -cm*x0
         dx1 <- cm*x0-ld*x1
         dy01 <- -cm*y01
         dy02 <- cm*y01-ka*y02
         dyf <-  ka*y02-(kf+delta)*yf
         dys <-  kf*yf - ks*ys
         
         # return the rate of change
         list(c(dx0, dx1, dy01, dy02, dyf, dys))
       }
  ) 
  
}


# run the model
out <- ode(y = state, times = times, func = prot, parms = parameters)

pd <- out[,"x1"]
pv <- out[,"yf"]+out[,"ys"]

plot(times,pd,type="l",col='blue',ylim=c(0,1))
lines(times,pv,type="l",col='red')



