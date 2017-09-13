# the formula for protective effect of MVDA

library(deSolve)

times <- seq(0, 3, by = (1/52))

#MODEL PARAMETERS
parameters <- c(m_deploy = 365, #12/1,          # rate of mass intervetion deployment
                ka = 52/2,       # 1/(time to full protective effect of vaccine after dose)
                kf = 5.62,       # fast rate of loss of vaccine protective effect
                ks = 0.4,        # slow rate of loss of vaccine protective effect
                delta = 12/5,    # 1/(time to transition between fast and slow rate for loss of vaccine protection)
                p_max_3 = 1 #0.9        # maximum starting vaccine protection level
                )

inity <- as.vector(parameters['p_max_3'])
intvStart <- .5

# MODEL INITIAL CONDITIONS
state <- c(y01_3=inity, y02_3 = 0, yf_3=0, ys_3=0, Y=0)

# set up a function to solve the equations
prot<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         dY <- 1
         # dx0 <- -m_deploy*x0
         # dx1 <- m_deploy*x0-ld*x1
         dy01_3 <- -m_deploy*y01_3*(Y>intvStart)
         dy02_3 <- m_deploy*y01_3*(Y>intvStart)-ka*y02_3
         dyf_3 <-  ka*y02_3-(kf+delta)*yf_3
         dys_3 <-  delta*yf_3 - ks*ys_3
         
         # return the rate of change
         list(c(dy01_3, dy02_3, dyf_3, dys_3, dY))
       }
  ) 
  
}


# run the model
out <- ode(y = state, times = times, func = prot, parms = parameters)


pv <- out[,"yf_3"]+out[,"ys_3"]

plot(times,pv,type="l",col='blue',ylim=c(0,1))




