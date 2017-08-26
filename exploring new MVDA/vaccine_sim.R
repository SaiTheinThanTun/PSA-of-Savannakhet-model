#testing collective protection 2 weeks after vaccination
#at the end of the scale up
library(manipulate)

manipulate({

n.event <- dur.scale.week #10
pop.size <- pop.size #1000 population size
no.p <- pop.size/n.event #100 #no. of people vaccinated each week
sim.week <- sim.week # 40 #no. of simulation weeks
eff <- eachMaxProtect #.90 #proportion of people immune out of all people vaccinated per week

delay.week <- 2
people.immune <- rep(0,times=n.event) #storage for no. of people immune each week

#timestep is per week
hl <- 90/7 #90 days


#eff*exp(-log(2)*(t/hl)^.8)
perc.people.protected <- NA

for(i in 1:sim.week){ #i for timestepping
  for(j in 1:n.event){ #j for storing values across vector
    if(i>=j+delay.week) people.immune[j] <- no.p*eff*exp(-log(2)*(max((i-j),0)/hl)^.8) #eff*exp(-log(2)*((i-j)/hl)^.8)
    else people.immune[j] <- 0
    #people.immune[j] <- eff*exp(-log(2)*(max((i-j),0)/hl)^.8)
  }
  perc.people.protected[i] <- sum(people.immune)/pop.size
}
plot(perc.people.protected, main=paste("percent with protection in ",pop.size," persons\n vaccinated over ",n.event," weeks, eachMaxProtect: ",eff,"\noverallMaxProtect: ",round(max(perc.people.protected),2), sep=""), type='l', , xlab='Time (weeks)')
},
dur.scale.week=slider(1,42, initial=10),
pop.size=slider(1000,5000),
sim.week=slider(40,600, initial=100),
eachMaxProtect=slider(.6,.9, initial=.9)
)