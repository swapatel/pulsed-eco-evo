# Run ODE model, written by Swati Patel, updated by Abigail D'Ovidio Long
# parameters from Emiljanowicz
# stage-dependent death based on Mermer 2021
library(deSolve)


# Overall SWD parameters
a_dev = c(1/1.4, 1/6, 1/5.8, 491.1/160) #development and reproduction
mu = c(0.627, 0.289, 0.162, 0.088) #death rate
s=c(0.0,0.0,0.9,0.0) #selection coefficients
h = c(0.5,0.5,0.5,0.5) # heterozygosity


#Methomyl Parameters
met_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.157,0.230,0.0001,0.01), #survival of susceptible (Methomyl)
  s=s, #selection coefficients
  h =h # heterozygosity
)


#time parameters and initialization
tau <- 20  # time between pulses
maxt <- 300
times<-seq(0, maxt, by=1)
treat.times <- seq(tau, maxt, by=tau)

init.cond <- c(x1=80, x2=80, x3=80, x4=80, y1=0.1, y2=0.1, y3=0.1, y4=0.1)



#ode solver 
met_ans <- ode(y=init.cond, times=times, func=odefn, parms=met_parms, events=list(func = treatfun, time = treat.times))


# plotting
cols <- c('sienna')

matplot(met_ans[,1], met_ans[,9], type="l",lty = 1, col = cols, xlab="Time", ylab="Resistance Frequency (Adults)", main = "Selection in Pupa for Methomyl")
