# Run ODE model, written by Swati Patel, updated by Abigail D'Ovidio Long
# parameters from Emiljanowicz
# stage-dependent death based on Mermer 2021
library(deSolve)


# Overall SWD parameters
a_dev = c(1/1.4, 1/6, 1/5.8, 491.1/160) #development and reproduction
mu = c(0.627, 0.289, 0.162, 0.088) #death rate
s=c(0.9,0.0,0.0,0.0) #selection coefficients
h = c(0.5,0.5,0.5,0.5) # heterozygosity


# Parameters for each pesticide

#Malathion
mal_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.149,0.276,0.020,0.100), #survival of susceptible (Malathion)
  s=s, #selection coefficients
  h =h # heterozygosity
)

#Zeta-cypermethrin
zeta_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.478,0.405,0.008,0.003), #survival of susceptible (Zeta-cypermethrin)
  s=s, 
  h =h 
)

#Spinosad
spino_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.172,0.241,0.003,0.052), #survival of susceptible (Spinosad)
  s=s, 
  h =h 
)

#Spinetoram
spinet_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.156,0.189,0.008,0.05), #survival of susceptible (Spinetoram)
  s=s, 
  h =h 
)

#Methomyl
met_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.157,0.230,0.0001,0.01), #survival of susceptible (Methomyl)
  s=s, 
  h =h 
)

#Cyantraniliprole
cya_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.272,0.308,0.050,0.05), #survival of susceptible (Cyantraniliprole)
  s=s, 
  h =h 
)

#Fenpropathrin
fen_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.370,0.335,0.050,0.05), #survival of susceptible (Fenpropathrin)
  s=s, 
  h =h 
)

#Phosmet
pho_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.089,0.085,0.005,0.010), #survival of susceptible (Phosmet)
  s=s, 
  h =h 
)

#Cyclaniliprole
cyc_parms<-list(
  a_dev = a_dev, 
  mu = mu, 
  delta = c(0.327,0.282,0.300,0.300), #survival of susceptible (Cyclaniliprole)
  s=s, 
  h =h 
)

# time parameters and initiation
tau <- 20  # time between pulses
maxt <- 300
times<-seq(0, maxt, by=1)
treat.times <- seq(tau, maxt, by=tau)

init.cond <- c(x1=80, x2=80, x3=80, x4=80, y1=0.1, y2=0.1, y3=0.1, y4=0.1)



#ode solver for each pesticide 
mal_ans <- ode(y=init.cond, times=times, func=odefn, parms=mal_parms, events=list(func = treatfun, time = treat.times))
zeta_ans <- ode(y=init.cond, times=times, func=odefn, parms=zeta_parms, events=list(func = treatfun, time = treat.times))
spino_ans <- ode(y=init.cond, times=times, func=odefn, parms=spino_parms, events=list(func = treatfun, time = treat.times))
spinet_ans <- ode(y=init.cond, times=times, func=odefn, parms=spinet_parms, events=list(func = treatfun, time = treat.times))
met_ans <- ode(y=init.cond, times=times, func=odefn, parms=met_parms, events=list(func = treatfun, time = treat.times))
cya_ans <- ode(y=init.cond, times=times, func=odefn, parms=cya_parms, events=list(func = treatfun, time = treat.times))
fen_ans <- ode(y=init.cond, times=times, func=odefn, parms=fen_parms, events=list(func = treatfun, time = treat.times))
pho_ans <- ode(y=init.cond, times=times, func=odefn, parms=pho_parms, events=list(func = treatfun, time = treat.times))
cyc_ans <- ode(y=init.cond, times=times, func=odefn, parms=cyc_parms, events=list(func = treatfun, time = treat.times))


# adult values v. time
ans_egg <- data.frame(mal_ans[,1], cya_ans[,6], cyc_ans[,6], fen_ans[,6], mal_ans[,6],
                      met_ans[,6], pho_ans[,6], spinet_ans[,6], spino_ans[,6],
                      zeta_ans[,6])
ans_larva <- data.frame(mal_ans[,1], cya_ans[,7], cyc_ans[,7], fen_ans[,7], mal_ans[,7],
                        met_ans[,7], pho_ans[,7], spinet_ans[,7], spino_ans[,7],
                        zeta_ans[,7])
ans_pupa <- data.frame(mal_ans[,1], cya_ans[,8], cyc_ans[,8], fen_ans[,8], mal_ans[,8],
                       met_ans[,8], pho_ans[,8], spinet_ans[,8], spino_ans[,8],
                       zeta_ans[,8])
ans_adult <- data.frame(mal_ans[,1], cya_ans[,9], cyc_ans[,9], fen_ans[,9], mal_ans[,9],
                        met_ans[,9], pho_ans[,9], spinet_ans[,9], spino_ans[,9],
                        zeta_ans[,9])


# plotting
cols <- c("red", "blue4", "black", "darkgreen", 'sienna', 'darkorange2',
          "hotpink1", "green2", "turquoise3")

layout(
  matrix(c(1, 2, 5,
           3, 4, 5), nrow = 2, byrow = TRUE),
  widths = c(1, 1, 0.35)
)

par(mar = c(4.5, 4, 2.5, 1))

# Egg selection
matplot(
  ans_egg[,1], ans_egg[,2:10],
  type = "l", lty = 1, col = cols,
  xlab = "Time", ylab = "Resistance in Adults",
  main = expression(s[1] == 0.9)
)

# Larva selection
matplot(
  ans_larva[,1], ans_larva[,2:10],
  type = "l", lty = 1, col = cols,
  xlab = "Time", ylab = "Resistance in Adults",
  main = expression(s[2] == 0.9)
)

# Pupa Selection
matplot(
  ans_pupa[,1], ans_pupa[,2:10],
  type = "l", lty = 1, col = cols,
  xlab = "Time", ylab = "Resistance in Adults",
  main = expression(s[3] == 0.9)
)

# Adult Selection
matplot(
  ans_adult[,1], ans_adult[,2:10],
  type = "l", lty = 1, col = cols,
  xlab = "Time", ylab = "Resistance in Adults",
  main = expression(s[4] == 0.9)
)

# Legend
par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  "center",
  legend = c(
    "Cyanatanilipole",
    "Cyclaniliprole",
    "Fenpropathrine",
    "Malathion",
    "Methomyl",
    "Phosmet",
    "Spinetoram",
    "Spinosad",
    "Zeta-cypermethrin"
  ),
  col = cols,
  lty = 1,
  bty = "n",
)
