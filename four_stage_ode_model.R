# 4 stage eco-evo pulsed ode, written by Swati Patel and 
# updated by Abigail D'Ovidio Long and Abigail Adjei
# used for figures 3 and 4

odefn <- function(t, n, parms) {
  with(as.list(c(parms, n)), {
    dx1 = -(mu[1] + a_dev[1])*x1 + a_dev[4]*x4
    dx2 = a_dev[1]*x1 - (mu[2] + a_dev[2])*x2
    dx3 = a_dev[2]*x2 - (mu[3] + a_dev[3])*x3
    dx4 = a_dev[3]*x3 - (mu[4])*x4
    
    dy1 = a_dev[4]*x4/x1*(y4-y1)
    dy2 = a_dev[1]*x1/x2*(y1-y2)
    dy3 = a_dev[2]*x2/x3*(y2-y3)
    dy4 = a_dev[3]*x3/x4*(y3-y4)
    
    return(list(c(dx1, dx2, dx3, dx4, dy1, dy2, dy3, dy4))) 
  })
}

# treatment function, each divided by (1-s) for resistant proportion
compute_marginal <- function(y_, d_, s_, h_){
  return(y_*d_/(1-s_) + (1-y_)*(1-s_*h_)*d_/(1-s_))
  # marginal fitness of resistant allele
}

#mean fitness
compute_meanfit <- function(y_, d_, s_, h_){
  d22 <- d_/(1-s_)
  d12 <- d_*(1-s_*h_)/(1-s_)
  d11 <- d_
  return(y_^2*d22 + 2*y_*(1-y_)*d12 + (1-y_)^2*d11)
}

# apply treatment to each stage
treatfun <- function(t, z, parms){
  with(as.list(parms, z), {
    #affects each stage density
    x <- z[1:4]
    y <- z[5:8]
    y1 <- y[1]
    y2 <- y[2]
    y3 <- y[3]
    y4 <- y[4]
    
    wbar1 <- compute_meanfit(y1, delta[1], s[1], h[1]) 
    wbar2 <- compute_meanfit(y2, delta[2], s[2], h[2])
    wbar3 <- compute_meanfit(y3, delta[3], s[3], h[3])
    wbar4 <- compute_meanfit(y4, delta[4], s[4], h[4])
    
    r1 <- compute_marginal(y1, delta[1], s[1], h[1])
    r2 <- compute_marginal(y2, delta[2], s[2], h[2])
    r3 <- compute_marginal(y3, delta[3], s[3], h[3])
    r4 <- compute_marginal(y4, delta[4], s[4], h[4])
    
    z[1] <- wbar1*z[1]
    z[2] <- wbar2*z[2]
    z[3] <- wbar3*z[3]
    z[4] <- wbar4*z[4]
    
    z[5] <- r1/wbar1*z[5]
    z[6] <- r2/wbar2*z[6]
    z[7] <- r3/wbar3*z[7]
    z[8] <- r4/wbar4*z[8]
    
    return(z)
  })
}
