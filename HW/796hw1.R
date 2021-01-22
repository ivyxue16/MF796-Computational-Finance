# MF796 Assignment 1
# Problem 5
# (b)
# initial parameters
s0 <- 100
r <- 0
beta <- 1
sigma <- 0.4
nd <- 250
simt <- 10000
K=100

# simulaton function
stksim <- function(s0=100, r=0,beta=1,sigma=0.4,nd=250,simt=10000){
  dt = 1/nd
  sim = data.frame(matrix(rep(0,(nd+1)*simt),nd+1,simt))
  sim[1,] = s0
  for(i in 2:(nd+1)){
    sim[i,] = sim[i-1,] + r*sim[i-1,] + sigma*(sim[i-1,]**beta)*rnorm(simt,0,dt**0.5)
  }
  sp = as.matrix(sim[251,])
  return(sp)
}

sp <- stksim()
sp <- t(sp)

class(sp)
dim(sp)
# European call option price from simulation

poff <- matrix(rep(0,simt))
dim(poff)
for(i in 1:simt){
  poff[i] <- max(sp[i]-K,0)
}
sim_price <- exp(-r*nd) * mean(poff)

paste("(b) European call option price from simulation is : ",sim_price)


# (c) European Call option pirce from BS Formula
t = 1
d1 = (log(s0/K) + (r + 0.5*(sigma**2))*t) / (sigma*sqrt(t))
d2 <- d1 - sigma*sqrt(t)

call_price <- s0 * pnorm(d1) - K*exp(-r*t)*pnorm(d2)
put_price <- K*exp(-r*t)*pnorm(-d2) - s0 * pnorm(-d1)

paste("(c) European Call option pirce from BS Formula is :", call_price)
paste("European Put option pirce from BS Formula is :", put_price)


# (d) delta
delta <- pnorm(d1)
paste("(d) Delta is:",delta)


# (e) delta-neutral hedge
shares <- -delta
paste("(e) shares of stock is" ,shares)

#(f)
payoff <- matrix(rep(0,simt))
for(i in 1:simt){
  payoff[i] <- max(sp[i]-100,0) + delta * (100-sp[i])
}

a <- mean(payoff)
paste("(f) The payoff of delta-neutral portfolio is ",a)
