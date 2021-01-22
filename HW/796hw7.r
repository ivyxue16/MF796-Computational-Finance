library(MASS)

## define a function that performs fft on Heston process
Heston_fft<-function( alpha, n, B, K, params ) {
  S0     = params[6]
  r      = params[7]
  T      = params[8]
  
  N = 2^n
  Eta = B / N
  Lambda_Eta = 2*pi/N
  Lambda = Lambda_Eta/Eta
  
  J = 1:N
  vj = (J-1)*Eta
  m = 1:N
  Beta = log(S0) - Lambda*N / 2
  km = Beta + (m-1)*Lambda
  
  ii <- complex(real=0, imaginary=1)
  #calculate values of characteristic function
  Psi_vj = rep(0, length(J))
  for (zz in 1:N) {
    u <- (vj[[zz]] - (alpha+1.0)*ii)
    numer <- Heston_cf(u, params) 
    denom <- ((alpha + ii* vj[[zz]]) * (alpha + 1.0 + ii*vj[[zz]]))
    
    Psi_vj[[zz]] <- (numer / denom)
  }
  
  #compute fft
  XX = (Eta/2)*Psi_vj*exp(-ii*Beta*vj)*(2-dirac(J-1))
  ZZ = fft(XX)
  
  #calculate option prices
  Multiplier = exp(-alpha*km)/pi
  ZZ2 <- Multiplier*ZZ
  Km <- exp(km)
  
  #discard strikes that are 0 or infinity to avoid errors in interpolation
  inds <- (!is.infinite(Km) & !is.nan(Km) & (Km > 1e-16) & (Km < 1e16))
  px_interp <- approx(Km[inds], Re(ZZ2[inds]), method = "linear", xout=K)
  
  fft_price = Re(exp(-r*T)*px_interp$y)
  return(fft_price)
}



## define a dirac delta function
dirac<-function(n) {
  y <- rep(0, length(n))
  y[n==0] = 1
  return(y)
}


## define a function that computes the characteristic function for variance gamma
Heston_cf<-function(u, params) {
  sigma <- params[1]
  eta0 <- params[2]
  kappa <- params[3]
  rho     <- params[4]
  theta     <- params[5]
  S0     <- params[6]
  r      = params[7]
  q      = params[8]
  T      = params[9]
  
  ii <- complex(real=0, imaginary=1)
  
  l = sqrt(sigma^2*(u^2+ii*u)+(kappa-ii*rho*sigma*u)^2)
  
  w = exp(ii*u*log(S0)+ii*u*(r-q)*T + 
            kappa*theta*T*(kappa-ii*rho*sigma*u)/sigma^2)/(cosh(l*T/2)+(kappa-ii*rho*sigma*u)/l*sinh(l*T/2))^(2*kappa*theta/sigma^2)
  
  y = w*exp(-(u^2+ii*u)*eta0/(l/tanh(l*T/2)+kappa-ii*rho*sigma*u))
  
  return(y)
}


## set parameters
kappa <- 3.51
theta <- 0.052
sigma <- 1.17
rho   <- -0.77
eta0  <- 0.034
S0 <- 282
r  <- 0.015
q  <- 0.0177
expT <- 1
K1 <- 285
K2 <- 315
N <- 20000
M <- 252

alpha <- 1
n <- 15
B <- 10000.0


## compute heston fft price
call_heston <- Heston_fft(alpha, n, B, K1, c(sigma,eta0,kappa,rho,theta,S0,r,q,expT))


# compute Heston simulation
dt = expT/M # time step size
Sigma <- dt*matrix(c(1,rho,rho,1),2,2) # the covariance matrix of the mutivariate normal distribution
s_price <- matrix(0,N,M) # stock prices matrix 
vol_matrix <- matrix(0,N,M) # vol matrix
s_price[,1] <- S0  # initial price
vol_matrix[,1] <- eta0  # initial volatility

# Monte Carlo Simulation
for (i in 2:M) {
  # mvrnorm(n=1,mu,Sigma) # library(MASS)
  z <- mvrnorm(n = N, rep(0, 2), Sigma)
  s_price[,i] <- s_price[,i-1]+(r-q)*s_price[,i-1]*dt+sqrt(pmax(vol_matrix[,i-1],0))*s_price[,i-1]*z[,1]
  vol_matrix[,i] <- vol_matrix[,i-1]+kappa*(theta-vol_matrix[,i-1])*dt+sigma*sqrt(pmax(vol_matrix[,i-1],0))*z[,2]
}

## price of options
UpandOutCall_payoff <- pmax(s_price[,M]-K1,0)*exp(-r*expT)
UpandOutCall_payoff[apply(s_price,1,max)>K2] <- 0
UpandOutCall_sim <- sum(UpandOutCall_payoff)/N

Euro_payoff <- pmax(s_price[,M]-K1,0)*exp(-r*expT)
Eurocall_sim <- exp(-r*expT)*sum(Euro_payoff)/N

# find the control variate c
c <- -cov(Euro_payoff,UpandOutCall_payoff)/var(Euro_payoff)
# compute the price of control varaite method
control_variate <- UpandOutCall_sim + c*(Eurocall_sim-call_heston)


# output
UpandOutCall_sim
Eurocall_sim

# draw simulation histogram to decide which is the best N
x <- matrix(0,200)

for(j in 1:200){
  
  # compute Heston simulation
  dt = expT/M # time step size
  Sigma <- dt*matrix(c(1,rho,rho,1),2,2) # the covariance matrix of the mutivariate normal distribution
  s_price <- matrix(0,N,M) # stock prices matrix 
  vol_matrix <- matrix(0,N,M) # vol matrix
  s_price[,1] <- S0  # initial price
  vol_matrix[,1] <- eta0  # initial volatility

# Monte Carlo Simulation
  for (i in 2:M) {
    # mvrnorm(n=1,mu,Sigma) # library(MASS)
    z <- mvrnorm(n = N, rep(0, 2), Sigma)
    s_price[,i] <- s_price[,i-1]+(r-q)*s_price[,i-1]*dt+sqrt(pmax(vol_matrix[,i-1],0))*s_price[,i-1]*z[,1]
    vol_matrix[,i] <- vol_matrix[,i-1]+kappa*(theta-vol_matrix[,i-1])*dt+sigma*sqrt(pmax(vol_matrix[,i-1],0))*z[,2]
  }
  UpandOutCall_payoff <- pmax(s_price[,M]-K1,0)*exp(-r*expT)
  UpandOutCall_payoff[apply(s_price,1,max)>K2] <- 0
  UpandOutCall_sim <- sum(UpandOutCall_payoff)/N
  Euro_payoff <- pmax(s_price[,M]-K1,0)*exp(-r*expT)
  Eurocall_sim <- exp(-r*expT)*sum(Euro_payoff)/N
  
  
  # find the control variate c
  c <- -cov(Euro_payoff,UpandOutCall_payoff)/var(Euro_payoff)
  # compute the price of control varaite method
  x[j] <- UpandOutCall_sim + c*(Eurocall_sim-call_heston)
}
  

hist(x,main="N=20000")



# compute the difference of price with FFT
Euro_err <- (Euro_call_sim/call_heston)-1

