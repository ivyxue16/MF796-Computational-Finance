#define a function that performs fft on Heston process
Heston_fft<-function( alpha, n, B, K, params ) {
  r      = params[7]
  T      = params[8]
  S0     = params[6]

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


#define a dirac delta function
dirac<-function(n) {
  y <- rep(0, length(n))
  y[n==0] = 1
  return(y)
}


#define a function that computes the characteristic function for variance gamma
Heston_cf<-function(u, params) {
  sigma <- params[1]
  eta0 <- params[2]
  kappa <- params[3]
  rho     <- params[4]
  theta     <- params[5]
  S0     <- params[6]
  r      = params[7]
  T      = params[8]

  ii <- complex(real=0, imaginary=1)

  l = sqrt(sigma^2*(u^2+ii*u)+(kappa-ii*rho*sigma*u)^2)

  w = exp(ii*u*log(S0)+ii*u*(r-0)*T +
            kappa*theta*T*(kappa-ii*rho*sigma*u)/sigma^2)/(cosh(l*T/2)+(kappa-ii*rho*sigma*u)/l*sinh(l*T/2))^(2*kappa*theta/sigma^2)

  y = w*exp(-(u^2+ii*u)*eta0/(l/tanh(l*T/2)+kappa-ii*rho*sigma*u))

  return(y)
}

#######Q1########
sigma = 0.2
eta0 = 0.08
kappa = 0.7
rho = -0.4
theta = 0.1
S0 = 250
r = 0.02
expT = 0.5

params <- c(sigma,eta0,kappa,rho,theta,S0,r,expT)

K = seq(80, 120, by=2.5)
t1 = Sys.time()
alpha <- 20
n <- 14
B <- 250.0
Heston_fft(alpha, n, B, K, params )[9]
t2 = Sys.time()
# print(call_px)
# put_px = vg_fft(-alpha, n, B, K, params )
# print(put_px)

K = seq(80, 120, by=2.5)
alpha <- 1
temp = c()
for (n in c(9,10,11,12,13)){
  for (B in c(100,150,200,250)){
    t1 = Sys.time()
    temp2 = (Heston_fft(alpha, n, B, K, params )[9]-11.41)^2
    t2 = Sys.time()
    temp = c(temp,1/(temp2*as.numeric(t2-t1)))
  }
}
plot(temp,type = "l",ylab = "Efficiency")

alpha <- 1
n <- 14
B <- 250
Heston_fft(alpha, n, B, 105, params )

alpha <- 1
temp = c()
for (n in c(9,10,11,12,13)){
  for (B in c(100,150,200,250)){
    t1 = Sys.time()
    temp2 = (Heston_fft(alpha, n, B, 105, params )-9.13)^2
    t2 = Sys.time()
    temp = c(temp,1/(temp2*as.numeric(t2-t1)))
  }
}
plot(temp,type = "l",ylab = "Efficiency")



#######Q2########
#Define the B-S model and the root searching function to get the implied volatility from "market price"

blsPrice = function(K, S, vol, r, T){
  d1 = (log(S/K)+(r+vol^2/2)*T)/(vol*sqrt(T))
  d2 = d1 - vol*sqrt(T)
  return(pnorm(d1)*S - K*exp(-r*T)*pnorm(d2))
}
blsPriceError<-function(vol, mkt_px, K, S, r, T)
{
  err <- blsPrice(K, S, vol, r, T)-mkt_px
  return(err)
}

blsImpVol<-function(K, S, px, r, T, tol=1e-8, maxiter=1000){
  vol = uniroot(blsPriceError, lower = -10, upper = 10,
                mkt_px = px, K = K, S = S, T = T,
                r = r,tol = tol, maxiter = maxiter)
  return(vol$root)
}


###(i)###
sigma = 0.4
eta0 = 0.09
kappa = 0.5
rho = 0.25
theta = 0.12
S0 = 150
r = 0.025
expT = 3/12

params <- c(sigma,eta0,kappa,rho,theta,S0,r,expT)

K = seq(80, 120, by=2.5)
alpha <- 1
n <- 14
B <- 250.0
price = Heston_fft(alpha, n, B, K, params )

Imp_vol = rep(0,length(K))

for (i in 1:length(K)){
  Imp_vol[i] = blsImpVol(K[i],S0,price[i],r,expT)

}

plot(K,Imp_vol,type = "l")



###(ii)###
sigma = 0.6
eta0 = 0.02
kappa = 1
rho = -0.25
theta = 0.5
S0 = 150
r = 0.025
K = 150
expT = seq(1:24)/12

alpha <- 1
n <- 14
B <- 250.0
price = rep(0,24)


for (i in 1:24){
  params <- c(sigma,eta0,kappa,rho,theta,S0,r,expT[i])
  price[i] = Heston_fft(alpha, n, B, K, params )
}


Imp_vol = rep(0,24)
test = rep(0,24)
for (i in 1:24){
  Imp_vol[i] = blsImpVol(K,S0,price[i],r,expT[i])
  test[i] = blsPrice(K,S0,theta,r,expT[i])
}

plot(expT,Imp_vol,type = "l")
