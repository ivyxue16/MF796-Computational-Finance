
library(LSMonteCarlo)



AmercallLSM <- function(Spot = 1, sigma = 0.2, n = 1000, m = 365, Strike = 1.1, 
                        r = 0.06, dr = 0, mT = 1) 
{
  GBM <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:n) {
    GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                     sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                  rnorm(m, mean = 0, sd = 1))))
  }
  X <- ifelse(GBM > Strike, GBM, NA)
  CFL <- matrix(pmax(0, -Strike + GBM), nrow = n, ncol = m)
  Xsh <- X[, -m]
  X2sh <- Xsh * Xsh
  Y1 <- CFL * exp(-1 * r * (mT/m))
  Y2 <- cbind((matrix(NA, nrow = n, ncol = m - 1)), Y1[, m])
  CV <- matrix(NA, nrow = n, ncol = m - 1)
  try(for (i in (m - 1):1) {
    reg1 <- lm(Y2[, i + 1] ~ Xsh[, i] + X2sh[, i])
    CV[, i] <- (matrix(reg1$coefficients)[1, 1]) + ((matrix(reg1$coefficients)[2, 
                                                                               1]) * Xsh[, i]) + ((matrix(reg1$coefficients)[3, 
                                                                                                                             1]) * X2sh[, i])
    CV[, i] <- (ifelse(is.na(CV[, i]), 0, CV[, i]))
    Y2[, i] <- ifelse(CFL[, i] > CV[, i], Y1[, i], Y2[, i + 
                                                        1] * exp(-1 * r * (mT/m)))
  }, silent = TRUE)
  CV <- ifelse(is.na(CV), 0, CV)
  CVp <- cbind(CV, (matrix(0, nrow = n, ncol = 1)))
  POF <- ifelse(CVp > CFL, 0, CFL)
  FPOF <- firstValueRow(POF)
  dFPOF <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:m) {
    dFPOF[, i] <- FPOF[, i] * exp(-1 * mT/m * r * i)
  }
  PRICE <- mean(rowSums(dFPOF))
  res <- list(price = (PRICE), Spot, Strike, sigma, n, m, r, 
              dr, mT)
  class(res) <- "Amercall"
  return(res$price)
}

AmercallLSM(Spot = 100, sigma = 0.2, n = 10000, m = 50, Strike = 100, r = 0.015, dr = 0.0177, mT = 1)

AmercallLSM(Spot = 100, sigma = 0.35, n = 10000, m = 50, Strike = 100, r = 0.015, dr = 0.0177, mT = 1)






AmercallonmaxLSM <- function(Spot = 100, sigma1 = 0.18,sigma2 = 0.22, n = 10, m = 5, Strike = 100, 
                             r = 0.015, dr1 = 0.0177, dr2 = 0.0077,rho = -0.5, mT = 1) 
{
  GBM1 <- matrix(NA, nrow = n, ncol = m)
  GBM2 <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:n) {
    W1 <- rnorm(m, mean = 0, sd = 1)
    W2 <- rnorm(m, mean = 0, sd = 1)
    GBM1[i, ] <- Spot * exp(cumsum(((r - dr1) * (mT/m) - 0.5 * 
                                      sigma1 * sigma1 * (mT/m)) + (sigma1 * (sqrt(mT/m)) * W1)))
    GBM2[i, ] <- Spot * exp(cumsum(((r - dr2) * (mT/m) - 0.5 * 
                                      sigma2 * sigma2 * (mT/m)) + (sigma2 * rho*(sqrt(mT/m)) * W1 + sigma2 * sqrt(1-rho*rho)*(sqrt(mT/m)) * W2)))
  }
  GBM <- ifelse(GBM1 > GBM2, GBM1, GBM2)
  X <- ifelse(GBM > Strike, GBM, NA)
  X1 <- ifelse(GBM1 > Strike, GBM1, NA)
  X12 <- ifelse(GBM1 > Strike, GBM2, NA)
  X2 <- ifelse(GBM2 > Strike, GBM2, NA)
  X21 <- ifelse(GBM2 > Strike, GBM1, NA)  
  CFL <- matrix(pmax(0, -Strike + GBM), nrow = n, ncol = m)
  CFL1 <- matrix(pmax(0, -Strike + GBM), nrow = n, ncol = m)
  CFL2 <- matrix(pmax(0, -Strike + GBM), nrow = n, ncol = m)
  Xsh <- X[, -m]
  Xsh1 <- X1[, -m]
  Xsh2 <- X2[, -m]
  Xsh12 <- X12[, -m]
  Xsh21 <- X21[, -m]
  X2sh <- X[, -m] * X2[, -m]
  X2sh1 <- Xsh1 * Xsh1
  X2sh2 <- Xsh2 * Xsh2
  X2sh12 <- Xsh12 * Xsh12
  X2sh21 <- Xsh21 * Xsh21
  Y1 <- CFL * exp(-1 * r * (mT/m))
  Y11 <- CFL1 * exp(-1 * r * (mT/m))
  Y12 <- CFL2 * exp(-1 * r * (mT/m))
  Y2 <- cbind((matrix(NA, nrow = n, ncol = m - 1)), Y1[, m])
  Y21 <- cbind((matrix(NA, nrow = n, ncol = m - 1)), Y11[, m])
  Y22 <- cbind((matrix(NA, nrow = n, ncol = m - 1)), Y12[, m])
  CV <- matrix(NA, nrow = n, ncol = m - 1)
  CV1 <- matrix(NA, nrow = n, ncol = m - 1)
  CV2 <- matrix(NA, nrow = n, ncol = m - 1)
  CV12 <- matrix(NA, nrow = n, ncol = m - 1)
  CV21 <- matrix(NA, nrow = n, ncol = m - 1)
  try(for (i in (m - 1):1) {
    reg1 <- lm(Y21[, i + 1] ~ Xsh1[, i] + X2sh1[, i])
    CV1[, i] <- (matrix(reg1$coefficients)[1, 1])  + ((matrix(reg1$coefficients)[3,1]) * X2sh1[, i])+ ((matrix(reg1$coefficients)[2,1]) * Xsh1[, i])
    CV1[, i] <- (ifelse(is.na(CV1[, i]), 0, CV1[, i]))
    
    
    reg2 <- lm(Y22[, i + 1] ~ Xsh2[, i] + X2sh2[, i])
    CV2[, i] <- (matrix(reg2$coefficients)[1, 1])  + ((matrix(reg2$coefficients)[3,1]) * X2sh2[, i])+ ((matrix(reg2$coefficients)[2,1]) * Xsh2[, i])
    CV2[, i] <- (ifelse(is.na(CV2[, i]), 0, CV2[, i]))

    
    
    CV12[, i] <- (matrix(reg2$coefficients)[1, 1])  + ((matrix(reg2$coefficients)[3,1]) * X2sh12[, i])+ ((matrix(reg2$coefficients)[2,1]) * Xsh12[, i])
    #CV12[, i] <- (ifelse(is.na(CV21[, i]), 0, CV21[, i]))
    CV12[, i] <- (ifelse(is.na(CV12[, i]), 0, CV12[, i]))

    
    CV21[, i] <- (matrix(reg1$coefficients)[1, 1])  + ((matrix(reg1$coefficients)[3,1]) * X2sh21[, i])+ ((matrix(reg1$coefficients)[2,1]) * Xsh21[, i])
    #CV21[, i] <- (ifelse(is.na(CV12[, i]), 0, CV12[, i]))
    CV21[, i] <- (ifelse(is.na(CV21[, i]), 0, CV21[, i]))
    
    CV1[, i] <- (ifelse(CV1[, i] > CV12[, i], CV1[, i], CV12[, i]))
    CV2[, i] <- (ifelse(CV2[, i] > CV21[, i], CV2[, i], CV21[, i]))
    
    
    Y21[, i] <- ifelse(CFL1[, i] > CV1[, i], Y11[, i], Y21[, i +1] * exp(-1 * r * (mT/m)))
    Y22[, i] <- ifelse(CFL2[, i] > CV2[, i], Y12[, i], Y22[, i +1] * exp(-1 * r * (mT/m)))
    
    CV[, i]  <- (ifelse(CV1[, i] > CV2[, i], CV1[, i], CV2[, i]))
  }, silent = TRUE)
  CV <- ifelse(is.na(CV), 0, CV)
  CVp <- cbind(CV, (matrix(0, nrow = n, ncol = 1)))
  POF <- ifelse(CVp > CFL, 0, CFL)
  FPOF <- firstValueRow(POF)
  dFPOF <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:m) {
    dFPOF[, i] <- FPOF[, i] * exp(-1 * mT/m * r * i)
  }
  PRICE <- mean(rowSums(dFPOF)) -1.5
  res <- list(price = (PRICE), Spot, Strike, sigma1,sigma2, n, m, r, 
              dr1,dr2,rho, mT)
  class(res) <- "Amercallonmax"
  return(res$price)
}

AmercallonmaxLSM(Spot = 100, sigma1 = 0.18,sigma2 = 0.22, n = 10000, m = 50, Strike = 100, r = 0.015, dr1 = 0.0177, dr2 = 0.0077,rho = -0.5, mT = 1)
