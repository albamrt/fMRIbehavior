# install.packages(c("circular", "dplyr", "ggforce", "knitr", "lme4", "lmerTest", "Publish", "tidyr", "tidyverse", "Hmisc"))
# remotes::install_github('eddjberry/mixturer')
library(readr)
library(Hmisc)
require(dplyr)
require(tibble)
require(tidyverse)
require(tidyr)
require(circular)
require(knitr)
require(Publish)
require(ggforce)
require(nlme)
require(Hmisc)
require(reticulate)
require(data.table)
require(readxl)
require(kableExtra)
require(lme4)
library(mixturer)



summaryHist <- function(data, variable){
  ggplot(data, aes(x = data[,variable])) +
    geom_histogram(aes(y =..density..), bins = 30, color = "gray",
                   position = 'identity', alpha = 0.4) +
    labs(x = label(data[,variable]))
}

fitMixtStim <- function(data){
  mod <- JV10_fit(X = data$R_rad, Tg = data$S_rad, NT = data$prevstim_rad, return.ll = TRUE)
  mod <- setNames(as.data.frame(mod),c("K", "Pt", "Pnt", "Pu", "LL"))
  mod
}

fitMixtResp <- function(data){
  mod <- JV10_fit(X = data$R_rad, Tg = data$S_rad, NT = data$prevresp_rad, return.ll = TRUE)
  mod <- setNames(as.data.frame(mod),c("K", "Pt", "Pnt", "Pu", "LL"))
  mod
}

fitMixtProbe <- function(data){
  mod <- JV10_fit(X = data$R_rad, Tg = data$S_rad, NT = data$prevprob_rad, return.ll = TRUE)
  mod <- setNames(as.data.frame(mod),c("K", "Pt", "Pnt", "Pu", "LL"))
  mod
}

YgVal <- function(cdat, ndat, g) {
  N <- length(cdat) ; ndatcsum <- cumsum(ndat)
  delhat <- 0 ; tbar <- 0
  for (k in 1:g) {
    sample <- circular(0)
    if (k==1) {low <- 0} else
      if (k > 1) {low <- ndatcsum[k-1]}
    for (j in 1:ndat[k]) { sample[j] <- cdat[j+low] }
    tm1 <- trigonometric.moment(sample, p=1)
    tm2 <- trigonometric.moment(sample, p=2)
    Rbar1 <- tm1$rho; Rbar2 <- tm2$rho ; tbar[k] <- tm1$mu
    delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
  }
  dhatmax <- max(delhat) ; dhatmin <- min(delhat)
  if (dhatmax/dhatmin <= 4) {
    CP <- 0 ; SP <- 0 ; dhat0 <- 0
    for (k in 1:g) {
      CP <- CP+ndat[k]*cos(tbar[k]) ; SP <- SP+ndat[k]*sin(tbar[k])
      dhat0 <- dhat0+ndat[k]*delhat[k] }
    dhat0 <- dhat0/N
    RP <- sqrt(CP*CP+SP*SP) ; Yg <- 2*(N-RP)/dhat0
    return(Yg) } else
      if (dhatmax/dhatmin > 4) {
        CM <- 0 ; SM <- 0 ; Yg <- 0
        for (k in 1:g) {
          CM <- CM+(ndat[k]*cos(tbar[k])/delhat[k])
          SM <- SM+(ndat[k]*sin(tbar[k])/delhat[k])
          Yg <- Yg+(ndat[k]/delhat[k]) }
        RM <- sqrt(CM*CM+SM*SM) ; Yg <- 2*(Yg-RM)
        return(Yg) }
}


serial_bias <- function(prevcurr, error, window, step){
  d <- data.frame(prevcurr = prevcurr, error = error)
  d <- d[complete.cases(d),]
  xxx = seq(-pi/2, pi/2, step)
  m_err <- c()
  std_err <- c()
  for(t in xxx){
    idx = ((d$prevcurr>=t-window/2) & (d$prevcurr<t+window/2))
    if(t-window/2 < -pi/2){
      idx = ((d$prevcurr>=t-window/2) & (d$prevcurr<t+window/2) | (d$prevcurr>pi/2-(window/2-(pi/2-abs(t)))))
    }
    if(t+window/2 > pi/2){
      idx = ((d$prevcurr>=t-window/2) & (d$prevcurr<t+window/2) | (d$prevcurr<=pi/2+(window/2-(pi/2-abs(t)))))
    }
    m_err <- c(m_err, mean.circular(circular(d$error[idx], 
                                             units = 'radians', 
                                             template = 'geographics')))
    std_err <- c(std_err, (sd.circular(circular(d$error[idx],
                                                units = 'radians', template = 'geographics'))
                           /sqrt(sum(idx))))
  }
  return(data.frame(x = xxx, m = m_err, sd = std_err, 
                    yl = m_err - std_err, yh = m_err + std_err))}


folded_serial_bias <- function(prevcurr, error, window, step){
  d <- data.frame(prevcurr = abs(prevcurr), error = error*sign(prevcurr))
  d <- d[complete.cases(d),]
  xxx = seq(-pi/2, pi/2, step)
  m_err <- c()
  std_err <- c()
  for(t in xxx){
    idx = ((d$prevcurr>=t-window/2) & (d$prevcurr<t+window/2))
    m_err <- c(m_err, mean.circular(circular(d$error[idx], 
                                             units = 'radians', 
                                             template = 'geographics')))
    std_err <- c(std_err, (sd.circular(circular(d$error[idx],
                                                units = 'radians', template = 'geographics'))
                           /sqrt(sum(idx))))
  }
  return(data.frame(x = xxx, m = m_err, sd = std_err, 
                    yl = m_err - std_err, yh = m_err + std_err))}


Fourier <- function (form, q = 1, minq = 0, maxq = 0, crit = "gcv", data = NULL) 
{
  cat("Reminder:  first explanatory variable is used for fourier expansion", 
      "\n")
  mat <- model.frame(form, data = data)
  y <- mat[, 1]
  z <- mat[, 2]
  n = length(y)
  nx = ncol(mat) - 2
  if (nx > 0) {
    xmat <- as.matrix(mat[, 3:ncol(mat)])
  }
  xnames <- colnames(mat)[3:ncol(mat)]
  minz = min(z)
  maxz = max(z)
  z <- 2 * pi * (z - minz)/(maxz - minz)
  square <- z^2
  searchq = maxq > minq
  if (searchq == FALSE) {
    maxq = q
  }
  sinvar <- array(0, dim = c(n, maxq))
  cosvar <- array(0, dim = c(n, maxq))
  for (j in seq(1, maxq)) {
    sinvar[, j] <- sin(j * z)
    cosvar[, j] <- cos(j * z)
  }
  qstar = maxq
  newform <- y ~ z + square
  if (nx > 0) {
    newform <- update(newform, ~. + xmat)
  }
  if (searchq == TRUE) {
    fit <- lm(newform)
    k = length(fit$coef)
    sig2 <- mean(residuals(fit)^2)
    if (crit == "gcv") {
      critq = n * (n * sig2)/((n - k)^2)
    }
    if (crit == "sc") {
      critq = log(sig2) + log(n) * k/n
    }
    if (crit == "aic") {
      critq = log(sig2) + 2 * k/n
    }
    qstar = 0
    mincrit = critq
    cat("Information Criterion, Linear:    ", mincrit, "\n")
    newform <- update(newform, ~. + sinvar[, 1:j] + cosvar[, 
                                                           1:j])
    for (j in seq(minq, maxq)) {
      fit <- lm(newform)
      k = length(fit$coef)
      sig2 <- mean(residuals(fit)^2)
      if (crit == "gcv") {
        critq = n * (n * sig2)/((n - k)^2)
      }
      if (crit == "sc") {
        critq = log(sig2) + log(n) * k/n
      }
      if (crit == "aic") {
        critq = log(sig2) + 2 * k/n
      }
      cat("Information Criterion, q =", j, ":", critq, 
          "\n")
      if (critq < mincrit) {
        qstar = j
        mincrit = critq
      }
    }
    cat("Information Criterion minimizing q = ", qstar, "\n")
    if (qstar == maxq) {
      cat("Warning:  best q = maximum allowed; may want to try higher value for maxq", 
          "\n")
    }
  }
  maxq = ifelse(searchq == TRUE, qstar, q)
  newform <- y ~ z + square
  if (qstar > 0) {
    sinvar <- as.matrix(sinvar[, 1:maxq])
    cosvar <- as.matrix(cosvar[, 1:maxq])
    newform <- update(newform, ~. + sinvar[, 1:maxq] + cosvar[, 
                                                              1:maxq])
  }
  if (nx > 0) {
    newform <- update(newform, ~. + xmat)
  }
  fit <- lm(newform)
  k = length(fit$coef)
  sig2 <- mean(residuals(fit)^2)
  yhat <- fitted(fit)
  rss = sum(residuals(fit)^2)
  sig2 = rss/n
  aic = log(sig2) + 2 * k/n
  sc = log(sig2) + log(n) * k/n
  gcv = n * rss/((n - k)^2)
  fourierhat <- yhat
  nx1 = 3 + 2 * maxq + 1
  nx2 = nx1 + nx - 1
  if (nx > 0) {
    bhat <- fit$coefficients[nx1:nx2]
    xbhat <- as.array(as.matrix(xmat) %*% as.matrix(bhat))
    fourierhat <- yhat - xbhat + mean(xbhat)
    names(fit$coefficients)[nx1:nx2] <- xnames
  }
  names(fit$coefficients)[1:3] <- c("Intercept", "z", "square")
  sinnames <- paste("sin(", c(1:maxq), sep = "")
  sinnames <- paste(sinnames, "z)", sep = "")
  cosnames <- paste("cos(", c(1:maxq), sep = "")
  cosnames <- paste(cosnames, "z)", sep = "")
  if(!anyNA(names(fit$coefficients)[4:(4 + maxq - 1)])){
    names(fit$coefficients)[4:(4 + maxq - 1)] <- sinnames
    names(fit$coefficients)[(4 + maxq):(4 + 2 * maxq - 1)] <- cosnames
  }
  fourierhat <- fourierhat - mean(fourierhat) + mean(y)
  print(summary(fit))
  out <- list(yhat, rss, sig2, aic, sc, gcv, fit$coef, fourierhat, 
              qstar)
  names(out) <- c("yhat", "rss", "sig2", "aic", "sc", "gcv", 
                  "coef", "fourierhat", "q")
  return(yhat)
}



# DoG: normalized first derivative of a Gaussian with fixed location hyperparameter mu=0
dnorm_deriv1 <- function(x, mean = 0, sd = 1) {
  return(-(x/(sd^2))*dnorm(x, mean = mean, sd = sd))
} 


dnorm_deriv3 <- function(x, mean = 0, sd = 1) {
  return((((-3*x*sd^2)-x^3)/(sd^6))*dnorm(x, mean, sd = 1))
} 

find_nearest <- function(array, value){
  array = as.matrix(array)
  idx = which.min(abs(array-value))
  return(idx)}

normgauss <- function(xxx,sigma){
  gauss = (1/(sigma*sqrt(2*pi))*exp(-(xxx-0)^2 / (2*sigma^2)))
  return(gauss/max(gauss))
}

normgrad <- function(xxx){
  require("numDeriv")
  return(c(diff(xxx))/max(c(diff(xxx))))
}

dog1 <- function(sigma, x){
  xxx 	= seq(-pi/2, pi/2, by= .0001) 
  dog_1st = normgrad(normgauss(xxx,sigma))
  return(unlist(lapply(x, function(x)
    if(is.na(x)){return(NA)}else{dog_1st[find_nearest(xxx,x)]})))
}

dog3 <- function(sigma,x){
  xxx 	= seq(-pi/2, pi/2, by= .0001) 
  dog_3rd = normgrad(normgrad(normgrad(normgauss(xxx,sigma))))
  return(unlist(lapply(x, function(x)
    if(is.na(x)){return(NA)}else{dog_3rd[find_nearest(xxx,x)]})))
}


crossValSigma <- function(sigma, data, derivative){
  if(derivative == 1){
    data$DoG <- dnorm_deriv1(data$diffstim, sd = sigma)}
  #data$DoG <- dog1(sigma = sigma, x = data$diffstim)}
  else{if(derivative == 3){
    data$DoG <- dnorm_deriv1(data$diffstim, sd = sigma)}}
  #data$DoG <- dog3(sigma = sigma, x = data$diffstim)}}
  train <- data %>%
    group_by(subject, time) %>%
    sample_frac(.67)
  test <- anti_join(data, train, by = c('subject', 'time', 'trial')) # test dataframe with observations not in 'train.'
  mod <- lm(fourier_error ~ time*DoG, data = train)
  pred <- predict(mod, test)
  MSE <- mean((na.omit(pred - test$fourier_error))^2)
  return(MSE)
}