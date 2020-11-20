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

fourier_fit <- function(x = NULL, n = NULL, up = 10L, plot = TRUE, add = FALSE, main = NULL, ...){
  N <- length(x)
  #The direct transformation
  #The first frequency is DC, the rest are duplicated
  dff = fft(x)
  #The time
  t = seq(from = 1, to = length(x))
  #Upsampled time
  nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
  #New spectrum
  ndff = array(data = 0, dim = c(length(nt), 1L))
  ndff[1] = dff[1] #Always, it's the DC component
  if(n != 0){
    ndff[2:(n+1)] = dff[2:(n+1)] #The positive frequencies always come first
    #The negative ones are trickier
    ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
  }
  #The inverses
  indff = fft(ndff/(N+1), inverse = TRUE)
  idff = fft(dff/(N+1), inverse = TRUE)
  if(plot){
    if(!add){
      plot(x = t, y = x, pch = 16L, xlab = "Bin", ylab = "Measurement",
           main = ifelse(is.null(main), paste(n, "harmonics"), main))
      lines(y = idff, x = t, col = adjustcolor(1L, alpha = 0.5))
    }
    lines(y = indff, x = nt, ...)
  }
  ret = data.frame(time = nt, y = Mod(indff))
  return(ret)
}