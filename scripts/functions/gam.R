# GAM (Helper) Functions


##################################################################
##                       EWS GAM Function                       ##
##################################################################

# This function calculates EWS using GAM. It returns a dataframe with:
# autoregressive coefficient, as well as its upper and lower limit (ll,
# ul), the intercept upper and lower limit, the mu upper and lower limit.
# There are also inbuilt checks for the GAM: test of the BIC,
# significance, edf, whether there is a rise in the model (by comparing
# the 35th fitted value, to the 35th fitted value from the end).
# 
# Scaling denotes whether you want to standardize the variable or not. For
# the circumplex variables (PA/NA hi/lo), scaling has been done prior to
# combining (averaging).
# 
# Mu denotes the mu/attractor of the proces (= intercept/(1-AR)).
# 
# Mean denotes the fitted values using a GAM with time as a predictor.
# This might represent an 'easier' version of mu (i.e. version that does
# not depend on the AR).
# 
# The sections that are commented out were written by M. J. Schreuder for
# a different paper, where the rise in AR was determined using the
# derivative of the two weeks before a transition. We will not apply that
# here.

set.seed(123) # for reproducible results
GAM.ews <- function(dw) {
  set.seed(123) 
  out <- data.frame("time" = dw$time, "check1.BIC.says.tv" = NA, "check2.sig.p.int" = NA, "check2.sig.p.AR" = NA, "check3.edf.tv" = NA, "check4a.rise" = NA,
                    "BIC.ar" = NA, "BIC.tvar" = NA, "pval.int.tvar" = NA, "pval.AR.tvar" = NA, "edf" = NA,
                    "intercept" = NA, 'll.int' = NA, 'ul.int' = NA, "AR" = NA, "ll"=NA, "ul"=NA,"mu" = NA, 'll.mu' = NA, 'ul.mu' = NA, 
                    "mean" = NA, "ll.mean" = NA, "ul.mean" = NA,
                    'gam.check' = NA)
  # Note that more columns are created in the function! 
  newd <- data.frame(time = dw$time, VarL = 1) 
  # fit a gam model 
  tvar.mod <- gam(Var ~ s(time,k=10,bs="tp") + s(time, by = VarL, k=10,bs= "tp"), data = dw) 
  ar.mod <- gam(Var ~ VarL, data = dw) # add normal AR
  out$pval.int.tvar <- summary(tvar.mod)[[8]][1]
  out$check2.sig.p.int <- summary(tvar.mod)[[8]][1]< .05
  out$pval.AR.tvar <-  summary(tvar.mod)[[8]][2]
  out$check2.sig.p.AR <-  summary(tvar.mod)[[8]][2] < .05
  # fit a second gam model to check the effect of the mean
  mod.mean <- gam(Var ~ s(time, bs = "tp", k = 10), data = dw) #for gam.check
  BICmin <- which.min(c(BIC(tvar.mod),BIC(ar.mod))) 
  out$BIC.ar <- BIC(ar.mod)
  out$BIC.tvar <- BIC(tvar.mod)
  if(BICmin == 2){
    best.model <- ar.mod
    out$check1.BIC.says.tv <- FALSE
  }else{
    best.model <- tvar.mod
    out$check1.BIC.says.tv <- TRUE
  }
  if(length(tvar.mod$fitted.values) > 71){
    out$check4a.rise <- (tvar.mod$fitted.values[35] < best.model$fitted.values[length(tvar.mod$fitted.values)-35])
  } else{
    out$check4a.rise <- "Cannot test"
  }
  Xp <- predict(tvar.mod, newd, type="lpmatrix", seWithMean = TRUE) 
  Xp.mean <- predict(mod.mean, newd, type="lpmatrix", seWithMean = TRUE)
  res.gam <- data.frame(intercept = Xp[,1:10] %*% coef(tvar.mod)[1:10], AR = Xp[,11:20] %*% coef(tvar.mod)[11:20], mean = Xp.mean %*% coef(mod.mean))
  out$intercept <- res.gam$intercept # this is how the intercept changes over time
  out$AR <- res.gam$AR # this is how the autocorrelation changes over time
  out$mean <- res.gam$mean # this is how the variable changes over time
  nboot <- 10000 # note that nboot is also used for calculating the confidence intervals around the derivative
  modr <- mvrnorm(nboot, coef(tvar.mod), tvar.mod$Vp+diag(20)*10^(-30)) # draw nboot samples from a normal distribution with mu = our estimates and sigma = the variance of our estimates. 
  modr.mean <- mvrnorm(nboot, coef(mod.mean), mod.mean$Vp+diag(10)*10^(-30)) 
  # The confidence intervals around the AR estimate. 
  # Simpson calls this an simultaneous confidence interval: we're estimating the CI of the entire smooth rather than a single point.
  phi.ci <- matrix(NA, nrow(dw), nboot) 
  beta_0.sm <- matrix(NA, nrow(dw), nboot) 
  mu.sm <- matrix(NA, nrow(dw), nboot)
  mean.sm <- matrix(NA, nrow(dw), nboot)
  for (m in 1 : nboot) {
    phi.ci[,m] <- Xp[,11:20] %*% modr[m,11:20] # calculate how autocorrelation changes over time for this set of estimates. This will be used for derivatives later!
    beta_0.sm[,m] <- Xp[,1:10] %*% modr[m,1:10] 
    mu.sm[,m] <- beta_0.sm[,m]/(1-phi.ci[,m]) # note: mu = intercept/(1 - AR), see Bringmann et al., 2017
    mean.sm[,m] <-  Xp.mean %*% modr.mean[m,]
  }
  out[,c('ll.int', 'ul.int')] <- cbind(apply(beta_0.sm,1,quantile,c(.025)), apply(beta_0.sm,1,quantile,c(.975))) # for each row (time point), get the boundaries of 95%. 
  out[,c('ll', 'ul')] <- cbind(apply(phi.ci,1,quantile,c(.025)), apply(phi.ci,1,quantile,c(.975))) # for each row (time point), get the boundaries of 95%. 
  out[,c('ll.mu', 'mu', 'ul.mu')] <- t(apply(mu.sm,1,quantile,c(.025,.5,.975))) # the mu, based on simulations, and its 95% CI.
  out[,c('ll.mean', 'ul.mean')] <- apply(mean.sm,1,quantile,c(.025,.975)) 
  # Calculate the derivative of the second smooth to check whether EWS are increasing or not.
  # Recall: if the derivative is not 0, this means that the line is either in- or decreasing significantly.
  # Using central finite differences (= the average of forward and backward finite differences)
  # eps <- 1e-7
  # newd2 <- newd
  # newd2$time <- newd$time + eps #shift time a little bit forward
  # Xfd <- predict(mod,newd2,type="lpmatrix",seWithMean = T) #use this to calculate first derivative
  # newd2$time <- newd$time - eps #shift time a little bit backward
  # Xfd2 <- predict(mod,newd2,type="lpmatrix",seWithMean = T) 
  #       
  # # Calculate centered difference approximation
  # X0 <- (Xfd - Xfd2)/(2*eps) #delta y/delta x. Note that X0 still contains the values of the linear predictor (i.e. still has to be multiplied by coef(mod)).
  #       
  # # Calculate confidence interval around derivative
  # # NB: if something goes wrong, it might have to do with the matrix multiplication (to transpose or not to transpose??)
  # fd <- X0[,11:20] %*% coef(mod)[11:20] #first derivative of the autocorrelation
  # sim.fd <- X0[,11:20] %*% t(modr)[11:20,] #simulated first derivative
  # sim.fd.int <- X0[,1:10] %*% t(modr)[1:10,] #simulated first derivative of the intercept
  # fd.CI <- apply(sim.fd, 1, quantile, probs = c(0.025, 0.975)) #upper and lower confidence interval (i.e. for each time point, the limits wherein 95% of the simulations lie)
  # 
  # # Do the same for the mean: 
  # newd2$time <- newd$time + eps #shift time a little bit forward
  # Xfd <- predict(mod.mean,newd2,type="lpmatrix",seWithMean = T) #use this to calculate first derivative
  # newd2$time <- newd$time - eps #shift time a little bit backward
  # Xfd2 <- predict(mod.mean,newd2,type="lpmatrix",seWithMean = T) 
  #       
  # # Calculate centered difference approximation
  # X0 <- (Xfd - Xfd2)/(2*eps) 
  # fd.mean <- X0 %*% coef(mod.mean)
  # sim.fd <- X0 %*% t(modr.mean)
  # fd.CI.mean <- apply(sim.fd, 1, quantile, probs = c(0.025, 0.975)) 
  # Calculate the first derivative of mu. Note that to do this we use the quotientrule.
  # Recall that mu = intercept/(1 - AR).
  # Hence, according to the quotientrule, the first derivative of mu equals ((f.acc * g) - (g.acc * f)) / (g^2), where f = intercept, f.acc = first derivative of intercept, g = (1-AR), and g.acc = first derivative of (1-AR).
  # f.acc <- X0[,1:10] %*% coef(mod)[1:10] #first derivative of the intercept
  # f <- out$intercept
  # g <- 1-out$AR
  # g.acc <- -1 * fd
  # fd.mu <- ((f.acc * g) - (g.acc * f)) / (g^2)
  # 
  # # Now simulate a confidence interval around the first derivative. 
  # # In other words, calculate nboot derivatives by resampling the coefficients from the original model!
  # fd.CI.mu <- matrix(NA, nrow = nboot, ncol = nrow(dw))
  # for (i in 1:nboot) {
  #  g.acc <- -1 * sim.fd[,i]
  #  fd.CI.mu[i,] <- ((sim.fd.int[,i] * phi.ci[,i]) - (g.acc * beta_0.sm[,i])) / (phi.ci[,i] ^2) }
  # 
  # fd.CI.mu <- apply(sim.fd, 1, quantile, probs = c(0.025, 0.975))
  # 
  # fd <- data.frame(fd = fd, lwr = fd.CI[1,], upr = fd.CI[2,], fd.mu = fd.mu, lwr.mu = fd.CI.mu[1,], upr.mu = fd.CI.mu[2,], fd.mean = fd.mean, lwr.mean = fd.CI.mean[1,], upr.mean = fd.CI.mean[2,])
  # fd$sig <- ifelse(fd$lwr > 0 | fd$upr < 0, 1, 0) #if the confidence interval does not include 0, sig equals 1
  # fd$incr.AR <- ifelse((fd$fd > 0 & fd$lwr > 0), 'yes', 'no')
  # fd$decr.AR <- ifelse((fd$fd < 0 & fd$upr < 0), 'yes', 'no')
  # fd$sig.mu <- ifelse(fd$lwr.mu > 0 | fd$upr.mu < 0, 1, 0)
  # fd$incr.mu <- ifelse((fd$fd.mu > 0 & fd$lwr.mu > 0), 'yes', 'no')
  # fd$decr.mu <- ifelse((fd$fd.mu < 0 & fd$upr.mu < 0), 'yes', 'no')
  # fd$sig.mean <- ifelse(fd$lwr.mean > 0 | fd$upr.mean < 0, 1, 0)
  # fd$incr.mean <- ifelse((fd$fd.mean > 0 & fd$lwr.mean > 0), 'yes', 'no')
  # fd$decr.mean <- ifelse((fd$fd.mean < 0 & fd$upr.mean < 0), 'yes', 'no')
  # out[,names(fd)] <- fd
  # Check the original GAM model: enough basis functions? 
  c <- k.check(tvar.mod)
  p.val <- c[2,4] # If p-value is low, the number of basis functions might have been too low...
  edf <- c[2,2] #...especially if the edf is close to 10!
  out$check3.edf.tv <- edf > 2 
  out$edf <- edf
  if (p.val < 0.1 & edf > 9) {
    out$gam.check <- 'problem'
    sprintf('there were too few basis functions')
  } else {
    out$gam.check <- 'no.problem'
  }
  #out$tvar.gam.mod <- list(tvar.mod)
  # # Denote segments, needed for plotting (quite complicated loop but reduces code lines)
  # for (i in c('incr', 'decr')) {
  #   for(j in c('.AR', '.mu', '.mean')) {
  #     out[,paste(paste('segment', i, sep='_'), j, sep = '')] <- 1
  #     for (p in 2:nrow(out)) {
  #     if (out[p, paste(i, j, sep = '')] != out[p-1, paste(i, j, sep = '')]) { # if the increase of this time point is not similar to the last one
  #       out[p:nrow(out), paste(paste('segment', i, sep='_'), j, sep = '')] <- out[p-1, paste(paste('segment', i, sep='_'), j, sep = '')] + 1 # a new segment has started
  #     }}}} 
  return(out) # out contains 24 columns
}

set.seed(123) # for reproducible results
GAM.sum <- function(dw, k) {
  set.seed(123) 
  out <- data.frame("time" = dw$time, "edf" = NA, "pval.s.time" = NA, "check.sig.s.time" = NA, 
                    "BIC.y" = NA, "BIC.lag" = NA, "BIC.ar" = NA, "best.model" = NA,
                    "mean" = NA, "ll.mean" = NA, "ul.mean" = NA,
                    'gam.check' = NA)
  # Note that more columns are created in the function! 
  
  # fit a gam model 
  gam_y <- gam(Var ~ s(time, k=k, bs="tp"), data = dw) 
  gam_lag <- gam(Var ~ VarL, data = dw)
  gam_ar <- gamm(Var ~ s(time), data = dw, correlation = corAR1(form = ~ time))
  
  # p-value smooth term
  out$pval.s.time <- summary(gam_y)[[8]][1]
  out$check.sig.s.time <- summary(gam_y)[[8]][1]< .05
  
  # compare models and save best model
  BICmin <- which.min(c(BIC(gam_y),BIC(gam_lag), BIC(gam_ar$lme))) 
  out$BIC.y <- BIC(gam_y)
  out$BIC.lag <- BIC(gam_lag)
  out$BIC.ar <- BIC(gam_ar$lme)
  out$best.model <- c("gam_y", "gam_lag", "gam_ar")[[BICmin]]
  best.model <- c("gam_y", "gam_lag", "gam_ar")[[BICmin]] %>% get # save best model
  
  # predicted values and base functions
  newd <- data.frame(time = dw$time) 
  Xp <- predict(gam_y, newd, type="lpmatrix", seWithMean = TRUE) 
  
  Xlag <- predict(gam_ar$gam, newd, type="lpmatrix", seWithMean = TRUE) 
  meanLag <- cbind(newd, pred = Xlag[,1:k] %*% coef(gam_ar$gam)[1:k])
  
  
  res.gam <- data.frame(mean = Xp[,1:k] %*% coef(gam_y)[1:k])
  out$mean <- res.gam$mean # this is how the variable changes over time
  nboot <- 10000 # note that nboot is also used for calculating the confidence intervals around the derivative
  modr.mean <- mvrnorm(nboot, coef(gam_y), gam_y$Vp+diag(k)*k^(-30)) # draw nboot samples from a normal distribution with mu = our estimates and sigma = the variance of our estimates. 
  # The confidence intervals around the AR estimate. 
  # Simpson calls this an simultaneous confidence interval: we're estimating the CI of the entire smooth rather than a single point.
  mean.sm <- matrix(NA, nrow(dw), nboot)
  for (m in 1 : nboot) {
    mean.sm[,m] <-  Xp %*% modr.mean[m,]
  }
  out[,c('ll.mean', 'ul.mean')] <- apply(mean.sm,1,quantile,c(.025,.975)) 
  
  # Check the original GAM model: enough basis functions? 
  c <- k.check(gam_y)
  p.val <- c[1,4] # If p-value is low, the number of basis functions might have been too low...
  edf <- c[1,2] #...especially if the edf is close to 10!
  out$edf <- edf
  
  if (p.val < 0.1 & edf > 9) {
    out$gam.check <- 'problem'
    message('there were too few basis functions')
  } else {
    out$gam.check <- 'no.problem'
  }
  
  return(out) # out contains 24 columns
}

# EWS GAM helper function 
ewsdf_function <- function(df, range) {
  ews <- NA
  
  if (any(df$incr.AR[df$time %in% range] == 'yes', na.rm = T)) {
    ews <- c(ews, 'AR') }
  if ('mu' %in% names(df)) {
    if (any(df$incr.mu[df$time %in% range] == 'yes', na.rm = T)) {
      ews <- c(ews, 'mu') } 
  }
  if ('mean' %in% names(df)) {
    if (any(df$incr.mean[df$time %in% range] == 'yes', na.rm = T)) {
      ews <- c(ews, 'mean') }
  }
  if ('sd' %in% names(df)) {
    if (any(df$incr.sd[df$time %in% range] == 'yes', na.rm = T)) {
      ews <- c(ews, 'sd') }
  }
  
  ews <- ews[!is.na(ews)]
  
  if (length(ews) == 0) {
    k <- 'none of the indicators' 
  } else {
    k <- paste(ews[1])
    ll <- length(ews)
    if (ll > 1) {
      for (f in 2:ll) {
        k <- paste(k, paste(ews[f]), sep = ', ')
      }}}
  
  return (k)
  
}
tautab <- function(df) {
  out <- data.frame(tau_AR = NA, P_tau_AR = NA, corr.P_AR = NA, tau_mean = NA, P_tau_mean = NA, corr.P_mean = NA, 
                    tau_sd = NA, P_tau_sd = NA, corr.P_sd = NA, incr.AR = NA, incr.mean = NA, incr.sd = NA)
  
  for (estimate in c('AR', 'sd', 'mean')) {
    ds <- subset(df, select = c('time', estimate))
    names(ds) <- c('time', 'estimate')
    out[,paste("tau", estimate, sep = '_')] <- Kendall(ds$estimate, ds$time)$tau 
    out[,paste("P_tau", estimate, sep = '_')]  <- Kendall(ds$estimate, ds$time)$sl # P value of tau (non-corrected)
    out[,paste("corr.P", estimate, sep = '_')] <- mmkh(ds$estimate[!is.na(ds$estimate)])[2] # corrected P value
    out[,paste("incr", estimate, sep = '.')] <- ifelse((out[,paste("corr.P", estimate, sep = '_')] < 0.05) & (out[,paste("tau", estimate, sep = '_')] > 0), 'yes', 'no') # if tau is significant (with Hamed-Rao correction) and greater than 0, the indicator (AR, sd, or mean) is increasing.
  }
  return(out)
}

# EWS GAM Plotting Functions
helpfun <- function(tp, mw = FALSE, winsize = NULL)  {
  #  if (tp[3] != 'no transition(s)') { # if there was one transition
  l <- as.numeric(tp[1:2])
  if (mw) { l[1] <- l[1] + winsize-1 } # if the moving window technique was used, the first AR estimate is the length of the window (i.e. 70).
  no.tp <- 1 #}
  return(list(l, no.tp))
}
transition_line <- function(plot, tp, l) {
  if (tp[3] == 'no transition(s)') { # if there was no transition
    # plot <- plot + geom_vline(xintercept = l[2], col = ifelse(tp[3] == 'depression', 'blue', ifelse(tp[3] == 'mania', 'red', 'black'))) }
    plot <- plot + geom_vline(xintercept = l[2], col = 'white', alpha = 1 ) }
  if (tp[3] != 'no transition(s)') { # if there was a transition
    plot <- plot + geom_vline(xintercept = l[2], 
                              col = ifelse(tp[3] == '1 week(s)', 'purple', 
                                           ifelse(tp[3] == '2 week(s)', 'seagreen3', 
                                                  ifelse(tp[3] == '3 week(s)', 'red2', 'orange')))) + 
      annotate("rect", ymin = -Inf, ymax = Inf, alpha = 0.2, xmin = l[2],
               xmax = ifelse(tp[3] == '1 week(s)', l[2]+(5*7),
                             ifelse(tp[3] == '2 week(s)', l[2]+(5*14),
                                    ifelse(tp[3] == '3 week(s)', l[2]+(5*21), l[2]+(5*28)))),
               fill = ifelse(tp[3] == '1 week(s)', 'purple',
                             ifelse(tp[3] == '2 week(s)', 'seagreen3',
                                    ifelse(tp[3] == '3 week(s)', 'red2', 'orange')))) + theme(legend.position = "none")
  }
  return(plot)
}
tp_fun <- function(plot, tp, data, mw = FALSE) {
  names(data) <- c('time', "estimate", "ll", "ul")
  l <- helpfun(tp, mw = mw)[[1]]
  no.tp <- helpfun(tp, mw = mw)[[2]]
  plot <- transition_line(plot, tp, l)
  for (n in no.tp) { # add the line and its confidence intervals
    nw.data <- data[c(which(data$time == l[n]) : which(data$time == l[n+1])),]
    plot <- plot + geom_line(data=nw.data, aes(x = as.numeric(time), y = as.numeric(estimate)), inherit.aes=FALSE)
    plot <- plot + geom_ribbon(data=nw.data, aes(x = as.numeric(time), ymin = as.numeric(ll), ymax = as.numeric(ul)), inherit.aes=FALSE, fill = 'black', alpha = 0.2)
  }
  return(plot)
}
sign.change <- function(plot, data) {
  names(data) <- c('time', 'estimate', 'incr', 'decr', 'segment_incr', 'segment_decr')
  if(any(data$incr[!is.na(data$incr)] == 'yes')) {
    for (k in unique(data$segment_incr[data$incr == 'yes'])) {
      plot <- plot + geom_line(data = data[data$segment_incr == k & data$incr == 'yes',], aes(time, estimate), col = 'red') }
  }
  if(any(data$decr[!is.na(data$decr)] == 'yes')) {
    for (k in unique(data$segment_decr[data$decr == 'yes'])) {
      plot <- plot + geom_line(data = data[data$segment_decr == k & data$decr == 'yes',], aes(time, estimate), col = 'blue') }
  }
  return(plot)
}
tau_fun <- function(plot, tp, data, mw = FALSE, winsize = NULL) {
  names(data) <- c('time', 'estimate')
  l <- helpfun(tp, mw = mw, winsize = winsize)[[1]]
  no.tp <- helpfun(tp, mw = mw, winsize = winsize)[[2]]
  for (n in no.tp) { # for each pretransition period
    # calculate tau based on three weeks before the transition (= 105 observations prior to transition)
    df <- data[c(which(data$time == c(l[n+1] - 105)) : which(data$time == l[n+1])),]
    tau <- Kendall(df$estimate, df$time)$tau # calculate tau
    p_tau <- Kendall(df$estimate, df$time)$sl # P value of tau (non-corrected)
    p.corrected <- mmkh(df$estimate)[2] # corrected P value
    y <- c(max(df$estimate, na.rm=T), min(df$estimate, na.rm=T))
    x <- c(min(df$time[!is.na(df$estimate)]), max(df$time[!is.na(df$estimate)]))
    x.coord <- c(min(df$time[!is.na(df$estimate)]), max(df$time[!is.na(df$estimate)]))
    half <- round(((max(df$time[!is.na(df$estimate)]) -  min(df$time[!is.na(df$estimate)]))/2) + min(df$time[!is.na(df$estimate)]))
    max.slope <- (y[1] - y[2]) / (x[2] - x[1]) # the maximum slope I can display (deltay  / delta x)
    new.slope <- tau * max.slope # convert tau to a slope I can plot: the max. slope of tau is not 1, but the max. given by the plot boundaries.
    col_tau <- ifelse(p.corrected < 0.05, 'red', ifelse(p_tau < 0.05, 'blue', 'black'))
    half.sg <- ifelse(new.slope > 0, half, half * -1) # this is to make sure that tau is 0 at the half of the graph.
    y.start <- mean(df$estimate[which(df$beep %in% c((half-5):(half+5)))], na.rm=T) - new.slope * (half - x.coord[1])
    y.end <- y.start + new.slope * (x.coord[2] - x.coord[1])
    plot <- plot + geom_segment(x = x.coord[1], y = y.start, xend = x.coord[2], yend = y.end, col = col_tau)
  }
  return(plot)
}
mw_fun <- function(plot, data, tp, winsize) {
  names(data) <- c('time', 'estimate')
  plot <- plot + theme_classic()
  # Add a line for the tau estimate
  plot <- tau_fun(plot, tp, data, mw = TRUE)
  l <- helpfun(tp, mw = TRUE, winsize = winsize)[[1]]
  no.tp <- helpfun(tp, mw = TRUE, winsize = winsize)[[2]]
  plot <- transition_line(plot, tp, l)
  # Add line depicting the EWS (either mean or SD)
  for (n in no.tp) { # add the line and its confidence intervals
    nw.data <- data[c(which(data$time == l[n]) : which(data$time == l[n+1])),]
    plot <- plot + geom_line(data=nw.data, aes(x = as.numeric(time), y = as.numeric(estimate)), inherit.aes=FALSE) }
  return(plot)
}
mw_plots <- function(df, tp) {
  p1 <- ggplot(data=df, aes(time, AR)) + theme_classic() + ylim(min(df$ll, na.rm=T), max(df$ul, na.rm=T)) + geom_hline(yintercept = -1) + geom_hline(yintercept = 1) + labs(y = 'autoregressive coefficient') + xlim(0, max(df$time, na.rm = T))
  p1 <- tp_fun(p1, tp, data = df[!is.na(df$AR),c('time', 'AR', 'll', 'ul')], mw = TRUE) # add vertical transition lines
  p1 <- tau_fun(p1, tp, df[c('time', 'AR')], mw = TRUE) # if the corrected p value of tau is significant, tau is printed in red. If the uncorrected p value of tau is significant, tau is printed in blue. If tau is not significant, black line.
  p2 <- ggplot(data=df, aes(time, mean)) + labs(y = 'mean') + xlim(0, max(df$time, na.rm = T))
  p2 <- mw_fun(p2, df[, c('time', 'mean')], tp) # function for the moving window estimates
  p3 <- ggplot(data=df, aes(time, sd)) + labs(y = 'standard deviation') + xlim(0, max(df$time, na.rm = T))
  p3 <- mw_fun(p3, df[, c('time', 'sd')], tp)
  return(list(p1, p2, p3))
}
gam_plots <- function(df, tp, orig.df) {
  
  p0 <- ggplot(data=df, aes(time, intercept)) + theme_classic() + xlim(0, nrow(df)) + ylim(min(df$ll.int, na.rm=T), max(df$ul.int, na.rm=T)) + geom_hline(yintercept = -1) + geom_hline(yintercept = 1) + labs(y = 'intercept')
  p0 <- tp_fun(p0, tp, data = df[!is.na(df$AR),c('time', 'intercept', 'll.int', 'ul.int')], mw = FALSE) # add vertical transition lines
  p0 <- p0 + geom_point(data = orig.df, aes(time, Var), col = 'grey', alpha = 0.5, inherit.aes = F)
  
  p1 <- ggplot(data=df, aes(time, AR)) + theme_classic() + xlim(0, nrow(df)) + ylim(min(df$ll, na.rm=T), max(df$ul, na.rm=T)) + geom_hline(yintercept = -1) + geom_hline(yintercept = 1) + labs(y = 'autoregressive coefficient')
  p1 <- tp_fun(p1, tp, data = df[!is.na(df$AR),c('time', 'AR', 'll', 'ul')], mw = FALSE) # add vertical transition lines
  #  p1 <- sign.change(p1, data = df[,c('time', 'AR', 'incr.AR', 'decr.AR', 'segment_incr.AR', 'segment_decr.AR')]) # add colored segments for significant in/decreases in EWS
  p1 <- p1 + geom_point(data = orig.df, aes(time, Var), col = 'grey', alpha = 0.5, inherit.aes = F)
  p2 <- ggplot(data=df, aes(time, mu)) + theme_classic() + ylim(min(df$ll.mu), max(df$ul.mu)) + xlim(0, nrow(df)) + labs(y = 'mu')
  p2 <- tp_fun(p2, tp, data = df[!is.na(df$mu),c('time', 'mu', 'll.mu', 'ul.mu')], mw = FALSE)
  #  p2 <- sign.change(p2, dat = df[,c('time', 'mu', 'incr.mu', 'decr.mu', 'segment_incr.mu', 'segment_decr.mu')])
  p2 <- p2 + geom_point(data = orig.df, aes(time, Var), col = 'grey', alpha = 0.5, inherit.aes = F)
  
  return(list(p0, p1, p2)) }

# EWS GAM Table Function
table.pp <- function(tabl, var, pp){
  tabl %>% 
    plyr::mutate(Checklist = cell_spec(Checklist, "html", color = "white", bold = T,
                                       background = ifelse(Checklist == "PASS", "#A7CB1B", "#E1004C") )) %>%
    kable(format = 'html', escape = FALSE, booktabs = TRUE, align = "l",  
          caption = paste("<b>", paste(var),"</b>", "GAM check for participant", paste(pp)) ) %>% kable_styling(full_width = F) 
  
}
