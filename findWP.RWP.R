# t: time
# the following parameters come from the bayesian weibull model
# r.trt: is the shape parameter for treatment
# r.cont: shape parameter for control
# lambda.trt: scale parameter for treatment
# lambda.cont: scale parameter for control
integrandParam <- function(t, r.trt, r.cont, lambda.trt, lambda.cont){
  # shape.trt = r.trt
  # scale.trt = lambda.trt^(-1/r.trt)
  
  # shape.cont = r.cont
  # scale.cont = lambda.cont^(-1/r.cont)
  
  # We want: exp(log(1-Fa)+log(fc))
  # log(1-Fa)=pweibull(, lower.tail = FALSE, log.p = TRUE)
  # log(fc)= dweibull(, log = TRUE)
  
  return(exp(pweibull(t, shape = r.trt, scale = lambda.trt^(-1/r.trt), lower.tail = FALSE, log.p=TRUE) + 
               dweibull(t, shape = r.cont, scale = lambda.cont^(-1/r.cont), log = TRUE)))
}

# lower: lower bound of integral
# upper: upper bound of integral
# the following parameters come from the bayesian weibull model
# r.trt: is the shape parameter for treatment
# r.cont: shape parameter for control
# lambda.trt: scale parameter for treatment
# lambda.cont: scale parameter for control
wpIntegral <- function(r.trt, r.cont, lambda.trt, lambda.cont, lower, upper){
  out <- tryCatch({
    int <- integrate(integrandParam, lower=lower, upper=upper, r.trt=r.trt, r.cont=r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont)
    return(int$value)
  },silent = TRUE , error = function(err){
    print("")
    print("WP Integral failed")
    vals <- data.frame(r.trt = r.trt, r.cont= r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, lower = lower, upper = upper)
    print(vals)
    return(NA)
  }
  )
  return(out)
}

# tau: the max of the integral
# the following parameters come from the bayesian weibull model
# r.trt: is the shape parameter for treatment
# r.cont: shape parameter for control
# lambda.trt: scale parameter for treatment
# lambda.cont: scale parameter for control
findWPTieVal <- function(r.trt, r.cont, lambda.trt, lambda.cont, tau){
  out <- tryCatch({
    # r2 is the ranking function 2
    # WP2 = 1*P(r2=1)+0.5*P(r2=0.5)
    # WP2 = P(Ya>Yc, Yc<=tau) + 0.5*P(min(Ya,Yc)>tau)
    
    # p1 = P(Ya>Yc, Yc<=tau)
    # int[ 0 to tau] [1-Fa(t)]fc(t)
    p1 <- wpIntegral(r.trt=r.trt, r.cont=r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, lower = 0, upper = tau, seed=seed)
    
    # p2 = 0.5*P(min(Ya,Yc)>tau)
    # 0.5*[1-Fa(tau)][1-Fc(tau)]
    p2  <- 0.5*pweibull(tau, shape = r.trt, scale = lambda.trt^(-1/r.trt), lower.tail = FALSE)*pweibull(tau, shape = r.cont, scale = lambda.cont^(-1/r.cont), lower.tail = FALSE)
    
    return(p1+p2)
  },silent = TRUE , error = function(err){
    print("")
    print("WP Tie Val Failed")
    vals <- data.frame(r.trt = r.trt, r.cont= r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, tau=tau)
    print(vals)
    return(NA)
  }
  )
  return(out)
}

# the following parameters come from the bayesian weibull model
# r.trt: is the shape parameter for treatment
# r.cont: shape parameter for control
# lambda.trt: scale parameter for treatment
# lambda.cont: scale parameter for control
findWPVal <- function(r.trt, r.cont, lambda.trt, lambda.cont){
  out <- tryCatch({
    return(wpIntegral(r.trt=r.trt, r.cont=r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, lower = 0, upper = Inf, seed=seed))
    
  },silent = TRUE , error = function(err){
    print("")
    print("WP Val Failed")
    vals <- data.frame(r.trt = r.trt, r.cont= r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, tau=tau)
    print(vals)
    return(NA)
  }
  )
  return(out)
}

# This function will use right-censored time-to-event data and calculate win probability (WP) and restricted win probability (RWP) 
# using a Bayesian Weibull Model
############################################
# time: the event or censoring time
# group: a 0/1 indicator variable for group
# event: indicator for if time is event time (1) or censoring time (0)
# tau: the time used for win probability with ties as last time
# alpha: two-sided type I error [default 0.05]
# seed: seed used; needed for error messages
# burn: TRUE if want to include burn-in
# burn.N: length of burn-in
# iterations: number iterations for the bayesian model

# returned values
# RWP.PP: posterior probability for RWP summary measure
# RWP.Median: median for RWP summary measure
# RWP.CI: 95% credible interval for RWP summary measure
# WP.PP: posterior probability for WP summary measure
# WP.Median: median for WP summary measure
# WP.CI: 95% credible interval for WP summary measure
findWP.RWP <- function(time, group, event, burn.N = 1000, iterations = 5000, 
                       tau, alpha = 0.05, seed = 1, burn = TRUE){
  
  set.seed(seed)
  dat <- data.frame(time = time, group = group, event = event)
  dat$t <- ifelse(dat$event == 1, dat$time, NA)
  dat$t.cen <- dat$time
  dat$is.censored <- ifelse(dat$event == 1, 0, 1)
  
  dat.trt <- subset(dat, dat$group == 1)
  dat.cont <- subset(dat, dat$group == 0)
  
  N.trt = length(dat.trt$t); t.trt = dat.trt$t; t.cen.trt = dat.trt$t.cen; is.censored.trt = dat.trt$is.censored
  N.cont = length(dat.cont$t); t.cont = dat.cont$t; t.cen.cont = dat.cont$t.cen; is.censored.cont = dat.cont$is.censored

  out <- tryCatch({
    
    # List with all elements required in model
    d.jags <- list(N.trt = N.trt, t.trt = t.trt, t.cen.trt = t.cen.trt, N.cont = N.cont, t.cont = t.cont, t.cen.cont = t.cen.cont,
                   is.censored.trt= is.censored.trt, is.censored.cont = is.censored.cont)
    
    # initialize values
    i.jags <- function(){
      list(eta.trt = rnorm(1, 0, 0.0001), r.trt=1, eta.cont=rnorm(1, 0, 0.0001), r.cont = 1)
    }
    
    # Vector of monitored/saved parameters
    p.jags <- c("lambda.trt", "r.trt", "lambda.cont", "r.cont")
    
    # compile the JAGS model
    m <- jags.model(data = d.jags, file = "/TwoIndWeibull_JAGS_Paper.txt", inits = i.jags, n.chains  = 3)
    
    # Burn-in
    if(burn){
      update(m, burn.N)
    }
    
    
    # run model for 50,000 and keep 1 in 10 for proper thinning
    res <- coda.samples(m, variable.names = p.jags, n.iter = iterations, n.tim = 10)
    
    result <- as.mcmc(do.call(rbind, res))
    
    lambda.cont <- result[,1]
    lambda.trt <- result[,2]
    r.cont <- result[,3]
    r.trt <- result[,4]
    
    # Calculate win probability with ties
    wp.tie <-  mapply(findWPTieVal, r.trt=r.trt, r.cont=r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, tau=tau, seed=seed)
    num.int.fail <- sum(length(which(is.na(wp.tie))))
    wp.tie <- na.omit(wp.tie)
    
    # Calculate Posterior probability
    wpTiePP <- (1/length(wp.tie))*sum(1*(wp.tie>0.5))
    wpTie.median <- quantile(wp.tie, probs=c(0.5))
    wpTieCI <- quantile(wp.tie, probs=c(0.025,0.975))
    
    # Calculate Win Proability
    wp <- mapply(findWPVal, r.trt=r.trt, r.cont=r.cont, lambda.trt=lambda.trt, lambda.cont= lambda.cont, seed = seed)
    num.int.fail.wp <- sum(length(which(is.na(wp))))
    wp <- na.omit(wp)
    
    # Calculate Posterior probability
    wpPP <- (1/length(wp))*sum(1*(wp>0.5))
    wp.median <- quantile(wp, probs=c(0.5))
    wpCI <- quantile(wp, probs=c(0.025,0.975))
    
    paste("WP had ", num.int.fail.wp, "draws with failed integration.", " RWP had ", num.int.fail, "draws with failed integration.", sep = "")
    paste("")
    return(list(RWP.PP = wpTiePP, RWP.Median = wpTie.median, RWP.CI = wpTieCI, WP.PP = wpPP, WP.Median = wp.median, WP.CI = wpCI))
    
  }, silent = TRUE , error = function(w){
    print("")
    print(paste("Find Win Tie failed: N= ", N, "; seed = ", seed, sep = ""))
    return(c(NA, NA))
  }
  )
  return(out)
}