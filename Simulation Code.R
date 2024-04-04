########################################
# Win Probability MSI Code 
# This file will simulate survival data from the given scenario and calculate power for the following methods: Win Probability (WP), 
# Restricted Win Probability (RWP), RMST, Frequentist WR, Log-rank test. It will also calculate additional performance measures for WP and RWP/
#
# Author: Michelle Sonnenberger
# Date: 4/1/2024
# Scenario Information:
# Use R 4.0.5
# 
########################################


###################################################
# Simulation set up
###################################################

# Choose scenario: "ph", "weibull", "log.early", "log.late"
scenario <- "log.early"

# choose which True WP value to use: 1,2,3 with (0.606, 0.556, 0.526)
i <- 1


# Scenario Set-Up
set.wp.vals <- c(0.606, 0.556, 0.526)

# same across all scenarios
median.surv <- c(9)
#N <- c(seq(10, 45, 5), seq(50, 300, 50))
alpha <- 0.05
burn.N <- 500
iterations <- 5000
nsim <- 2000
recruitment.time <- 9
follow.up <- 12

set.wp <- set.wp.vals[i]

all.seeds <- 10*nsim*seq(13,25,by = 1)+1

if (scenario == "ph"){
    N <- c(100, 200, 400, 600, 800, 1200, 1600, 2400, 3200, 4200)/2
    HRs <- c(0.65, 0.8, 0.9)
    lambdas <- (-1/median.surv)*log(0.5)
    change <- 0
    shape.trt <- 1
    shape.cont <- 1
    start.seeds <- all.seeds[1:3]
    burn <- FALSE
    X <- HRs[i]
    startSeed <- start.seeds[i]
    slope = NA; int= NA
} else if (scenario == "log.early"){
  N <- c(100, 200, 400, 800, 1600, 2800, 4400, 6000, 7500)/2  
  start.seeds <- all.seeds[7:9]
  
  a.vals <- c(0.0102, 0.05101, 0.1409)
  change <- 0
  shape.trt <- 1
  shape.cont <- 1
  lambdas <- (-1/median.surv)*log(0.5)
  burn <- TRUE
  start.seeds <- all.seeds[7:9]
  X <- a.vals[i]
  startSeed <- start.seeds[i]
  slope = NA; int= NA
} else if (scenario == "log.late"){
  N <- c(100, 200, 400, 800, 1200, 1600, 2400, 3200, 4000, 5000)/2

  a.vals <- c(-0.0962, -0.0344, -0.0153)
  change <- 0
  shape.trt <- 1
  shape.cont <- 1
  lambdas <- (-1/median.surv)*log(0.5)
  burn <- TRUE
  start.seeds <- all.seeds[10:12]
  X <- a.vals[i]
  startSeed <- start.seeds[i]
  slope = NA; int= NA
  
}

# Function to run a single simulation
# seed: the seed for set.seed
# N: the sample size for a single group
# X: either HR for PH scenario, or value to be multiplied by lambda for WP
# lambda: the lambda value for all scenarios
# tau: the time set for RMST, WMST, WP with ties [default is follow-up time]
# slope: for log-linear scenario, the slope of the log(hazard) for treatment group
# int: for log-linear scenario, the intercept of the log(hazard) for treatment group
# shape.trt: the shape parameter for treatment group [default is 1]
# shape.cont: the shape parameter for control group [default is 1]
# scenario: the scenario which is being run:
#   ph: proportional hazard scenario. Exponential models
#   log.late: log-linear hazard for treatment. Exponential control. late difference
#   log.early: log-linear hazard for treatmnet. Exponential control. early difference
# rec.time: recruitment time
# follow.up.time: follow-up time
# burn: T/F if want to include burn-in for bayesian model [default is FALSE]
# burn.N: if burn=TRUE, length of burn-in [default is 500]
# iterations: number of iterations for bayesian model [default is 1000]
doSim <- function(seed, N, X = 1, lambda, rec.time, follow.up.time, tau = follow.up.time, slope = NA, int= NA, shape.trt = 1, shape.cont = 1,scenario = c("ph", "log.late", "log.early", "weibull", "null"), burn = FALSE, burn.N = 500, iterations = 1000){
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
  # seed: the seed used to generate the data [only used for error output]
  # the following parameters come from the bayesian weibull model
  # r.trt: is the shape parameter for treatment
  # r.cont: shape parameter for control
  # lambda.trt: scale parameter for treatment
  # lambda.cont: scale parameter for control
  wpIntegral <- function(r.trt, r.cont, lambda.trt, lambda.cont, lower, upper, seed){
    out <- tryCatch({
      int <- integrate(integrandParam, lower=lower, upper=upper, r.trt=r.trt, r.cont=r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont)
      return(int$value)
    },silent = TRUE , error = function(err){
      print("")
      print(paste("WP Integral failed: N= ", N.trt, "; seed = ", seed, sep = ""))
      vals <- data.frame(r.trt = r.trt, r.cont= r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, lower = lower, upper = upper)
      print(vals)
      return(NA)
    }
    )
    return(out)
  }
  
  # tau: the max of the integral
  # seed: the seed used to generate the data [only used for error output]
  # the following parameters come from the bayesian weibull model
  # r.trt: is the shape parameter for treatment
  # r.cont: shape parameter for control
  # lambda.trt: scale parameter for treatment
  # lambda.cont: scale parameter for control
  findWPTieVal <- function(r.trt, r.cont, lambda.trt, lambda.cont, tau, seed){
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
      print(paste("WP Tie Val Failed: N= ", N.trt, "; seed = ", seed, sep = ""))
      vals <- data.frame(r.trt = r.trt, r.cont= r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, tau=tau)
      print(vals)
      return(NA)
    }
    )
    return(out)
  }
  
  # seed: the seed used to generate the data [only used for error output]
  # the following parameters come from the bayesian weibull model
  # r.trt: is the shape parameter for treatment
  # r.cont: shape parameter for control
  # lambda.trt: scale parameter for treatment
  # lambda.cont: scale parameter for control
  findWPVal <- function(r.trt, r.cont, lambda.trt, lambda.cont, seed){
    out <- tryCatch({
      return(wpIntegral(r.trt=r.trt, r.cont=r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, lower = 0, upper = Inf, seed=seed))

    },silent = TRUE , error = function(err){
      print("")
      print(paste("WP Val Failed: N= ", N.trt, "; seed = ", seed, sep = ""))
      vals <- data.frame(r.trt = r.trt, r.cont= r.cont, lambda.trt=lambda.trt, lambda.cont=lambda.cont, tau=tau)
      print(vals)
      return(NA)
    }
    )
    return(out)
  }
  
  # N: sample size per group
  # t.trt/t.cont: Event time or NA if censored
  # t.cen.trt/t.cen.cont: event or censoring time
  # is.censored.trt/is.censored.cont: 0 if event, 1 if right censored
  # shape.trt/shape.cont: the shape for the initialization of the weibull model
  # tau: the time used for win probability with ties as last time
  # alpha: two-sided type I error [default 0.05]
  # seed: seed used; needed for error messages
  # lambda: the scale of the control weibull used for initialzing the model
  # X: scale for treatment is lambda*X used for initialzing the model
  # burn: TRUE if want to include burn-in
  # burn.N: length of burn-in
  # iterations: number iterations for the bayesian model
  findWinTie <- function(N, t.trt, t.cen.trt, t.cont, t.cen.cont, burn, burn.N, iterations, 
                         tau, alpha = 0.05, seed, lambda, X, is.censored.trt,
                         is.censored.cont, shape.trt, shape.cont){
    out <- tryCatch({
      
      # List with all elements required in model
      d.jags <- list(N.trt = N, t.trt = t.trt, t.cen.trt = t.cen.trt, N.cont = N, t.cont = t.cont, t.cen.cont = t.cen.cont,
                     is.censored.trt= is.censored.trt, is.censored.cont = is.censored.cont)
      
      # Initial values at the truth
      if (burn){
        i.jags <- function(){
          list(eta.trt = rnorm(1, 0, 0.0001), r.trt=1, eta.cont=rnorm(1, 0, 0.0001), r.cont = 1)
        }
      } else {
        i.jags <- function(){
          list(eta.trt = exp(lambda*X), r.trt=shape.trt, eta.cont=exp(lambda), r.cont = shape.cont)
        }
      }
      

      
      # Vector of monitored/saved parameters
      p.jags <- c("lambda.trt", "r.trt", "lambda.cont", "r.cont")
      
      # compile the JAGS model
      # MSI Code
      m <- jags.model(data = d.jags, file = "/TwoIndWeibull_JAGS_Paper.txt", inits = i.jags, n.chains  = 3)
 
      # burn-in: 1000 simulations
      # Not needed since we are initiating our parameters at the truth
      if(burn){
        update(m, burn.N)
      }
      #update(m, burn)
      
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
      
      return(c(wpTiePP, num.int.fail, wpTie.median, wpTieCI, wpPP, num.int.fail.wp, wp.median, wpCI))
      
    }, silent = TRUE , error = function(w){
      print("")
      print(paste("Find Win Tie failed: N= ", N, "; seed = ", seed, sep = ""))
      return(c(NA, NA))
    }
    )
    return(out)
  }
  
  
  library(KMsurv)
  library(survival)
  library(rjags)
  library(WR) #frequentist win ratio
  library(survRM2)
  library(clinfun)
  #library(survminer)
  library(mvtnorm)
  library(nph)
  library(survWMST)
  library(simsurv)
  set.seed(seed)
  all.out <- tryCatch({
    if (scenario == "ph"){
      u.c <- runif(N, max = 1)
      u.t <- runif(N, max = 1)
      
      #t.cont.event <- -(1/lambda)*log(1-u.c)
      #t.trt.event <- -(1/(HR*lambda))*log(1-u.t)
      
      t.cont.event <- qexp(u.c, rate = lambda, lower.tail = TRUE)
      t.trt.event <- qexp(u.t, rate = lambda*X, lower.tail = TRUE)
      
    } else if (scenario == "log.late"){
      get.time.log.linear.late <- function(u, a, lambda){
        if(u >= exp(0.5*lambda*(1/a))){
          return((1/a)*log(1-a*(1/lambda)*log(u)))
        } else{
          return((1/a)*(1+log(0.5))-2*(1/lambda)*log(u))
        }
      }
      
      u.c <- runif(N, max = 1)
      u.t <- runif(N, max = 1)
      
      t.cont.event <- qexp(u.c, rate = lambda, lower.tail = TRUE)
      t.trt.event <- unlist(lapply(u.t, get.time.log.linear.late, a = X, lambda = lambda))

    } else if (scenario == "log.early"){
      get.time.log.linear.early <- function(u, a, lambda){
        k <- 0.6
        change.pt <- -(1/a)*log(k)
        surv.trt.early4 <- function(t, a, lambda, k){
          change.pt <- -(1/a)*log(k)
          if (t< change.pt){
            return(exp(-(k*lambda*(1/a)*(exp(a*t)-1))))
          } else{
            return(exp(-(k*lambda*(1/a)*((1/k)-1)+lambda*(t-change.pt))))
          }
        }
        
        change.surv <- surv.trt.early4(t = change.pt, a=a, lambda=lambda, k = k)
        
        if(u > change.surv){ # i.e. t < change.pt
          return((1/a)*log(1-(1/k)*(1/lambda)*a*log(u)))
        } else{
          return(change.pt-(1/lambda)*(log(u)+k*lambda*(1/a)*((1/k)-1)))
        }
      }
      
      u.c <- runif(N, max = 1)
      u.t <- runif(N, max = 1)
      
      t.cont.event <- qexp(u.c, rate = lambda, lower.tail = TRUE)
      t.trt.event <- unlist(lapply(u.t, get.time.log.linear.early, a = X, lambda = lambda))
      
    }
    
    cen.time <- runif(N*2, min = follow.up.time, max = rec.time + follow.up.time)
    
    sim.dat <- data.frame(time.event = c(t.cont.event, t.trt.event), time.censor = cen.time, group = c(rep(0, N), rep(1, N)))
    sim.dat$status <- ifelse(sim.dat$time.event > sim.dat$time.censor, 0, 1)
    sim.dat$time <- ifelse(sim.dat$time.event > sim.dat$time.censor, sim.dat$time.censor, sim.dat$time.event)
    
    n.events <- sum(sim.dat$status)
    
    #log-rank test
    temp.log.rank <- survdiff(Surv(time, status) ~ group, data=sim.dat)
    
    # RMST
    temp.tau <- min(tau, min(max(sim.dat$time[which(sim.dat$group == 0)]),max(sim.dat$time[which(sim.dat$group == 1)])))
    temp.rmst <- rmst2(sim.dat$time, sim.dat$status, sim.dat$group, tau = temp.tau, covariates = NULL, alpha = alpha)
    
    # Frequentist WR
    sim.dat$id <- c(1:length(sim.dat$time))
    win.ratio <- WRrec(sim.dat$id, sim.dat$time, sim.dat$status, sim.dat$group, strata = NULL, naive = TRUE)
    
    # Win probability and Win Probability with Ties
    
    sim.dat$t <- ifelse(sim.dat$status == 1, sim.dat$time, NA)
    sim.dat$t.cen <- sim.dat$time
    sim.dat$is.censored <- ifelse(sim.dat$status == 1, 0, 1)
    
    dat.trt <- subset(sim.dat, sim.dat$group == 1)
    dat.cont <- subset(sim.dat, sim.dat$group == 0)
    
    wp <-   findWinTie(N = N, t.trt = dat.trt$t, t.cen.trt = dat.trt$t.cen, t.cont = dat.cont$t, t.cen.cont = dat.cont$t.cen, burn=burn, burn.N=burn.N, iterations = iterations, 
                                   tau = tau, alpha = 0.05, seed = seed, lambda = lambda, X = X, is.censored.trt = dat.trt$is.censored,
                                   is.censored.cont = dat.cont$is.censored, shape.trt = shape.trt, shape.cont = shape.cont)
    
    # WP with ties
    wp.Tie.int.fail <- wp[2]
    wp.Tie.median <- wp[3]
    wp.Tie.ci <- wp[4:5]
    wp.Tie.PP <- wp[1]
    
    # WP
    wp.int.fail <- wp[7]
    wp.median <- wp[8]
    wp.ci <- wp[9:10]
    wp.PP <- wp[6]


      # Calculate Power
      # One-sided error
      
      # log-rank
      temp.log.rank.pval <- pnorm(sqrt(temp.log.rank$chisq), lower.tail = FALSE)
      log.rank.sig <- ifelse(temp.log.rank.pval< alpha/2, 1, 0)
    
      # RMST
      est <- temp.rmst$unadjusted.result[1,1]
      se <- abs((est - temp.rmst$unadjusted.result[1,3])/qnorm(1 - alpha/2))
      rmst.pval <- pnorm(abs(est)/se, lower.tail = FALSE) 
      rmst.sig <- ifelse(rmst.pval< alpha/2, 1, 0)
      
      # Frequentist win ratio
      wr.pval <- pnorm(abs(win.ratio$log.WR/win.ratio$se), lower.tail = FALSE)
      wr.sig <- ifelse(wr.pval< alpha/2, 1, 0)
      
      # WP with Ties Determine significance
      wpTie.sig <- ifelse(wp[1] > 1-alpha/2, 1, 0)

      # WP Determine significance
      wp.sig <- ifelse(wp[6] > 1-alpha/2, 1, 0)
    
    
    # create dataset that will save information from each iteration
    save.dat <- sim.dat
    save.dat$seed <- seed
    save.dat$follow.up <- follow.up.time
    save.dat$recruitment <- rec.time
    save.dat$tau <- tau
    
    return(list(log.rank.sig=log.rank.sig, rmst.sig=rmst.sig, FWR.sig=wr.sig, RWP.sig=wpTie.sig, n.events=n.events, RWP.int.fail=wp.Tie.int.fail, 
                RWP.median=wp.Tie.median, low.RWP.ci = wp.Tie.ci[1], up.RWP.ci= wp.Tie.ci[2], seed = seed, tau = temp.tau, wp.sig = wp.sig, 
                wp.int.fail = wp.int.fail, wp.median = wp.median, wp.ci = wp.ci, RWP.PP = wp.Tie.PP, wp.PP = wp.PP))

  },warning = function(w){
    #nFail <- nFail + 1
    print("")
    print(paste("Simulation failed: N= ", N, "; seed = ", seed, sep = ""))
    return(NULL)
  }
  
  )
  return (all.out)
}
  
single.sim <- doSim(seed = iter, N = n, X = X, lambda = lambda, rec.time = recruitment.time, follow.up.time = follow.up, 
                        tau = follow.up, slope = slope, int= int, shape.trt = shape.trt, shape.cont = shape.cont,
                        scenario = scenario, burn = burn, burn.N = burn.N, iterations = iterations) 


# Function that finds simulation performance measures: bias, mse, ci coverage, and ci width from simulation results (estimate, ci, and true value)
simPerfMeas <- function(estimate, ci.low, ci.up, truth){
  bias <- mean(estimate) - truth
  
  mse <- (1/length(estimate))*sum((estimate-truth)^2)
  
  ci.cov <- mean( ci.low < truth & ci.up > truth) 
  
  ci.width <- mean(ci.up - ci.low)
  
  return(c(bias, mse, ci.cov, ci.width))
}

# calculate true WP
calcTrueWP <- function(lambda, X, shape.trt = 1, shape.cont = 1, scenario){
  integrandParamTruth <- function(t, lambda.trt, lambda.cont, shape.trt = 1, shape.cont = 1){
    r.scale.trt <- lambda.trt^(-1/shape.trt)
    r.scale.cont <- lambda.cont^(-1/shape.cont)
    
    return(exp(pweibull(t, shape = shape.trt, scale = r.scale.trt, lower.tail = FALSE, log.p=TRUE) + 
                 dweibull(t, shape = shape.cont, scale = r.scale.cont, log = TRUE)))
  }
  
  if(scenario == "log.late"){
    surv.trt.start.late <- function(t, a, lambda){
      return(exp(-lambda * (1/a)*(exp(a*t)-1)))
    }
    
    surv.trt.end.late <- function(t, a, lambda){
      # at log(0.5)/a h(t) 0.5lambda and the survival function change
      return(exp(-0.5*lambda*(t-(1/a)*(log(0.5)+1))))
    }
    
    integrandParam.start.late <- function(t, lambda, a){
      return(exp(log(surv.trt.start.late(t = t, a = a, lambda = lambda)) + 
                   dexp(t, rate=lambda, log = TRUE)))
    }
    
    integrandParam.end.late <- function(t, lambda, a){
      return(exp(log(surv.trt.end.late(t = t, a = a, lambda = lambda)) + 
                   dexp(t, rate=lambda, log = TRUE)))
    }
    
    find.WP.Late.Effect <- function(lambda, a){
      int.start <- integrate(integrandParam.start.late, lower=0, upper=(1/a)*log(0.5), a=a, lambda = lambda)
      int.end <- integrate(integrandParam.end.late, lower=(1/a)*log(0.5), upper=Inf, a=a, lambda = lambda)
      return(int.start$value + int.end$value)
    }
    
    return(find.WP.Late.Effect(lambda = lambda, a = X))
  } 
  
  if(scenario == "log.early"){
    find.WP.Early.Effect4 <- function(lambda, a, k){
      surv.trt.early.start4 <- function(t, a, lambda, k){
        return(exp(-(k*lambda*(1/a)*(exp(a*t)-1))))
      }
      
      surv.trt.early.end4 <- function(t, a, lambda, k){
        # at 36 h(t) =lambda and the survival function change
        change.pt <- -(1/a)*log(k)
        return(exp(-(k*lambda*(1/a)*((1/k)-1)+lambda*(t-change.pt))))
      }
      
      integrandParam.start.early4 <- function(t, lambda, a, k){
        return(exp(log(surv.trt.early.start4(t = t, a = a, lambda = lambda, k = k)) + 
                     dexp(t, rate=lambda, log = TRUE)))
      }
      
      integrandParam.end.early4 <- function(t, lambda, a, k){
        return(exp(log(surv.trt.early.end4(t = t, a = a, lambda = lambda, k = k)) + 
                     dexp(t, rate=lambda, log = TRUE)))
      }
      
      change.pt <- (1/a)*log(1/k)
      int.start <- integrate(integrandParam.start.early4, lower=0, upper=change.pt, a=a, lambda = lambda, k = k)
      int.end <- integrate(integrandParam.end.early4, lower=change.pt, upper=Inf, a=a, lambda = lambda, k = k)
      return(int.start$value + int.end$value)
    }
    
    return(find.WP.Early.Effect4(lambda = lambda, a = X, k = 0.6))
  } else{
    int.true <- integrate(integrandParamTruth, lower=0, upper=Inf, lambda.trt=lambda*X, lambda.cont=lambda, shape.trt = shape.trt, shape.cont = shape.cont)
    true.wp <- int.true$value
  }
  return(true.wp)
}

# calculate True RWP for included scenarios
calcTrueWPTies <- function(lambda, X, shape.trt = 1, shape.cont = 1, scenario, tau){
  integrandParamTruth <- function(t, lambda.trt, lambda.cont, shape.trt = 1, shape.cont = 1){
    r.scale.trt <- lambda.trt^(-1/shape.trt)
    r.scale.cont <- lambda.cont^(-1/shape.cont)
    
    return(exp(pweibull(t, shape = shape.trt, scale = r.scale.trt, lower.tail = FALSE, log.p=TRUE) + 
                 dweibull(t, shape = shape.cont, scale = r.scale.cont, log = TRUE)))
  }
  
  if(scenario == "log.late"){
    surv.trt.start.late <- function(t, a, lambda){
      return(exp(-lambda * (1/a)*(exp(a*t)-1)))
    }
    
    surv.trt.end.late <- function(t, a, lambda){
      # at log(0.5)/a h(t) 0.5lambda and the survival function change
      return(exp(-0.5*lambda*(t-(1/a)*(log(0.5)+1))))
    }
    
    integrandParam.start.late <- function(t, lambda, a){
      return(exp(log(surv.trt.start.late(t = t, a = a, lambda = lambda)) + 
                   dexp(t, rate=lambda, log = TRUE)))
    }
    
    integrandParam.end.late <- function(t, lambda, a){
      return(exp(log(surv.trt.end.late(t = t, a = a, lambda = lambda)) + 
                   dexp(t, rate=lambda, log = TRUE)))
    }
    
    find.WPTies.Late.Effect <- function(lambda, a){
      change.pt <- (1/a)*log(0.5)
      if (change.pt < tau){
        # tau > change.pt in survival function
        # integral is broken up into two pieces around changept
        
        # P(Ya>Yc, Yc<=tau)
        # int [0 to tau] [1-Fa]fc
        int.start <- integrate(integrandParam.start.late, lower=0, upper=change.pt, a=a, lambda = lambda)
        int.end <- integrate(integrandParam.end.late, lower=change.pt, upper=tau, a=a, lambda = lambda)
        p1.int <- int.start$value + int.end$value
        
        # treatment survival at tau will follow end function
        # 0.5P(min(Ya,Yc)>tau)
        # [Sa(tau)][1-Fc(tau)]
        p2.true <- 0.5*surv.trt.end.late(t = tau, a = a, lambda = lambda)*
          pexp(tau, rate = lambda, lower.tail = FALSE)
        wp.val <- p1.int + p2.true
      } else{
        # tau < change.pt in survival function
        # integral is one function
        
        # P(Ya>Yc, Yc<=tau)
        # int [0 to tau] [1-Fa]fc
        p1.int <- integrate(integrandParam.start.late, lower=0, upper=tau, a=a, lambda = lambda)$value
        
        # treatment survival at tau will follow start function
        # 0.5P(min(Ya,Yc)>tau)
        # [Sa(tau)][1-Fc(tau)]
        p2.true <- 0.5*surv.trt.start.late(t = tau, a = a, lambda = lambda)*
          pexp(tau, rate = lambda, lower.tail = FALSE)
        wp.val <- p1.int + p2.true
      }
      
      return(wp.val)
    }
    
    return(find.WPTies.Late.Effect(lambda = lambda, a = X))
  } 
  
  if(scenario == "log.early"){
    surv.trt.early.start4 <- function(t, a, lambda, k){
      return(exp(-(k*lambda*(1/a)*(exp(a*t)-1))))
    }
    
    surv.trt.early.end4 <- function(t, a, lambda, k){
      # at 36 h(t) =lambda and the survival function change
      change.pt <- -(1/a)*log(k)
      return(exp(-(k*lambda*(1/a)*((1/k)-1)+lambda*(t-change.pt))))
    }
    
    integrandParam.start.early4 <- function(t, lambda, a, k){
      return(exp(log(surv.trt.early.start4(t = t, a = a, lambda = lambda, k = k)) + 
                   dexp(t, rate=lambda, log = TRUE)))
    }
    
    integrandParam.end.early4 <- function(t, lambda, a, k){
      return(exp(log(surv.trt.early.end4(t = t, a = a, lambda = lambda, k = k)) + 
                   dexp(t, rate=lambda, log = TRUE)))
    }
    
    find.WPTies.Early.Effect4 <- function(lambda, a){
      start.hr = 0.6
      change.pt <- -(1/a)*log(start.hr)
      if (change.pt < tau){
        # tau > change.pt in survival function
        # integral is broken up into two pieces around changept
        
        # P(Ya>Yc, Yc<=tau)
        # int [0 to tau] [1-Fa]fc
        int.start <- integrate(integrandParam.start.early4, lower=0, upper=change.pt, a=a, lambda = lambda, k = start.hr)
        int.end <- integrate(integrandParam.end.early4, lower=change.pt, upper=tau, a=a, lambda = lambda, k = start.hr)
        p1.int <- int.start$value + int.end$value
        
        # treatment survival at tau will follow end function
        # 0.5P(min(Ya,Yc)>tau)
        # [Sa(tau)][1-Fc(tau)]
        p2.true <- 0.5*surv.trt.early.end4(t = tau, a = a, lambda = lambda, k = start.hr)*
          pexp(tau, rate = lambda, lower.tail = FALSE)
        wp.val <- p1.int + p2.true
      } else{
        # tau < change.pt in survival function
        # integral is one function
        
        # P(Ya>Yc, Yc<=tau)
        # int [0 to tau] [1-Fa]fc
        p1.int <- integrate(integrandParam.start.early4, lower=0, upper=tau, a=a, lambda = lambda, k = start.hr)$value
        
        # treatment survival at tau will follow start function
        # 0.5P(min(Ya,Yc)>tau)
        # [Sa(tau)][1-Fc(tau)]
        p2.true <- 0.5*surv.trt.early.start4(t = tau, a = a, lambda = lambda, k = start.hr)*
          pexp(tau, rate = lambda, lower.tail = FALSE)
        wp.val <- p1.int + p2.true
      }
      
      return(wp.val)
    }
    
    return(find.WPTies.Early.Effect4(lambda = lambda, a = X))
    
  } else{
    # P(Ya>Yc, Yc<=tau)
    # int [0 to tau] [1-Fa]fc
    p1.int <- integrate(integrandParamTruth, lower=0, upper=tau, lambda.trt=lambda*X, lambda.cont=lambda, shape.trt = shape.trt, shape.cont = shape.cont)
    p1.true <- p1.int$value
    
    # 0.5P(min(Ya,Yc)>tau)
    # [1-Fa(tau)][1-Fc(tau)]
    p2.true <- 0.5*pweibull(tau, shape = shape.trt, scale = (lambda*X)^(-1/shape.trt), lower.tail = FALSE)*
      pweibull(tau, shape = shape.cont, scale = lambda^(-1/shape.cont), lower.tail = FALSE)
    
    true.wp <- p1.true + p2.true
    
    return(true.wp)
  }
  
}