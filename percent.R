percent <-  function(data, B, t, m, mu.true) {
  # Bootstrapping
  resamples <-
    replicate(B, sample(data, size = length(data), replace = TRUE))
  
  # Analyze each bootstrapped resamples
  temp <- apply(resamples, 2, function(x)
    prelim(x, t, m, mu.true))
  
  # Results from the generating model ---------------------------------------
  
  # Mean inter-event times
  mu.gen <- sapply(temp, '[[', 15)
  
  # Pr(X < mu)
  pr.gen <- sapply(temp, '[[', 16)
  
  # Hazard rates
  haz.gen <- sapply(temp, '[[', 17)
  
# Percentile bootstrapped C.I for the estimates from generating model -----
  
  # Mean inter-event times
  mu.gen.ci <-  quantile(mu.gen, probs = c(0.025, 0.5, 0.975))
  
  # Pr(X < mu)
  pr.gen.ci <-  quantile(pr.gen, probs = c(0.025, 0.5, 0.975))
  
  # Hazard rates
  haz.gen.ci <-
    apply(haz.gen, 1, function(x)
      quantile(x, probs = c(0.025, 0.5, 0.975)))
  
  # Results from each single best model -------------------------------------
  
  # Mean inter-event times
  mu.best <- sapply(temp, '[[', 12)
  
  # Pr(X < mu)
  pr.best <- sapply(temp, '[[', 13)
  
  # Hazard rates
  haz.best <- sapply(temp, '[[', 14)
  
  # Bootstapped propotion weights -------------------------------------------
  
  # Index for a single best model
  id <- sapply(temp, '[[', 11)
  
  # Bootstarpoed proportion weights
  prop.weights <- sapply(1:6, function(i)
    length(which(i == id)) / B)
  
  # Model-averaged percentile bootstrapped confidence intervals --------------
  
  # Mean inter-event times
  mu.best.ci <-  quantile(mu.best, probs = c(0.025, 0.5, 0.975))
  
  # Pr(X < mu)
  pr.best.ci <-  quantile(pr.best, probs = c(0.025, 0.5, 0.975))
  
  # Hazard rates
  haz.best.ci <-
    apply(haz.best, 1, function(x)
      quantile(x, probs = c(0.025, 0.5, 0.975)))
  
  # Output ------------------------------------------------------------------
  
  out <-
    list(
      'prop.weights' = prop.weights,
      'mu.gen' = as.numeric(mu.gen.ci[2]),
      'mu.gen.lower' = as.numeric(mu.gen.ci[1]),
      'mu.gen.upper' = as.numeric(mu.gen.ci[3]),
      'mu.best' = as.numeric(mu.best.ci[2]),
      'mu.best.lower' = as.numeric(mu.best.ci[1]),
      'mu.best.upper' = as.numeric(mu.best.ci[3]),
      'pr.gen' = as.numeric(pr.gen.ci[2]),
      'pr.gen.lower' = as.numeric(pr.gen.ci[1]),
      'pr.gen.upper' = as.numeric(pr.gen.ci[3]),
      'pr.best' = as.numeric(pr.best.ci[2]),
      'pr.best.lower' = as.numeric(pr.best.ci[1]),
      'pr.best.upper' = as.numeric(pr.best.ci[3]),
      'haz.gen' = as.numeric(haz.gen.ci[2,]),
      'haz.gen.lower' = as.numeric(haz.gen.ci[1,]),
      'haz.gen.upper' = as.numeric(haz.gen.ci[3,]),
      'haz.best' = as.numeric(haz.best.ci[2,]),
      'haz.best.lower' = as.numeric(haz.best.ci[1,]),
      'haz.best.upper' = as.numeric(haz.best.ci[3,])
    )
  
  return(out)
  
}


