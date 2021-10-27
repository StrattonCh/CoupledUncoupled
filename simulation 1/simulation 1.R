# clear workspace
rm(list = ls())

# packages
packs <- c("dplyr", "nimble", "htmltools", "ggplot2", "sf", "Rcpp", "RcppArmadillo", "inline", "mvtnorm", "readr", "parallel", "xtable", "rstan", "coda", "vegan", "tidyverse", "lubridate", "bcmaps", "knitr", "kableExtra", "ggalluvial")
sapply(packs, require, character.only = T)
rm(packs)

# convenience
`%notin%` <- Negate("%in%")
options(mc.cores = parallel::detectCores())

# from model 0
load("simulation 1/parameter_estimates.rdata")

# functions
# simulate data
sim_dat <- function(
  nsites = 100, nspecies = 8, nvisits = 4, seed = NULL, 
  psi = runif(nspecies, .4, .9), 
  lambda = abs(rnorm(nspecies, 0, 100)),
  theta = t(apply(18*diag(nspecies)+2, 1, function(x) rdirch(1, x))),
  phi = runif(nspecies),
  ambig_frac = .75
){
  
  # optional seed
  if(!is.null(seed)) set.seed(seed)
  
  # build empty df
  df <- tibble(
    site = rep(1:nsites, each = nspecies * nspecies * nvisits),
    visit = rep(rep(1:nvisits, each = nspecies * nspecies), nsites),
    true_spp = rep(1:nspecies, nsites * nvisits * nspecies),
    id_spp = rep(rep(1:nspecies, each = nspecies), nsites * nvisits)
  )
  
  df2 <- left_join(
    df, 
    tibble(
      true_spp = 1:nspecies,
      lambda = lambda,
      psi = psi
    ), by = "true_spp"
  ) %>%
    left_join(
      ., 
      tibble(
        true_spp = rep(1:nspecies, nspecies),
        id_spp = rep(1:nspecies, each = nspecies),
        theta = c(theta)
      ), by = c("true_spp", "id_spp")
    )
  
  # latent z state
  df3 <- df2 %>%
    select(site, true_spp, psi) %>%
    distinct() %>%
    mutate(z = rbinom(n(), size = 1, prob = .data$psi)) %>%
    left_join(
      df2, 
      ., 
      by = c("site", "true_spp", "psi")
    ) %>%
    mutate(
      count = rpois(n(), z * lambda * theta)
    )
  full_df <- df3
  
  # allocate some unambiguous calls
  if(all(phi == 1)){
    df4 <- df3 %>%
      select(site, visit) %>%
      distinct %>%
      group_by(site) %>%
      sample_frac(ambig_frac) %>% 
      mutate(true_spp = NA) %>%
      ungroup()
    df5 <- bind_rows(
      anti_join(df3, df4, by = c("site", "visit")),
      inner_join(
        df3 %>% select(-true_spp),
        df4,
        by = c("site", "visit")
      )
    ) %>%
      arrange(site, visit) %>%
      group_by(site, visit, true_spp, id_spp) %>%
      summarize(count = sum(count)) %>%
      mutate(type = ifelse(is.na(true_spp), "ambiguous", "unambiguous"))
    
    out <- df5
  } else{
    # cases 1 and 2 first
    df3point5 <- df3 %>%
      filter(visit %in% c(1:(nvisits/2)))
    df4 <- df3point5 %>%
      select(site, visit) %>%
      distinct %>%
      group_by(site) %>%
      sample_frac(.5) %>% 
      mutate(true_spp = NA) %>%
      ungroup()
    df5 <- bind_rows(
      anti_join(df3point5, df4, by = c("site", "visit")),
      inner_join(
        df3point5 %>% select(-true_spp),
        df4,
        by = c("site", "visit")
      )
    ) %>%
      arrange(site, visit) %>%
      group_by(site, visit, true_spp, id_spp) %>%
      summarize(count = sum(count)) %>%
      mutate(type = ifelse(is.na(true_spp), "ambiguous", "unambiguous")) %>%
      mutate(lik_type = ifelse(is.na(true_spp), "2", "1")) %>%
      mutate(phi = ifelse(type == "ambiguous", 0, 1))
    
    # case 3
    df6 <- df3 %>%
      filter(visit %in% c((nvisits/2 + 1):nvisits))
    
    if(all(phi == "k") | (length(phi) == nspecies)){
      if(all(phi == "k")) phi <- runif(nspecies)
      
      df7 <- df6 %>%
        left_join(
          ., 
          tibble(true_spp = 1:nspecies, phi = phi),
          by = "true_spp"
        ) %>%
        mutate(
          count_unambig = round(phi * count),
          count_ambig = count - count_unambig
        )
      
      df8 <- df7 %>%
        select(site, visit, true_spp, id_spp, count_ambig) %>%
        group_by(site, visit, id_spp) %>%
        summarize(ambig_total = sum(count_ambig)) %>%
        mutate(
          type = "ambiguous",
          lik_type = "3b"
        ) %>%
        mutate(
          count = ambig_total
        ) %>%
        select(-ambig_total)
      
      df9 <- df7 %>%
        select(site, visit, true_spp, id_spp, count, phi, count_unambig) %>%
        mutate(type = "unambiguous", lik_type = "3a") %>%
        mutate(count = count_unambig) %>%
        select(-count_unambig)
      
      df10 <- bind_rows(
        df8, df9
      ) %>% arrange(site, visit, true_spp, id_spp)
    }
    
    
    out <- bind_rows(
      df5, df10
    ) %>% arrange(site, visit, true_spp, id_spp)
    
  }
  
  
  out2 <- list(df = out, params = list(psi = psi, theta = theta, lambda = lambda, phi = phi))
  return(out2)
}

# fit model
cd_model_nophi <- function(cd_data, save = T, filename = "test.rds", parallel = F, niter = 5000, alpha0, inits_list){
  
  fit_model <- function(seed = 1, code, data, constants, inits_list, niter, nchains, thin = 1){
    library(nimble)
    
    # handle inits in each chain
    inits_ <- list()
    inits_$psi <- inits_list$psi + runif(10, -.05, .05)
    inits_$theta <- inits_list$theta + matrix(runif(100, -.05, .05), 10, 10)
    inits_$theta <- abs(inits_$theta)  / rowSums(abs(inits_$theta))
    inits_$lambda <- inits_list$lambda + runif(10, -inits_list$lambda, inits_list$lambda)
    
    # R model
    model <- nimbleModel(code, constants, data)
    
    # C model
    model_c <- compileNimble(model)
    
    # R mcmc
    model_conf <- configureMCMC(model)
    
    # R mcmc
    mcmc <- buildMCMC(model_conf)
    
    # C mcmc
    mcmc_c <- compileNimble(mcmc, project = model_c)
    
    # run model
    out <- runMCMC(
      mcmc_c, 
      niter = niter, 
      nchains = nchains, 
      thin = thin, 
      init = inits_,
      setSeed = seed
    )
    
    # out
    return(out)
  }
  
  # housekeeping
  nspecies <- length(unique(cd_data$id_spp))
  
  # no covariates
  site_df <- cd_data %>%
    ungroup %>%
    select(
      site
    ) %>%
    distinct()
  
  lik1_df <- cd_data %>%
    filter(type == "unambiguous")
  lik2_df <- cd_data %>%
    filter(type == "ambiguous")
  
  code <- nimbleCode({
    # priors
    for(k in 1:nspecies){
      psi[k] ~ dbeta(1, 1)
      lambda[k] ~ T(dnorm(0, sd = 100), 0, Inf)
      theta[k, 1:nspecies] ~ ddirch(alpha = alpha0[k, 1:nspecies])
    }
    
    # likelihood - site level occupancy
    for(site in 1:nsites){
      for(k in 1:nspecies){
        z[site,k] ~ dbern(psi[k])
      }
    }
    
    # likelihood 1
    for(row in 1:n1){
      y1[row] ~ dpois(
        z[site1[row],true_spp1[row]] *
          lambda[true_spp1[row]] *
          theta[true_spp1[row],id_spp1[row]]
      )
    }
    
    # likelihood 2
    for(row in 1:n2){
      for(true in 1:nspecies){
        mean2[row, true] <- z[site2[row],true]*lambda[true]*theta[true,id_spp2[row]]
      }
      y2[row] ~ dpois(sum(mean2[row, 1:nspecies]))
    }
  })
  
  # model housekeeping
  # init_func <- function(){
  #   out <- list(
  #     psi = runif(10),
  #     lambda = abs(rnorm(10, 0, 100)),
  #     theta = t(apply(alpha0, 1, function(x) rdirch(1, x))),
  #     z = matrix(1, nrow = 100, ncol = 10)
  #   )
  #   
  #   return(out)
  # }
  
  if(parallel) {
    library(parallel)
    this_cluster <- makeCluster(3)
    fit <- parLapply(
      cl = this_cluster,
      X = 1:3,
      fun = fit_model,
      code = code,
      data = list(
        # likelihood 1
        y1 = lik1_df$count,
        
        # likelihood 2
        y2 = lik2_df$count,
        
        # theta prior
        alpha0 = alpha0
      ),
      constants = list(
        nsites = nrow(site_df),
        nspecies = length(unique(cd_data$id_spp)),
        
        # likelihood 1
        n1 = nrow(lik1_df),
        site1 = lik1_df$site,
        true_spp1 = lik1_df$true_spp,
        id_spp1 = lik1_df$id_spp,
        
        # likelihood 2
        n2 = nrow(lik2_df),
        site2 = lik2_df$site,
        id_spp2 = lik2_df$id_spp
      ),
      inits_list = inits_list,
      niter = niter,
      nchains = 1,
      thin = 1
    )
    stopCluster(this_cluster)
    
  } else{
    fit <- fit_model(
      seed = 1:3,
      code = code,
      data = list(
        # likelihood 1
        y1 = lik1_df$count,
        
        # likelihood 2
        y2 = lik2_df$count,
        
        # theta prior
        alpha0 = alpha0
      ),
      constants = list(
        nsites = nrow(site_df),
        nspecies = length(unique(cd_data$id_spp)),
        
        # likelihood 1
        n1 = nrow(lik1_df),
        site1 = lik1_df$site,
        true_spp1 = lik1_df$true_spp,
        id_spp1 = lik1_df$id_spp,
        
        # likelihood 2
        n2 = nrow(lik2_df),
        site2 = lik2_df$site,
        id_spp2 = lik2_df$id_spp
      ),
      inits_list = inits_list,
      niter = niter,
      nchains = 3,
      thin = 1
    )
  }
  
  if(save) saveRDS(fit, file = filename)
  
  out <- fit
  return(out)
}

cd_model_nophi_uncoupledaux <- function(cd_data, save = T, filename = "test.rds", parallel = F, niter = 5000, alpha0, inits_list){
  
  fit_model <- function(seed = 1, code, data, constants, inits_list, niter, nchains, thin = 1){
    library(nimble)
    
    # handle inits in each chain
    inits_ <- list()
    inits_$psi <- inits_list$psi + runif(10, -.05, .05)
    inits_$theta <- inits_list$theta + matrix(runif(100, -.05, .05), 10, 10)
    inits_$theta <- abs(inits_$theta)  / rowSums(abs(inits_$theta))
    inits_$lambda <- inits_list$lambda + runif(10, -inits_list$lambda, inits_list$lambda)
    
    # R model
    model <- nimbleModel(code, constants, data)
    
    # C model
    model_c <- compileNimble(model)
    
    # R mcmc
    model_conf <- configureMCMC(model)
    
    # R mcmc
    mcmc <- buildMCMC(model_conf)
    
    # C mcmc
    mcmc_c <- compileNimble(mcmc, project = model_c)
    
    # run model
    out <- runMCMC(
      mcmc_c, 
      niter = niter, 
      nchains = nchains, 
      thin = thin, 
      init = inits_,
      setSeed = seed
    )
    
    # out
    return(out)
  }
  
  # housekeeping
  nspecies <- length(unique(cd_data$id_spp))
  
  # get aux counts
  aux_counts <- cd_data %>% 
    filter(type == "unambiguous") %>%
    group_by(true_spp, id_spp) %>%
    summarize(total = sum(count)) %>%
    ungroup %>%
    select(total) %>%
    unlist() %>% unname() %>%
    matrix(., nrow = nspecies, ncol = nspecies, byrow = T)
  
  # no covariates
  site_df <- cd_data %>%
    ungroup %>%
    select(
      site
    ) %>%
    distinct()
  
  lik1_df <- cd_data %>%
    filter(type == "unambiguous")
  lik2_df <- cd_data %>%
    filter(type == "ambiguous")
  
  code <- nimbleCode({
    for(k in 1:nspecies){
      # priors
      psi[k] ~ dbeta(1, 1)
      lambda[k] ~ T(dnorm(0, sd = 100), 0, Inf)
      theta[k, 1:nspecies] ~ ddirch(alpha = alpha0[k, 1:nspecies])
      
      # aux counts to inform theta
      aux_counts[k, 1:nspecies] ~ dmulti(size = aux_size[k], prob = theta[k, 1:nspecies])
    }
    
    # likelihood - site level occupancy
    for(site in 1:nsites){
      for(k in 1:nspecies){
        z[site,k] ~ dbern(psi[k])
      }
    }
    
    # likelihood 1
    # for(row in 1:n1){
    #   y1[row] ~ dpois(
    #     z[site1[row],true_spp1[row]] *
    #       lambda[true_spp1[row]] *
    #       theta[true_spp1[row],id_spp1[row]]
    #   )
    # }
    
    # likelihood 2
    for(row in 1:n2){
      for(true in 1:nspecies){
        mean2[row, true] <- z[site2[row],true]*lambda[true]*theta[true,id_spp2[row]]
      }
      y2[row] ~ dpois(sum(mean2[row, 1:nspecies]))
    }
  })
  
  # model housekeeping
  init_func <- function(){
    
    out <- list(
      psi = runif(10),
      lambda = abs(rnorm(10, 0, 100)),
      theta = t(apply(alpha0, 1, function(x) rdirch(1, x))),
      z = matrix(1, nrow = 100, ncol = 10)
    )
    
    return(out)
  }
  
  if(parallel) {
    library(parallel)
    this_cluster <- makeCluster(3)
    fit <- parLapply(
      cl = this_cluster,
      X = 1:3,
      fun = fit_model,
      code = code,
      data = list(
        # # likelihood 1
        # y1 = lik1_df$count,
        
        # aux counts
        aux_counts = aux_counts,
        
        # likelihood 2
        y2 = lik2_df$count,
        
        # theta prior
        alpha0 = alpha0
      ),
      constants = list(
        nsites = nrow(site_df),
        nspecies = length(unique(cd_data$id_spp)),
        
        # # likelihood 1
        # n1 = nrow(lik1_df),
        # site1 = lik1_df$site,
        # true_spp1 = lik1_df$true_spp,
        # id_spp1 = lik1_df$id_spp,
        
        # aux_counts
        aux_size = rowSums(aux_counts),
        
        # likelihood 2
        n2 = nrow(lik2_df),
        site2 = lik2_df$site,
        id_spp2 = lik2_df$id_spp
      ),
      niter = niter,
      inits_list = inits_list,
      nchains = 1,
      thin = 1
    )
    stopCluster(this_cluster)
    
  } else{
    fit <- fit_model(
      seed = 1:3,
      code = code,
      data = list(
        # # likelihood 1
        # y1 = lik1_df$count,
        
        # aux counts
        aux_counts = aux_counts,
        
        # likelihood 2
        y2 = lik2_df$count,
        
        # theta prior
        alpha0 = alpha0
      ),
      constants = list(
        nsites = nrow(site_df),
        nspecies = length(unique(cd_data$id_spp)),
        
        # # likelihood 1
        # n1 = nrow(lik1_df),
        # site1 = lik1_df$site,
        # true_spp1 = lik1_df$true_spp,
        # id_spp1 = lik1_df$id_spp,
        
        # aux_counts
        aux_size = rowSums(aux_counts),
        
        # likelihood 2
        n2 = nrow(lik2_df),
        site2 = lik2_df$site,
        id_spp2 = lik2_df$id_spp
      ),
      inits_list = inits_list,
      niter = niter,
      nchains = 3,
      thin = 1
    )
  }
  
  if(save) saveRDS(fit, file = filename)
  
  out <- fit
  return(out)
}

cd_model_nophi_prioronly <- function(cd_data, save = T, filename = "test.rds", parallel = F, niter = 5000, alpha0, inits_list){
  
  fit_model <- function(seed = 1, code, data, constants, inits_list, niter, nchains, thin = 1, alpha0){
    library(nimble)
    
    # handle inits in each chain
    inits_ <- list()
    inits_$psi <- inits_list$psi + runif(10, -.05, .05)
    inits_$theta <- inits_list$theta + matrix(runif(100, -.05, .05), 10, 10)
    inits_$theta <- abs(inits_$theta)  / rowSums(abs(inits_$theta))
    inits_$lambda <- inits_list$lambda + runif(10, -inits_list$lambda, inits_list$lambda)
    
    # R model
    model <- nimbleModel(code, constants, data)
    
    # C model
    model_c <- compileNimble(model)
    
    # R mcmc
    model_conf <- configureMCMC(model)
    
    # R mcmc
    mcmc <- buildMCMC(model_conf)
    
    # C mcmc
    mcmc_c <- compileNimble(mcmc, project = model_c)
    
    # run model
    out <- runMCMC(
      mcmc_c, 
      niter = niter, 
      nchains = nchains, 
      thin = thin, 
      init = inits_,
      setSeed = seed
    )
    
    # out
    return(out)
  }
  
  # housekeeping
  nspecies <- length(unique(cd_data$id_spp))
  
  # no covariates
  site_df <- cd_data %>%
    ungroup %>%
    select(
      site
    ) %>%
    distinct()
  
  lik1_df <- cd_data %>%
    filter(type == "unambiguous")
  lik2_df <- cd_data %>%
    filter(type == "ambiguous")
  
  code <- nimbleCode({
    # priors
    for(k in 1:nspecies){
      psi[k] ~ dbeta(1, 1)
      lambda[k] ~ T(dnorm(0, sd = 100), 0, Inf)
      theta[k, 1:nspecies] ~ ddirch(alpha = alpha0[k, 1:nspecies])
    }
    
    # likelihood - site level occupancy
    for(site in 1:nsites){
      for(k in 1:nspecies){
        z[site,k] ~ dbern(psi[k])
      }
    }
    
    # # likelihood 1
    # for(row in 1:n1){
    #   y1[row] ~ dpois(
    #     z[site1[row],true_spp1[row]] *
    #       lambda[true_spp1[row]] *
    #       theta[true_spp1[row],id_spp1[row]]
    #   )
    # }
    
    # likelihood 2
    for(row in 1:n2){
      for(true in 1:nspecies){
        mean2[row, true] <- z[site2[row],true]*lambda[true]*theta[true,id_spp2[row]]
      }
      y2[row] ~ dpois(sum(mean2[row, 1:nspecies]))
    }
  })
  
  # # model housekeeping
  # init_func <- function(){
  #   
  #   out <- list(
  #     psi = runif(10),
  #     lambda = abs(rnorm(10, 0, 100)),
  #     theta = t(apply(alpha0, 1, function(x) rdirch(1, x))),
  #     z = matrix(1, nrow = 100, ncol = 10)
  #   )
  #   
  #   return(out)
  # }
  
  if(parallel) {
    library(parallel)
    this_cluster <- makeCluster(3)
    fit <- parLapply(
      cl = this_cluster,
      X = 1:3,
      fun = fit_model,
      code = code,
      data = list(
        # # likelihood 1
        # y1 = lik1_df$count,
        
        # likelihood 2
        y2 = lik2_df$count,
        
        # theta prior
        alpha0 = alpha0
      ),
      constants = list(
        nsites = nrow(site_df),
        nspecies = length(unique(cd_data$id_spp)),
        
        # # likelihood 1
        # n1 = nrow(lik1_df),
        # site1 = lik1_df$site,
        # true_spp1 = lik1_df$true_spp,
        # id_spp1 = lik1_df$id_spp,
        
        # likelihood 2
        n2 = nrow(lik2_df),
        site2 = lik2_df$site,
        id_spp2 = lik2_df$id_spp
      ),
      niter = niter,
      nchains = 1,
      inits_list = inits_list,
      thin = 1,
      alpha0 = alpha0
    )
    stopCluster(this_cluster)
    
  } else{
    fit <- fit_model(
      seed = 1:3,
      code = code,
      data = list(
        # # likelihood 1
        # y1 = lik1_df$count,
        
        # likelihood 2
        y2 = lik2_df$count,
        
        # theta prior
        alpha0 = alpha0
      ),
      constants = list(
        nsites = nrow(site_df),
        nspecies = length(unique(cd_data$id_spp)),
        
        # # likelihood 1
        # n1 = nrow(lik1_df),
        # site1 = lik1_df$site,
        # true_spp1 = lik1_df$true_spp,
        # id_spp1 = lik1_df$id_spp,
        
        # likelihood 2
        n2 = nrow(lik2_df),
        site2 = lik2_df$site,
        id_spp2 = lik2_df$id_spp
      ),
      inits_list = inits_list,
      niter = niter,
      nchains = 3,
      thin = 1,
      alpha0 = alpha0
    )
  }
  
  if(save) saveRDS(fit, file = filename)
  
  out <- fit
  return(out)
}

# summarize fit
nimble_summary <- function(fit, warmup = nrow(fit[[1]])/2, thin = 1){
  # convert to coda for normal summary
  fit_warmup <- lapply(fit, function(x) x[(warmup+1):nrow(x),])
  coda_samples <- as.mcmc.list(lapply(fit_warmup, function(x) as.mcmc(
    x, start = warmup+1, end = nrow(fit), thin = thin
  )))
  
  sum <- summary(coda_samples)
  params <- dimnames(sum$statistics)[[1]]
  tmp_sum <- cbind(sum$statistics, sum$quantiles)
  
  # get r hat / n_eff
  mat <- matrix(NA, nrow = nrow(tmp_sum), ncol = 3)
  colnames(mat) <- c("Rhat", "ess_bulk", "ess_tail")
  for(i in 1:nrow(tmp_sum)){
    tmp <- sapply(fit, function(x) x[,i])
    mat[i,] <- c(Rhat(tmp), ess_bulk(tmp), ess_tail(tmp))
  }
  
  # out 
  out <- cbind(tmp_sum, mat)
  return(out)
}

run_sims <- function(nsims, nsites = 100, nvisits = 8, lambda_true, theta_true, psi_true, niter = 2500){
  message(Sys.time())
  
  # storage
  out <- list()
  
  # run sims
  start <- Sys.time()
  for(sim in 1:nsims){
    # simulate data
    sim_data <- sim_dat(
      nsites = nsites, nspecies = 10, nvisits = nvisits, seed = sim + 100, 
      psi = psi_true, lambda = lambda_true, theta = theta_true,
      phi = 1, ambig_frac = .75
    )
    data <- sim_data$df
    
    # different priors
    alpha0_flat <- matrix(1, 10, 10)
    alpha0_informative <- matrix(1, 10, 10);diag(alpha0_informative) <- 30
    alpha0_reallyflat <- matrix(.1, 10, 10)
    
    # inits
    inits_list_base <- sim_data$params
    
    # model 1 - coupled, uniform prior
    Sys.time()
    model1 <- cd_model_nophi(
      cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_flat,
      save = T, filename = paste0("simulation 1/fits/model1_", sim, ".rds"),
      inits_list = inits_list_base
    )
    Sys.time()
    mod_sum <- nimble_summary(model1)
    sum_df1 <- tibble(
      model = 1,
      sim_num  = sim,
      param = rownames(mod_sum),
      mean = mod_sum[,1],
      lwr = mod_sum[,5],
      upr = mod_sum[,9],
      rhat = mod_sum[,10],
      ess = mod_sum[,11],
      truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    ) %>%
      mutate(
        capture = case_when(
          truth >= lwr & truth <= upr ~ 1,
          TRUE ~ 0
        ),
        width = upr - lwr
      )
    
    # model 2 - uncoupled - aux, uniform prior
    Sys.time()
    model2 <- cd_model_nophi_uncoupledaux(
      cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_flat,
      save = T, filename = paste0("simulation 1/fits/model2_", sim, ".rds"),
      inits_list = inits_list_base
    )
    Sys.time()
    mod_sum <- nimble_summary(model2)
    sum_df2 <- tibble(
      model = 2,
      sim_num  = sim,
      param = rownames(mod_sum),
      mean = mod_sum[,1],
      lwr = mod_sum[,5],
      upr = mod_sum[,9],
      rhat = mod_sum[,10],
      ess = mod_sum[,11],
      truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    ) %>%
      mutate(
        capture = case_when(
          truth >= lwr & truth <= upr ~ 1,
          TRUE ~ 0
        ),
        width = upr - lwr
      )
    
    # # model 3 - uncoupled - prior, uniform prior
    # Sys.time()
    # model3 <- cd_model_nophi_prioronly(
    #   cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_flat,
    #   save = T, filename = paste0("simulation 1/fits/model3_", sim, ".rds"),
    #   init_list = sim_data$params[c(1:3)]
    # )
    # Sys.time()
    # mod_sum <- nimble_summary(model3)
    # sum_df3 <- tibble(
    #   model = 3,
    #   sim_num  = sim,
    #   param = rownames(mod_sum),
    #   mean = mod_sum[,1],
    #   lwr = mod_sum[,5],
    #   upr = mod_sum[,9],
    #   rhat = mod_sum[,10],
    #   ess = mod_sum[,11],
    #   truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    # ) %>%
    #   mutate(
    #     capture = case_when(
    #       truth >= lwr & truth <= upr ~ 1,
    #       TRUE ~ 0
    #     ),
    #     width = upr - lwr
    #   )
    
    # model 4 - coupled, informative prior
    Sys.time()
    model4 <- cd_model_nophi(
      cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_informative,
      save = T, filename = paste0("simulation 1/fits/model4_", sim, ".rds"),
      inits_list = inits_list_base
    )
    Sys.time()
    mod_sum <- nimble_summary(model4)
    sum_df4 <- tibble(
      model = 4,
      sim_num  = sim,
      param = rownames(mod_sum),
      mean = mod_sum[,1],
      lwr = mod_sum[,5],
      upr = mod_sum[,9],
      rhat = mod_sum[,10],
      ess = mod_sum[,11],
      truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    ) %>%
      mutate(
        capture = case_when(
          truth >= lwr & truth <= upr ~ 1,
          TRUE ~ 0
        ),
        width = upr - lwr
      )
    
    # model 5 - uncoupled - aux, informative prior
    Sys.time()
    model5 <- cd_model_nophi_uncoupledaux(
      cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_informative,
      save = T, filename = paste0("simulation 1/fits/model5_", sim, ".rds"),
      inits_list = inits_list_base
    )
    Sys.time()
    mod_sum <- nimble_summary(model5)
    sum_df5 <- tibble(
      model = 5,
      sim_num  = sim,
      param = rownames(mod_sum),
      mean = mod_sum[,1],
      lwr = mod_sum[,5],
      upr = mod_sum[,9],
      rhat = mod_sum[,10],
      ess = mod_sum[,11],
      truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    ) %>%
      mutate(
        capture = case_when(
          truth >= lwr & truth <= upr ~ 1,
          TRUE ~ 0
        ),
        width = upr - lwr
      )
    
    # model 6 - uncoupled - prior, informative prior
    Sys.time()
    model6 <- cd_model_nophi_prioronly(
      cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_informative,
      save = T, filename = paste0("simulation 1/fits/model6_", sim, ".rds"),
      inits_list = inits_list_base
    )
    Sys.time()
    mod_sum <- nimble_summary(model6)
    sum_df6 <- tibble(
      model = 6,
      sim_num  = sim,
      param = rownames(mod_sum),
      mean = mod_sum[,1],
      lwr = mod_sum[,5],
      upr = mod_sum[,9],
      rhat = mod_sum[,10],
      ess = mod_sum[,11],
      truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    ) %>%
      mutate(
        capture = case_when(
          truth >= lwr & truth <= upr ~ 1,
          TRUE ~ 0
        ),
        width = upr - lwr
      )
    
    # model 7 - coupled, really flat prior
    Sys.time()
    model7 <- cd_model_nophi(
      cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_reallyflat,
      save = T, filename = paste0("simulation 1/fits/model7_", sim, ".rds"),
      inits_list = inits_list_base
    )
    Sys.time()
    mod_sum <- nimble_summary(model7)
    sum_df7 <- tibble(
      model = 7,
      sim_num  = sim,
      param = rownames(mod_sum),
      mean = mod_sum[,1],
      lwr = mod_sum[,5],
      upr = mod_sum[,9],
      rhat = mod_sum[,10],
      ess = mod_sum[,11],
      truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    ) %>%
      mutate(
        capture = case_when(
          truth >= lwr & truth <= upr ~ 1,
          TRUE ~ 0
        ),
        width = upr - lwr
      )
    
    # model 8 - uncoupled - aux, really flat prior
    Sys.time()
    model8 <- cd_model_nophi_uncoupledaux(
      cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_reallyflat,
      save = T, filename = paste0("simulation 1/fits/model8_", sim, ".rds"),
      inits_list = inits_list_base
    )
    Sys.time()
    mod_sum <- nimble_summary(model8)
    sum_df8 <- tibble(
      model = 8,
      sim_num  = sim,
      param = rownames(mod_sum),
      mean = mod_sum[,1],
      lwr = mod_sum[,5],
      upr = mod_sum[,9],
      rhat = mod_sum[,10],
      ess = mod_sum[,11],
      truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    ) %>%
      mutate(
        capture = case_when(
          truth >= lwr & truth <= upr ~ 1,
          TRUE ~ 0
        ),
        width = upr - lwr
      )
    
    # # model 9 - uncoupled - prior, really flat prior
    # Sys.time()
    # model9 <- cd_model_nophi_prioronly(
    #   cd_data = data, parallel = T, niter = niter, alpha0 = alpha0_reallyflat,
    #   save = T, filename = paste0("simulation 1/fits/model9_", sim, ".rds"),
    #   init_list = sim_data$params[c(1:3)]
    # )
    # Sys.time()
    # mod_sum <- nimble_summary(model9)
    # sum_df9 <- tibble(
    #   model = 9,
    #   sim_num  = sim,
    #   param = rownames(mod_sum),
    #   mean = mod_sum[,1], 
    #   lwr = mod_sum[,5],
    #   upr = mod_sum[,9],
    #   rhat = mod_sum[,10],
    #   ess = mod_sum[,11],
    #   truth = c(sim_data$params$lambda, sim_data$params$psi, c(sim_data$params$theta))
    # ) %>%
    #   mutate(
    #     capture = case_when(
    #       truth >= lwr & truth <= upr ~ 1, 
    #       TRUE ~ 0
    #     ),
    #     width = upr - lwr
    #   )
    
    # combine
    tmp <- bind_rows(sum_df1, sum_df2, sum_df4, sum_df5, sum_df6, sum_df7, sum_df8)
    # tmp <- sum_df6
    # tmp <- bind_rows(sum_df7, sum_df8, sum_df9)
    saveRDS(tmp, file = "simulation 1/tmp.rds")
    
    out[[sim]] <- tmp
    
    current <- Sys.time()
    dif <- current - start
    message(paste0("Simulation ", sim, " of ", nsims, " complete. Current runtime: ", round(as.numeric(dif), 3), " ", units(dif)))
  }
  
  return(out)
}

# run sims
sim1_new <- run_sims(
  nsims = 100, nsites = 100, nvisits = 8, lambda_true = lambda_true, theta_true = theta_true, psi_true = psi_true, niter = 2500
)
saveRDS(sim1_new, file = "simulation 1/sim1_new.rds")

# summarize sims
# read in file
sim1 <- readRDS("simulation 1/sim1_new.rds") %>%
  data.table::rbindlist(.) %>%
  as_tibble

# add coverage
plot_tbl <- sim1 %>%
  group_by(model, param) %>%
  mutate(coverage = mean(capture)) %>%
  mutate(model = factor(model)) 

# convergence issues
plot_tbl %>%
  group_by(sim_num, model) %>%
  summarize(max_rhat = max(rhat)) %>%
  filter(max_rhat >= 1.2) %>%
  ungroup %>%
  dplyr::select(model) %>% table

plot_tbl <- anti_join(
  plot_tbl, 
  plot_tbl %>%
    group_by(sim_num, model) %>%
    summarize(max_rhat = max(rhat)) %>%
    filter(max_rhat >= 1.2),
  by = c("sim_num", "model")
)

load("simulation 1/parameter_estimates.rdata")

# lambda
params <- tibble(
  species_char = names(lambda_true),
  lambda = lambda_true
)
## credibility interval width - all species
{
  p <- plot_tbl %>%
    filter(grepl("lambda", param)) %>%
    mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, mean_width, max, min, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 8, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    mutate(
      species_char = case_when(
        species == 1 ~ "EPFU",
        species == 2 ~ "LACI",
        species == 3 ~ "LANO",
        species == 4 ~ "MYCA",
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU",
        species == 8 ~ "MYVO",
        species == 9 ~ "MYYU",
        species == 10 ~ "other"
      )
    ) %>%
    distinct %>%
    ggplot() +
    geom_pointrange(aes(x = model_, y = mean_width, ymin = min, ymax = max, col = coverage)) +
    facet_wrap(~ species_char) +
    theme_bw() +
    scale_color_gradient(low = "red", high = "blue", limit = c(.5,1))+
    theme(legend.position = "bottom") + 
    # theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Relative activity credibility interval width", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_lambda.png", p, width = 2800, units = "px")
  
}

## credibility interval width - subset species
{
  p <- plot_tbl %>%
    filter(grepl("lambda", param)) %>%
    mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, mean_width, max, min, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 8, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    filter(species %in% c(5:7, 10)) %>%
    mutate(
      species_char = case_when(
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU",
        species == 10 ~ "other"
      )
    ) %>%
    distinct %>%
    ggplot() +
    geom_pointrange(aes(x = model_, y = mean_width, ymin = min, ymax = max, col = coverage)) +
    facet_wrap(~ species_char) +
    theme_bw() +
    scale_color_gradient(low = "red", high = "blue", limit = c(.5,1))+
    theme(legend.position = "bottom") + 
    # theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Relative activity credibility interval width", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_lambda_sub.png", p, width = 2800, units = "px")
}

# raw - subset species
{
  plot_tbl2 <- plot_tbl %>%
    filter(grepl("lambda", param)) %>%
    # mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 8, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    filter(species %in% c(2, 5:7)) %>%
    mutate(
      species_char = case_when(
        species == 2 ~ "LACI",
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU"
      )
    ) %>%
    distinct
  
  plot_tbl3 <- plot_tbl %>%
    filter(grepl("lambda", param)) %>%
    # mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 8, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    filter(species %in% c(2, 5:7)) %>%
    mutate(
      species_char = case_when(
        species == 2 ~ "LACI",
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU"
      )
    ) %>%
    distinct %>%
    group_by(model_, species_char) %>%
    summarize(mean = mean(mean), lwr = mean(lwr), upr = mean(upr), coverage = coverage)
  
  p <- ggplot() +
    geom_linerange(
      data = plot_tbl2,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage),
      position = position_jitter(width = .35),
      alpha = .10
    ) +
    geom_pointrange(
      data = plot_tbl3,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage)
    ) +
    facet_wrap(~ species_char, scales = "free_y") +
    geom_hline(
      data = params %>% 
        filter(species_char %in% c("LACI", "MYCI", "MYEV", "MYLU")), 
      aes(yintercept = lambda),
      linetype = "dotdash"
    ) + 
    theme_bw() + 
    # scale_x_discrete(expand = expansion(mult = 1.1)) +
    scale_color_gradient(low = "red", high = "blue", limit = c(.5,1))+
    theme(legend.position = "bottom") + 
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Posterior mean relative activity", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_lambda_sub_raw.png", p, width = 2800, units = "px")
}

# raw - all species
{
  plot_tbl2 <- plot_tbl %>%
    filter(grepl("lambda", param)) %>%
    # mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 8, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    mutate(
      species_char = case_when(
        species == 1 ~ "EPFU",
        species == 2 ~ "LACI",
        species == 3 ~ "LANO",
        species == 4 ~ "MYCA",
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU",
        species == 8 ~ "MYVO",
        species == 9 ~ "MYYU",
        species == 10 ~ "other"
      )
    ) %>%
    distinct
  
  plot_tbl3 <- plot_tbl %>%
    filter(grepl("lambda", param)) %>%
    # mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 8, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    mutate(
      species_char = case_when(
        species == 1 ~ "EPFU",
        species == 2 ~ "LACI",
        species == 3 ~ "LANO",
        species == 4 ~ "MYCA",
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU",
        species == 8 ~ "MYVO",
        species == 9 ~ "MYYU",
        species == 10 ~ "other"
      )
    ) %>%
    distinct %>%
    group_by(model_, species_char) %>%
    summarize(mean = mean(mean), lwr = mean(lwr), upr = mean(upr), coverage = coverage)
  
  p <- ggplot() +
    geom_linerange(
      data = plot_tbl2,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage),
      position = position_jitter(width = .35),
      alpha = .10
    ) +
    geom_pointrange(
      data = plot_tbl3,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage)
    ) +
    geom_hline(
      data = params,
      aes(yintercept = lambda),
      linetype = "dotdash"
    ) + 
    facet_wrap(~ species_char, scales = "free_y") +
    theme_bw() +
    # scale_x_discrete(expand = expansion(mult = 1.1)) +
    scale_color_gradient(low = "red", high = "blue", limit = c(.5,1))+
    theme(legend.position = "bottom") + 
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Posterior mean relative activity", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_lambda_raw.png", p, width = 2800, units = "px")
}

# psi
params <- tibble(
  species_char = names(lambda_true),
  psi = psi_true
)
## cred interval width - all
{
  p <- plot_tbl %>%
    filter(grepl("psi", param)) %>%
    mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, mean_width, max, min, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 5, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    mutate(
      species_char = case_when(
        species == 1 ~ "EPFU",
        species == 2 ~ "LACI",
        species == 3 ~ "LANO",
        species == 4 ~ "MYCA",
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU",
        species == 8 ~ "MYVO",
        species == 9 ~ "MYYU",
        species == 10 ~ "other"
      )
    ) %>%
    distinct %>%
    ggplot() +
    geom_pointrange(aes(x = model_, y = mean_width, ymin = min, ymax = max, col = coverage)) +
    facet_wrap(~ species_char) +
    theme_bw() +
    scale_color_gradient(low = "red", high = "blue", limit = c(.5,1))+
    theme(legend.position = "bottom") + 
    # theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Occupancy probability credibility interval width", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_psi.png", p, width = 2800, units = "px")
}

# raw - all species
{
  plot_tbl2 <- plot_tbl %>%
    filter(grepl("psi", param)) %>%
    # mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 5, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    mutate(
      species_char = case_when(
        species == 1 ~ "EPFU",
        species == 2 ~ "LACI",
        species == 3 ~ "LANO",
        species == 4 ~ "MYCA",
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU",
        species == 8 ~ "MYVO",
        species == 9 ~ "MYYU",
        species == 10 ~ "other"
      )
    ) %>%
    distinct
  
  plot_tbl3 <- plot_tbl %>%
    filter(grepl("psi", param)) %>%
    # mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species = sapply(param, function(x) stringr::str_sub(x, start = 5, end = nchar(x)-1)) %>% as.numeric
    ) %>%
    mutate(
      species_char = case_when(
        species == 1 ~ "EPFU",
        species == 2 ~ "LACI",
        species == 3 ~ "LANO",
        species == 4 ~ "MYCA",
        species == 5 ~ "MYCI",
        species == 6 ~ "MYEV",
        species == 7 ~ "MYLU",
        species == 8 ~ "MYVO",
        species == 9 ~ "MYYU",
        species == 10 ~ "other"
      )
    ) %>%
    distinct %>%
    group_by(model_, species_char) %>%
    summarize(mean = mean(mean), lwr = mean(lwr), upr = mean(upr), coverage = coverage)
  
  p <- ggplot() +
    geom_linerange(
      data = plot_tbl2,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage),
      position = position_jitter(width = .35),
      alpha = .10
    ) +
    geom_pointrange(
      data = plot_tbl3,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage)
    ) +
    geom_hline(
      data = params,
      aes(yintercept = psi),
      linetype = "dotdash"
    ) + 
    facet_wrap(~ species_char, scales = "fixed") +
    theme_bw() +
    # scale_x_discrete(expand = expansion(mult = 1.1)) +
    scale_color_gradient(low = "red", high = "blue", limit = c(.5,1))+
    theme(legend.position = "bottom") + 
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Posterior mean relative activity", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_psi_raw.png", p, width = 2800, units = "px")
}

# theta
params <- tibble(
  species_true_char = rep(names(psi_true), length(psi_true)),
  species_id_char = rep(names(psi_true), each = length(psi_true)),
  theta = c(theta_true)
)

# theta
## cred interval - sub 
{
  p <- plot_tbl %>%
    filter(grepl("theta", param)) %>%
    mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, mean_width, max, min, coverage) %>%
    distinct %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species_true = sapply(
        param, function(x) stringr::str_sub(
          x, start = 7, apply(str_locate(x, pattern = ","), 1, unique)-1
        )
      ) %>% as.numeric,
      species_id = sapply(
        param, function(x) stringr::str_sub(
          x, start = apply(str_locate(x, pattern = ","), 1, unique)+2, nchar(x)-1
        )
      ) %>% as.numeric
    ) %>%
    mutate(
      species_true_char = case_when(
        species_true == 1 ~ "EPFU",
        species_true == 2 ~ "LACI",
        species_true == 3 ~ "LANO",
        species_true == 4 ~ "MYCA",
        species_true == 5 ~ "MYCI",
        species_true == 6 ~ "MYEV",
        species_true == 7 ~ "MYLU",
        species_true == 8 ~ "MYVO",
        species_true == 9 ~ "MYYU",
        species_true == 10 ~ "other"
      )
    ) %>%
    mutate(
      species_id_char = case_when(
        species_id == 1 ~ "EPFU",
        species_id == 2 ~ "LACI",
        species_id == 3 ~ "LANO",
        species_id == 4 ~ "MYCA",
        species_id == 5 ~ "MYCI",
        species_id == 6 ~ "MYEV",
        species_id == 7 ~ "MYLU",
        species_id == 8 ~ "MYVO",
        species_id == 9 ~ "MYYU",
        species_id == 10 ~ "other"
      )
    ) %>%
    left_join(
      tibble(
        species_true = rep(c(2,5,6,7), each = 4),
        species_id = rep(c(2,5,6,7), 4)
      ),
      ., 
      by = c("species_true", "species_id")
    ) %>%
    ggplot() +
    geom_pointrange(aes(x = model_, y = mean_width, ymin = min, ymax = max, col = coverage)) +
    facet_grid(species_true_char ~ species_id_char) +
    theme_bw() +
    scale_color_gradient(low = "red", high = "blue", limit = c(0,1))+
    theme(legend.position = "bottom") + 
    # theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Classification probability credibility interval width", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_theta_sub.png", p, width = 2800, units = "px")
}

## cred interval - all
{
  p <- plot_tbl %>%
    filter(grepl("theta", param)) %>%
    mutate(min = min(width), max = max(width), mean_width = mean(width)) %>%
    select(model, param, mean_width, max, min, coverage) %>%
    distinct %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species_true = sapply(
        param, function(x) stringr::str_sub(
          x, start = 7, apply(str_locate(x, pattern = ","), 1, unique)-1
        )
      ) %>% as.numeric,
      species_id = sapply(
        param, function(x) stringr::str_sub(
          x, start = apply(str_locate(x, pattern = ","), 1, unique)+2, nchar(x)-1
        )
      ) %>% as.numeric
    ) %>%
    mutate(
      species_true_char = case_when(
        species_true == 1 ~ "EPFU",
        species_true == 2 ~ "LACI",
        species_true == 3 ~ "LANO",
        species_true == 4 ~ "MYCA",
        species_true == 5 ~ "MYCI",
        species_true == 6 ~ "MYEV",
        species_true == 7 ~ "MYLU",
        species_true == 8 ~ "MYVO",
        species_true == 9 ~ "MYYU",
        species_true == 10 ~ "other"
      )
    ) %>%
    mutate(
      species_id_char = case_when(
        species_id == 1 ~ "EPFU",
        species_id == 2 ~ "LACI",
        species_id == 3 ~ "LANO",
        species_id == 4 ~ "MYCA",
        species_id == 5 ~ "MYCI",
        species_id == 6 ~ "MYEV",
        species_id == 7 ~ "MYLU",
        species_id == 8 ~ "MYVO",
        species_id == 9 ~ "MYYU",
        species_id == 10 ~ "other"
      )
    ) %>%
    ggplot() +
    geom_pointrange(aes(x = model_, y = mean_width, ymin = min, ymax = max, col = coverage), size = .15) +
    facet_grid(species_true_char ~ species_id_char) +
    theme_bw() +
    scale_color_gradient(low = "red", high = "blue", limit = c(0,1))+
    theme(legend.position = "bottom") + 
    # theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Classification probability credibility interval width", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_theta.png", p, width = 2800, units = "px")
}

## raw - sub
{
  plot_tbl2 <- plot_tbl %>%
    filter(grepl("theta", param)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    distinct %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species_true = sapply(
        param, function(x) stringr::str_sub(
          x, start = 7, apply(str_locate(x, pattern = ","), 1, unique)-1
        )
      ) %>% as.numeric,
      species_id = sapply(
        param, function(x) stringr::str_sub(
          x, start = apply(str_locate(x, pattern = ","), 1, unique)+2, nchar(x)-1
        )
      ) %>% as.numeric
    ) %>%
    mutate(
      species_true_char = case_when(
        species_true == 1 ~ "EPFU",
        species_true == 2 ~ "LACI",
        species_true == 3 ~ "LANO",
        species_true == 4 ~ "MYCA",
        species_true == 5 ~ "MYCI",
        species_true == 6 ~ "MYEV",
        species_true == 7 ~ "MYLU",
        species_true == 8 ~ "MYVO",
        species_true == 9 ~ "MYYU",
        species_true == 10 ~ "other"
      )
    ) %>%
    mutate(
      species_id_char = case_when(
        species_id == 1 ~ "EPFU",
        species_id == 2 ~ "LACI",
        species_id == 3 ~ "LANO",
        species_id == 4 ~ "MYCA",
        species_id == 5 ~ "MYCI",
        species_id == 6 ~ "MYEV",
        species_id == 7 ~ "MYLU",
        species_id == 8 ~ "MYVO",
        species_id == 9 ~ "MYYU",
        species_id == 10 ~ "other"
      )
    ) %>%
    left_join(
      tibble(
        species_true = rep(c(2,5,6,7), each = 4),
        species_id = rep(c(2,5,6,7), 4)
      ),
      ., 
      by = c("species_true", "species_id")
    ) 
  
  plot_tbl3 <- plot_tbl %>%
    filter(grepl("theta", param)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    distinct %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species_true = sapply(
        param, function(x) stringr::str_sub(
          x, start = 7, apply(str_locate(x, pattern = ","), 1, unique)-1
        )
      ) %>% as.numeric,
      species_id = sapply(
        param, function(x) stringr::str_sub(
          x, start = apply(str_locate(x, pattern = ","), 1, unique)+2, nchar(x)-1
        )
      ) %>% as.numeric
    ) %>%
    mutate(
      species_true_char = case_when(
        species_true == 1 ~ "EPFU",
        species_true == 2 ~ "LACI",
        species_true == 3 ~ "LANO",
        species_true == 4 ~ "MYCA",
        species_true == 5 ~ "MYCI",
        species_true == 6 ~ "MYEV",
        species_true == 7 ~ "MYLU",
        species_true == 8 ~ "MYVO",
        species_true == 9 ~ "MYYU",
        species_true == 10 ~ "other"
      )
    ) %>%
    mutate(
      species_id_char = case_when(
        species_id == 1 ~ "EPFU",
        species_id == 2 ~ "LACI",
        species_id == 3 ~ "LANO",
        species_id == 4 ~ "MYCA",
        species_id == 5 ~ "MYCI",
        species_id == 6 ~ "MYEV",
        species_id == 7 ~ "MYLU",
        species_id == 8 ~ "MYVO",
        species_id == 9 ~ "MYYU",
        species_id == 10 ~ "other"
      )
    ) %>%
    left_join(
      tibble(
        species_true_char = rep(c("LACI", "MYCI", "MYEV", "MYLU"), each = 4),
        species_id_char = rep(c("LACI", "MYCI", "MYEV", "MYLU"), 4)
      ),
      ., 
      by = c("species_true_char", "species_id_char")
    )  %>%
    group_by(model_, species_true_char, species_id_char) %>%
    summarize(mean = mean(mean), lwr = mean(lwr), upr = mean(upr), coverage = coverage) %>%
    distinct
  
  p <- ggplot() +
    geom_linerange(
      data = plot_tbl2,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage),
      position = position_jitter(width = .35),
      alpha = .05
    ) +
    geom_pointrange(
      data = plot_tbl3 %>% ungroup,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage),
      size = .15
    ) +
    geom_hline(
      data = left_join(
        tibble(
          species_true_char = rep(c("LACI", "MYCI", "MYEV", "MYLU"), each = 4),
          species_id_char = rep(c("LACI", "MYCI", "MYEV", "MYLU"), 4)
        ),
        params, 
        by = c("species_true_char", "species_id_char")
      ),
      aes(yintercept = theta),
      linetype = "dotdash"
    ) + 
    facet_grid(species_true_char ~ species_id_char, scale = "free") +
    theme_bw() +
    # scale_x_discrete(expand = expansion(mult = 1.1)) +
    scale_color_gradient(low = "red", high = "blue", limit = c(0,1))+
    theme(legend.position = "bottom") + 
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Posterior mean classification probabilities", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_theta_sub_raw.png", p, width = 2800, units = "px")
}

## raw - all
{
  plot_tbl2 <- plot_tbl %>%
    filter(grepl("theta", param)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    distinct %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species_true = sapply(
        param, function(x) stringr::str_sub(
          x, start = 7, apply(str_locate(x, pattern = ","), 1, unique)-1
        )
      ) %>% as.numeric,
      species_id = sapply(
        param, function(x) stringr::str_sub(
          x, start = apply(str_locate(x, pattern = ","), 1, unique)+2, nchar(x)-1
        )
      ) %>% as.numeric
    ) %>%
    mutate(
      species_true_char = case_when(
        species_true == 1 ~ "EPFU",
        species_true == 2 ~ "LACI",
        species_true == 3 ~ "LANO",
        species_true == 4 ~ "MYCA",
        species_true == 5 ~ "MYCI",
        species_true == 6 ~ "MYEV",
        species_true == 7 ~ "MYLU",
        species_true == 8 ~ "MYVO",
        species_true == 9 ~ "MYYU",
        species_true == 10 ~ "other"
      )
    ) %>%
    mutate(
      species_id_char = case_when(
        species_id == 1 ~ "EPFU",
        species_id == 2 ~ "LACI",
        species_id == 3 ~ "LANO",
        species_id == 4 ~ "MYCA",
        species_id == 5 ~ "MYCI",
        species_id == 6 ~ "MYEV",
        species_id == 7 ~ "MYLU",
        species_id == 8 ~ "MYVO",
        species_id == 9 ~ "MYYU",
        species_id == 10 ~ "other"
      )
    ) 
  
  plot_tbl3 <- plot_tbl %>%
    filter(grepl("theta", param)) %>%
    select(model, param, sim_num, mean, lwr, upr, coverage) %>%
    mutate(model = as.numeric(as.character(model))) %>%
    distinct %>%
    mutate(model = as.numeric(as.character(model))) %>%
    mutate(
      model_ = case_when(
        model == 1 ~ "Cou - Unif",
        model == 2 ~ "Uncou_Aux - Unif",
        model == 3 ~ "Uncou_Prior - Unif",
        model == 4 ~ "Cou - Inf",
        model == 5 ~ "Uncou_Aux - Inf",
        model == 6 ~ "Uncou_Prior - Inf",
        model == 7 ~ "Cou - RD",
        model == 8 ~ "Uncou_Aux - RD",
        TRUE ~ "Uncou_Prior - RD"
      )
    ) %>%
    mutate(
      species_true = sapply(
        param, function(x) stringr::str_sub(
          x, start = 7, apply(str_locate(x, pattern = ","), 1, unique)-1
        )
      ) %>% as.numeric,
      species_id = sapply(
        param, function(x) stringr::str_sub(
          x, start = apply(str_locate(x, pattern = ","), 1, unique)+2, nchar(x)-1
        )
      ) %>% as.numeric
    ) %>%
    mutate(
      species_true_char = case_when(
        species_true == 1 ~ "EPFU",
        species_true == 2 ~ "LACI",
        species_true == 3 ~ "LANO",
        species_true == 4 ~ "MYCA",
        species_true == 5 ~ "MYCI",
        species_true == 6 ~ "MYEV",
        species_true == 7 ~ "MYLU",
        species_true == 8 ~ "MYVO",
        species_true == 9 ~ "MYYU",
        species_true == 10 ~ "other"
      )
    ) %>%
    mutate(
      species_id_char = case_when(
        species_id == 1 ~ "EPFU",
        species_id == 2 ~ "LACI",
        species_id == 3 ~ "LANO",
        species_id == 4 ~ "MYCA",
        species_id == 5 ~ "MYCI",
        species_id == 6 ~ "MYEV",
        species_id == 7 ~ "MYLU",
        species_id == 8 ~ "MYVO",
        species_id == 9 ~ "MYYU",
        species_id == 10 ~ "other"
      )
    ) %>%
    group_by(model_, species_true_char, species_id_char) %>%
    summarize(mean = mean(mean), lwr = mean(lwr), upr = mean(upr), coverage = coverage) %>%
    distinct
  
  p <- ggplot() +
    geom_linerange(
      data = plot_tbl2,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage),
      position = position_jitter(width = .35),
      alpha = .05
    ) +
    geom_pointrange(
      data = plot_tbl3 %>% ungroup,
      aes(x = model_, y = mean, ymin = lwr, ymax = upr, col = coverage),
      size = .15
    ) +
    geom_hline(
      data = params,
      aes(yintercept = theta),
      linetype = "dotdash"
    ) + 
    facet_grid(species_true_char ~ species_id_char, scale = "free") +
    theme_bw() +
    # scale_x_discrete(expand = expansion(mult = 1.1)) +
    scale_color_gradient(low = "red", high = "blue", limit = c(0,1))+
    theme(legend.position = "bottom") + 
    theme(axis.text.x = element_text(angle=65,hjust=1)) +
    labs(y = "Posterior mean classification probabilities", x = "Model")
  
  p
  ggsave("simulation 1/plots/sim1_theta_raw.png", p, width = 2800, units = "px")
}





