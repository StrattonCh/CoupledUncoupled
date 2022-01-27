# clear workspace
rm(list = ls())

# packages
packs <- c("dplyr", "nimble", "htmltools", "ggplot2", "sf", "Rcpp", "RcppArmadillo", "inline", "mvtnorm", "readr", "parallel", "xtable", "rstan", "coda", "vegan", "tidyverse", "lubridate", "bcmaps", "knitr", "kableExtra", "ggalluvial")
sapply(packs, require, character.only = T)
rm(packs)

# convenience
`%notin%` <- Negate("%in%")
options(mc.cores = parallel::detectCores())

# load data
bc_final <- readRDS("bc_final.rds") %>% 
  arrange(grts, sample_night, quad)

# functions
cd_format <- function(data, seed = NULL, propsites_unconf = 0, propvisits_unconf = .5, remove_other = TRUE){
  # function to format data for CD model
  # each row represents a single call
  # not super generalized right now - can work on 
  
  sample_down <- function(.data, frac) {
    sample_n(.data, floor({{frac}} * n()) )
  }
  
  # settings
  options(dplyr.summarise.inform = FALSE)
  if(!is.null(seed)) set.seed(seed)
  
  # response data grab first and last visit
  tmp <- data %>%
    select(grts, year, quad, sample_night) %>%
    distinct() %>%
    group_by(grts, quad, year) %>%
    arrange(grts, year, sample_night, quad) %>%
    mutate(visit = 1:n()) %>%
    filter(visit %in% c(min(visit), max(visit))) %>%
    select(-visit) %>%
    left_join(
      .,
      data %>% select(grts, quad, man, kal, sample_night, year)
    )
  
  # lump uncommon species
  tmp2 <- tmp %>% 
    ungroup() %>% 
    select(kal) %>%
    group_by(kal) %>% 
    summarize(count = n()) %>%
    filter(count < 1000)
  tmp3 <- tmp %>%
    mutate(
      man = ifelse(man %in% unlist(tmp2[,1]), "other", man),
      kal = ifelse(kal %in% unlist(tmp2[,1]), "other", kal)
    )
  
  if(remove_other){
    tmp3 <- tmp3 %>%
      filter(
        man != "other", kal != "other" 
      )
    
  }
  
  # insert NAs for manual call for some visits - unique to BC data
  ## grab sites that are completely unvetted first
  
  if(propsites_unconf == 0) {
    tmp4 <- tmp3 %>%
      ungroup %>%
      select(grts, year, quad, sample_night) %>%
      distinct %>%
      group_by(grts, year) %>%
      sample_down(propvisits_unconf) %>%
      mutate(man = NA) %>%
      ungroup()
    
    
  } else{
    # sites left completely ambiguous
    tmp3point1 <- tmp3 %>%
      ungroup %>%
      select(grts, year) %>%
      distinct %>%
      sample_down(propsites_unconf) %>%
      left_join(
        .,
        tmp3 %>% select(grts, year, quad, sample_night),
        by = c("grts", "year")
      ) %>%
      distinct() %>%
      mutate(man = NA)  %>%
      arrange(grts, year, quad, sample_night)
    
    # some visits from remaining sites
    tmp3point2 <- tmp3 %>%
      ungroup %>%
      anti_join(., tmp3point1, by = c("grts", "year")) %>%
      select(grts, year, quad, sample_night) %>%
      distinct %>%
      group_by(grts, year) %>%
      sample_down(propvisits_unconf) %>%
      mutate(man = NA) %>%
      ungroup()
    
    tmp4 <- bind_rows(
      tmp3point1, tmp3point2
    ) %>%
      arrange(grts, year, sample_night, quad)
  }
  
  tmp5 <- bind_rows(
    anti_join(tmp3, tmp4, by = c("grts", "year", "quad", "sample_night")),
    inner_join(
      tmp3 %>% select(-man),
      tmp4,
      by = c("grts", "year", "quad", "sample_night")
    )
  ) %>%
    arrange(grts, year, sample_night, quad)
  
  # if(propsites_unconf == 0) {
  #   tmp4 <- tmp3 %>%
  #     ungroup %>%
  #     select(grts, year, quad, sample_night) %>%
  #     distinct %>%
  #     group_by(grts, year) %>%
  #     sample_down(propvisits_unconf) %>%
  #     mutate(man = NA) %>%
  #     ungroup()
  #   
  #   
  # } else{
  #   # sites left completely ambiguous
  #   tmp3point1 <- tmp3 %>%
  #     ungroup %>%
  #     select(grts, year) %>% 
  #     distinct %>%
  #     sample_down(propsites_unconf) %>%
  #     left_join(
  #       ., 
  #       tmp3 %>% select(grts, year, quad, sample_night), 
  #       by = c("grts", "year")
  #     ) %>%
  #     distinct() %>%
  #     mutate(man = NA)  %>%
  #     arrange(grts, year, quad, sample_night)
  #   
  #   # some visits from remaining sites
  #   tmp3point2 <- tmp3 %>%
  #     ungroup %>%
  #     anti_join(., tmp3point1, by = c("grts", "year")) %>%
  #     select(grts, year, quad, sample_night) %>%
  #     distinct %>%
  #     group_by(grts, year) %>%
  #     sample_down(propvisits_unconf) %>%
  #     mutate(man = NA) %>%
  #     ungroup()
  #   
  #   tmp4 <- bind_rows(
  #     tmp3point1, tmp3point2
  #   ) %>%
  #     arrange(grts, year, sample_night, quad)
  # }
  
  tmp5 <- bind_rows(
    anti_join(tmp3, tmp4, by = c("grts", "year", "quad", "sample_night")),
    inner_join(
      tmp3 %>% select(-man),
      tmp4,
      by = c("grts", "year", "quad", "sample_night")
    )
  ) %>%
    arrange(grts, year, sample_night, quad)
  
  # create grts cell visit number
  tmp6 <- tmp5 %>%
    select(grts, quad, sample_night, year) %>%
    distinct() %>%
    group_by(grts, year) %>%
    mutate(visit = 1:n()) %>%
    left_join(
      .,
      tmp5 %>% select(grts, quad, man, kal, sample_night, year)
    ) %>%
    ungroup() %>%
    mutate(
      man = factor(man),
      kal = factor(kal)
    )
  
  # summarize counts by visit
  tmp7 <- tmp6 %>%
    group_by(grts, year, visit) %>%
    group_split()
  
  tmp8 <- list()
  key <- expand.grid(levels(tmp6$man), levels(tmp6$kal)) %>%
    as_tibble() %>%
    rename(man = Var1, kal = Var2) %>%
    mutate(count = 0)
  for(i in 1:length(tmp7)){
    df <- tmp7[[i]]
    if(all(is.na(df$man))){
      df2 <- df %>%
        group_by(kal)%>%
        summarize(
          count = n()
        ) 
      tmp8[[i]] <- bind_rows(
        df2, 
        anti_join(key %>% select(-man) %>% distinct, df2, by = "kal")
      ) %>%
        mutate(visit = df$visit[1]) %>%
        left_join(
          ., 
          df %>% select(-kal, -man) %>% distinct(), 
          by = "visit"
        ) %>%
        mutate(true_spp = NA) %>%
        rename(id_spp = kal) %>%
        select(grts, year, quad, visit, sample_night, id_spp, true_spp, count) %>%
        mutate(type = "ambiguous")
      
    } else{
      df2 <- df %>%
        group_by(kal, man)%>%
        summarize(
          count = n()
        )
      tmp8[[i]] <- bind_rows(
        df2, 
        anti_join(key, df2, by = c("kal", "man"))
      ) %>%
        mutate(visit = df$visit[1]) %>%
        left_join(
          ., 
          df %>% select(-kal, -man) %>% distinct(), 
          by = "visit"
        ) %>%
        rename(id_spp = kal, true_spp = man) %>%
        select(grts, year, quad, visit, sample_night, id_spp, true_spp, count) %>%
        mutate(type = "unambiguous")
    }
  }
  tmp9 <- do.call("rbind", tmp8)
  
  # add back in covariates
  out <- left_join(
    tmp9, 
    data %>% select(-man, -kal, -sb, -man_dirty) %>% distinct(),
    by = c("grts", "quad", "sample_night", "year")
  ) %>%
    ungroup()
  
  return(out)
}

cd_model_nophi_covs <- function(cd_data, save = T, filename = "test.rds", parallel = F, niter = 5000, nchains = 3, nburnin, hierarchical = F, site_model, detection_model){
  
  fit_model <- function(seed = 1, code, data, constants, inits, niter, nchains, thin = 1, nburnin = 0, hierarchical = F){
    library(nimble)
    
    # R model
    model <- nimbleModel(code, constants, data)
    
    # C model
    model_c <- compileNimble(model)
    
    # R mcmc
    model_conf <- configureMCMC(model)
    if(hierarchical) model_conf$addMonitors(c("beta", "alpha"))
    
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
      init = inits,
      setSeed = seed,
      nburnin = nburnin
    )
    
    # out
    return(out)
  }
  
  # housekeeping
  nspecies <- length(unique(cd_data$id_spp))
  
  # site covariates
  site_form <- as.formula(site_model)
  site_covs <- all.vars(site_form)
  site_covs_ <- ifelse(grepl("year", site_covs), site_covs[-which(site_covs == "year")], site_covs)
  site_df <- cd_data %>%
    ungroup %>%
    mutate(year = factor(year)) %>%
    select(
      site, site_covs
    ) %>%
    distinct() %>%
    mutate_at(.vars = site_covs_, .funs = ~(scale(.) %>% as.vector))
  X <- model.matrix(as.formula(site_model), data = site_df) %>% as.matrix
  
  # activity covariates
  act_form <- as.formula(detection_model)
  act_covs <<- all.vars(act_form)
  lik1_df <- cd_data %>%
    ungroup %>%
    filter(type == "unambiguous") %>%
    mutate_at(vars(act_covs), ~(scale(.) %>% as.vector))
  W1 <- model.matrix(as.formula(detection_model), data = lik1_df)
  lik2_df <- cd_data %>%
    ungroup %>%
    filter(type == "ambiguous") %>%
    mutate_at(vars(act_covs), ~(scale(.) %>% as.vector))
  W2 <- model.matrix(as.formula(detection_model), data = lik2_df)
  
  # initialize theta
  theta_inits_count <- cd_data %>%
    ungroup %>%
    filter(type == "unambiguous") %>%
    select(true_spp, id_spp, count) %>%
    group_by(true_spp, id_spp) %>%
    summarize(total = sum(count)) %>%
    ungroup %>%
    select(total) %>%
    unlist %>%
    unname %>%
    matrix(., nrow = nspecies, ncol = nspecies, byrow = T)
  theta_inits <- (theta_inits_count + 1) / rowSums(theta_inits_count + 1)
  
  # initialize beta
  initial_site_df <- cd_data %>% 
    group_by(site, id_spp) %>%
    summarize(
      total_count = sum(count)
    ) %>%
    mutate(naive_occ = ifelse(total_count == 0, 0, 1)) %>%
    ungroup %>%
    left_join(
      ., 
      cd_data %>% ungroup %>% select(site, site_covs),
      by = "site"
    ) %>%
    distinct() %>%
    mutate_at(vars(site_covs), ~(scale(.) %>% as.vector))
  
  beta_inits <- matrix(NA, nrow = ncol(X), ncol = nspecies)
  for(i in 1:nspecies){
    df <- initial_site_df %>% filter(id_spp == i)
    beta_inits[,i] <- coef(glm(
      as.formula(paste0("naive_occ ~ ", paste0(site_covs, collapse = " + "))), 
      family = "binomial", 
      data = df %>% mutate(year = factor(year))
    ))
  }
  
  # initialize alpha
  initial_act_df <- cd_data %>% 
    group_by(site, visit, id_spp) %>%
    summarize(
      total_count = sum(count)
    ) %>%
    mutate(naive_occ = ifelse(total_count == 0, 0, 1)) %>%
    ungroup %>%
    left_join(
      ., 
      cd_data %>% ungroup %>% select(site, visit, act_covs),
      by = c("site", "visit")
    ) %>%
    distinct() %>%
    mutate_at(vars(act_covs), ~(scale(.) %>% as.vector)) %>%
    filter(naive_occ != 0)
  
  alpha_inits <- matrix(NA, nrow = ncol(W1), ncol = nspecies)
  for(i in 1:nspecies){
    df <- initial_act_df %>% filter(id_spp == i)
    alpha_inits[,i] <- coef(glm(
      as.formula(paste0("total_count ~ ", paste0(act_covs, collapse = " + "))), 
      family = "poisson", 
      data = df
    ))
  }
  
  if(hierarchical){
    code <- nimbleCode({
      # priors
      for(p in 1:p_beta){
        mu_beta[p] ~ dnorm(0, sd = 2)
        sigma2_beta[p] ~ T(dnorm(0, sd = 2), 0, Inf)
      }
      for(p in 1:p_beta){
        for(k in 1:nspecies){
          beta[p, k] ~ dnorm(mu_beta[p], var = sigma2_beta[p])
        }
      }
      for(p in 1:p_alpha){
        mu_alpha[p] ~ dnorm(0, sd = 2)
        sigma2_alpha[p] ~ T(dnorm(0, sd = 2), 0, Inf)
      }
      for(p in 1:p_alpha){
        for(k in 1:nspecies){
          alpha[p, k] ~ dnorm(mu_alpha[p], var = sigma2_alpha[p])
        }
      }
      
      for(k in 1:nspecies){
        theta[k, 1:nspecies] ~ ddirch(alpha = alpha0[k, 1:nspecies])
      }
      
      # likelihood - site level occupancy
      for(site in 1:nsites){
        for(k in 1:nspecies){
          logit(psi[site, k]) <- (beta[1:p_beta, k] %*% X[site, 1:p_beta])[1,1]
          z[site, k] ~ dbern(psi[site, k])
        }
      }
      
      # likelihood 1
      for(row in 1:n1){
        # calculate lambda
        log(lambda[row]) <- (alpha[1:p_alpha, true_spp1[row]] %*% W1[row, 1:p_alpha])[1,1]
        
        # likelihood 
        y1[row] ~ dpois(
          z[site1[row],true_spp1[row]] *
            lambda[row] *
            theta[true_spp1[row],id_spp1[row]]
        )
      }
      
      # likelihood 2
      for(row in 1:n2){
        for(true in 1:nspecies){
          # calculate lambda
          log(lambda2[row, true]) <- (alpha[1:p_alpha, true] %*% W2[row, 1:p_alpha])[1,1]
          mean2[row, true] <- z[site2[row],true]*lambda2[row, true]*theta[true,id_spp2[row]]
        }
        y2[row] ~ dpois(sum(mean2[row, 1:nspecies]))
      }
      
    })
    
  } else{
    code <- nimbleCode({
      # priors
      for(k in 1:nspecies){
        for(p in 1:p_beta){
          beta[p, k] ~ dnorm(0, sd = 2)
        }
        for(p in 1:p_alpha){
          alpha[p, k] ~ dnorm(0, sd = 2)
        }
        theta[k, 1:nspecies] ~ ddirch(alpha = alpha0[k, 1:nspecies])
      }
      
      # likelihood - site level occupancy
      for(site in 1:nsites){
        for(k in 1:nspecies){
          logit(psi[site, k]) <- (beta[1:p_beta, k] %*% X[site, 1:p_beta])[1,1]
          z[site, k] ~ dbern(psi[site, k])
        }
      }
      
      # likelihood 1
      for(row in 1:n1){
        # calculate lambda
        log(lambda[row]) <- (alpha[1:p_alpha, true_spp1[row]] %*% W1[row, 1:p_alpha])[1,1]
        
        # likelihood 
        y1[row] ~ dpois(
          z[site1[row],true_spp1[row]] *
            lambda[row] *
            theta[true_spp1[row],id_spp1[row]]
        )
      }
      
      # likelihood 2
      for(row in 1:n2){
        for(true in 1:nspecies){
          # calculate lambda
          log(lambda2[row, true]) <- (alpha[1:p_alpha, true] %*% W2[row, 1:p_alpha])[1,1]
          mean2[row, true] <- z[site2[row],true]*lambda2[row, true]*theta[true,id_spp2[row]]
        }
        y2[row] ~ dpois(sum(mean2[row, 1:nspecies]))
      }
      
    })
    
  }
  
  # model housekeeping
  alpha0 <- matrix(.1, nspecies, nspecies)
  if(parallel) {
    library(parallel)
    this_cluster <- makeCluster(3)
    fit <- parLapply(
      cl = this_cluster,
      X = 1:nchains,
      fun = fit_model,
      code = code,
      data = list(
        # likelihood 1
        y1 = lik1_df$count,
        
        # likelihood 2
        y2 = lik2_df$count,
        
        # theta prior
        alpha0 = alpha0,
        
        # covariates
        X = X,
        W1 = W1,
        W2 = W2
      ),
      constants = list(
        p_beta = ncol(X),
        p_alpha = ncol(W1),
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
      niter = niter,
      nchains = 1,
      thin = 1,
      nburnin = nburnin,
      hierarchical = hierarchical,
      inits = list(theta = theta_inits, beta = beta_inits, alpha = alpha_inits)
    )
    stopCluster(this_cluster)
    
  } else{
    fit <- fit_model(
      seed = 1:nchains,
      code = code,
      data = list(
        # likelihood 1
        y1 = lik1_df$count,
        
        # likelihood 2
        y2 = lik2_df$count,
        
        # theta prior
        alpha0 = alpha0,
        
        # covariates
        X = X,
        W1 = W1,
        W2 = W2
      ),
      constants = list(
        p_beta = ncol(X),
        p_alpha = ncol(W1),
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
      niter = niter,
      nchains = nchains,
      thin = 1,
      nburnin = nburnin,
      hierarchical = hierarchical,
      inits = list(theta = theta_inits, beta = beta_inits, alpha = alpha_inits)
    )
  }
  
  if(save) saveRDS(fit, file = filename)
  
  out <- fit
  return(list(fit = fit, site_df = site_df, lik1_df = lik1_df, lik2_df = lik2_df))
}

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

run_sims <- function(scenarios, bc_final){
  # storage
  out <- list()
  out2 <- list()
  
  # loop through scenarios
  start <- Sys.time()
  message("Beginning simulations.")
  for(i in 1:nrow(scenarios)){
    # generate data
    tmp <- cd_format(
      data = bc_final,
      seed = scenarios[i,4] %>% as.numeric, 
      propsites_unconf = scenarios[i,2] %>% as.numeric,
      propvisits_unconf = scenarios[i,3] %>% as.numeric,
      remove_other = FALSE
    )
    
    data <- tmp %>%
      dplyr::group_by(grts, year) %>%
      mutate(
        site_num = cur_group_id(),
        id_spp_num = as.numeric(id_spp),
        true_spp_num = as.numeric(factor(true_spp))
      ) %>%
      select(grts, year, site_num, quad, visit, sample_night, id_spp, id_spp_num, true_spp, true_spp_num, everything()) %>%
      rename(
        site = site_num,
        id_spp_char = id_spp,
        id_spp = id_spp_num,
        true_spp_char = true_spp,
        true_spp = true_spp_num
      )
    
    # summarize data
    saveRDS(data, file = paste0("simulation 2a/data/data_", scenarios[i,5], ".rds"))
    tmp2 <- data %>%
      ungroup %>%
      dplyr::filter(type == "unambiguous") %>%
      group_by(id_spp_char) %>%
      summarize(total_calls = sum(count)) %>%
      mutate(
        scenario = scenarios[i,1] %>% as.numeric,
        seed = scenarios[i,4] %>% as.numeric
      ) %>%
      select(scenario, seed, everything())
    out2[[i]] <- tmp2
    
    # fit model
    fit <- cd_model_nophi_covs(
      cd_data = data,
      save = T, 
      filename = paste0("simulation 2a/models/fit_", scenarios[i,5], ".rds"), 
      parallel = T,
      niter = 5000,
      nchains = 3, 
      nburnin = 0,
      hierarchical = T,
      site_model = ~ year + mean_elevation + annual_precip + annual_mean_temp,
      detection_model = ~ min_air_temp + precip + lunar_illum 
    )
    fit <- fit[[1]]
    
    # summarize fit
    sum <- nimble_summary(fit)
    out[[i]] <- tibble(
      param = rownames(sum),
      mean = sum[,1],
      lwr = sum[,5],
      upr = sum[,9],
      rhat = sum[,10],
      ess = sum[,11],
      width = upr - lwr
    ) %>%
      mutate(
        scenario = scenarios[i,1] %>% as.numeric,
        seed = scenarios[i,4] %>% as.numeric
      ) %>%
      select(scenario, seed, everything())
    
    # output
    saveRDS(out, file = "simulation 2a/out.rds")
    saveRDS(out2, file = "simulation 2a/out2.rds")
    
    # progress
    current <- Sys.time()
    dif <- current - start
    message(paste0("Simulation ", i, " of ", nrow(scenarios), " complete. Current runtime: ", round(as.numeric(dif), 3), " ", units(dif)))
  }
  
  return(
    list(
      models = out,
      num_vetted = out2
    )
  )
}

# run sims
scenarios <- tibble(
  prop_sites_unconf = c(0, 0, 0, .25, .25, .25, .5, .5, .5),
  prop_visits_unconf = c(.25, .5, .75, .25, .5, .75, .25, .5, .75)
) %>%
  mutate(
    scenario = 1:n()
  ) %>%
  select(scenario, everything()) %>%
  mutate(freq = 3) %>%
  slice(rep(seq_len(n()), freq)) %>% 
  select(-freq) %>%
  mutate(seed = rep(c(1000, 2000, 3000), 9)) %>%
  mutate(fit = 1:n())
sim <- run_sims(scenarios = scenarios, bc_final = bc_final)
saveRDS(sim, file = "sim2a.rds")

# summarize sims
sim <- readRDS("sim2a.rds")
# sim <- readRDS("cori simulations/out.rds")
# sim2 <- readRDS("cori simulations/out2.rds")
sim_models <- do.call("rbind", sim$models)
# sim_models <- do.call("rbind", sim)
sim_vetted <- do.call("rbind", sim$num_vetted)
# sim_vetted <- do.call("rbind", sim2)
sim2a_addquant <- readRDS("sim2a_addquant.rds")
sim_models <- left_join(
  sim_models, 
  sim2a_addquant, 
  by = c("scenario", "seed", "param")
)


tmp <- sim_vetted %>%
  pivot_wider(names_from = "id_spp_char", values_from = "total_calls") %>%
  mutate(
    scenario = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  )

print(xtable::xtable(tmp, digits = 0), include.rownames = F)

rhat_tbl <- sim_models %>%
  group_by(scenario, seed) %>%
  mutate(scenario = factor(scenario), seed = factor(seed)) %>%
  summarize(max_rhat = max(rhat)) 

# add quantiles :(
# add_quant <- list()
# nimble_summary <- function(fit, warmup = nrow(fit[[1]])/2, thin = 1){
#   # convert to coda for normal summary
#   fit_warmup <- lapply(fit, function(x) x[(warmup+1):nrow(x),])
#   coda_samples <- as.mcmc.list(lapply(fit_warmup, function(x) as.mcmc(
#     x, start = warmup+1, end = nrow(fit), thin = thin
#   )))
# 
#   sum <- summary(coda_samples)
#   params <- dimnames(sum$statistics)[[1]]
#   tmp_sum <- cbind(sum$statistics, sum$quantiles)
# 
#   # get r hat / n_eff
#   mat <- matrix(NA, nrow = nrow(tmp_sum), ncol = 3)
#   colnames(mat) <- c("Rhat", "ess_bulk", "ess_tail")
#   for(i in 1:nrow(tmp_sum)){
#     tmp <- sapply(fit, function(x) x[,i])
#     mat[i,] <- c(Rhat(tmp), ess_bulk(tmp), ess_tail(tmp))
#   }
# 
#   # out
#   out <- cbind(tmp_sum, mat)
#   return(out)
# }
# for(i in 1:27){
#   fit <- readRDS(paste0("simulation 2a/models/fit_", i, ".rds"))
#   sum <- nimble_summary(fit)
#   add_quant[[i]] <- tibble(
#     param = rownames(sum),
#     q1 = sum[,6],
#     q3 = sum[,8],
#     model_num = i,
#     scenario = ceiling(i/3)
#   ) %>%
#     mutate(
#       seed = case_when(
#         i %% 3 == 1 ~ 1000,
#         i %% 3 == 2 ~ 2000,
#         i %% 3 == 0 ~ 3000
#       )
#     ) %>%
#     select(-model_num)
# }
# sim2a_addquant <- do.call("rbind", add_quant)
# saveRDS(sim2a_addquant, file = "sim2a_addquant.rds")

# occupancy prob plots
# get covariates
model.matrix(
  ~ year + mean_elevation + annual_precip + annual_mean_temp, 
  data = readRDS("simulation 2a/data/data_1.rds") %>%
    ungroup %>%
    mutate(year = factor(year)) %>%
    select(
      site, year, mean_elevation, annual_precip, annual_mean_temp
    ) %>%
    distinct()
) %>% colnames()

model.matrix( ~ min_air_temp + precip + lunar_illum , data = readRDS("simulation 2a/data/data_1.rds")) %>% colnames()

p <- sim_models %>%
  filter(grepl("beta", param)) %>%
  filter(!grepl("mu", param), !grepl("sigma", param)) %>%
  mutate(
    cov_num = substring(param, regexpr("[[]", param) + 1, regexpr(",", param) - 1),
    species = substring(param, regexpr(",", param) + 2, regexpr("[]]", param) - 1)
  ) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
  ) %>%
  mutate(
    scenario = factor(scenario),
    seed = factor(seed)
  ) %>%
  mutate(
    cov_clean = case_when(
      cov_num == 1 ~ "Intercept",
      cov_num == 2 ~ "Year2017",
      cov_num == 3 ~ "Year2018",
      cov_num == 4 ~ "Year2019",
      cov_num == 5 ~ "Year2020",
      cov_num == 6 ~ "Elevation",
      cov_num == 7 ~ "Annual precip",
      cov_num == 8 ~ "Annual temp"
    )
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
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = factor(scenario)),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  filter(cov_clean != "Intercept") %>%
  left_join(
    ., 
    rhat_tbl, 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  filter(species_char %in% c("LACI", "MYLU", "MYEV", "MYCI")) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = 1,
    position = position_dodge(.75)
  ) +
  facet_grid(cov_clean ~ species_char, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Vetting scenario",
    y = "Posterior distribution",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "black", high = "blue") +
  coord_flip()
ggsave("simulation 2a/plots/sim2a_occ.png", p, width = 2800, units = "px")

p <- sim_models %>%
  filter(grepl("beta", param)) %>%
  filter(!grepl("mu", param), !grepl("sigma", param)) %>%
  mutate(
    cov_num = substring(param, regexpr("[[]", param) + 1, regexpr(",", param) - 1),
    species = substring(param, regexpr(",", param) + 2, regexpr("[]]", param) - 1)
  ) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
  ) %>%
  mutate(
    scenario = factor(scenario),
    seed = factor(seed)
  ) %>%
  mutate(
    cov_clean = case_when(
      cov_num == 1 ~ "Intercept",
      cov_num == 2 ~ "Year2017",
      cov_num == 3 ~ "Year2018",
      cov_num == 4 ~ "Year2019",
      cov_num == 5 ~ "Year2020",
      cov_num == 6 ~ "Elevation",
      cov_num == 7 ~ "Annual precip",
      cov_num == 8 ~ "Annual temp"
    )
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
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = factor(scenario)),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  filter(cov_clean != "Intercept") %>%
  left_join(
    ., 
    rhat_tbl, 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = 1,
    position = position_dodge(.75)
  ) +
  facet_grid(cov_clean ~ species_char, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Vetting scenario",
    y = "Posterior distribution",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "red", high = "blue") +
  coord_flip()
ggsave("simulation 2a/plots/sim2a_occ_all.png", p, width = 2800, units = "px")

# relative activity plots
p <- sim_models %>%
  filter(grepl("alpha", param)) %>%
  filter(!grepl("mu", param), !grepl("sigma", param)) %>%
  mutate(
    cov_num = substring(param, regexpr("[[]", param) + 1, regexpr(",", param) - 1),
    species = substring(param, regexpr(",", param) + 2, regexpr("[]]", param) - 1)
  ) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
  ) %>%
  mutate(
    scenario = factor(scenario),
    seed = factor(seed)
  ) %>%
  mutate(
    cov_clean = case_when(
      cov_num == 1 ~ "Intercept",
      cov_num == 2 ~ "Nightly air temp",
      cov_num == 3 ~ "Nightly precip",
      cov_num == 4 ~ "Lunar illumination"
    )
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
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = factor(scenario)),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  filter(cov_clean != "Intercept") %>%
  left_join(
    ., 
    rhat_tbl, 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  filter(species_char %in% c("LACI", "MYLU", "MYEV", "MYCI")) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = 1,
    position = position_dodge(.75)
  ) +
  facet_grid(cov_clean ~ species_char, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Vetting scenario",
    y = "Posterior distribution",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "red", high = "blue") +
  coord_flip()
ggsave("simulation 2a/plots/sim2a_relact.png", p, width = 2800, units = "px")

p <- sim_models %>%
  filter(grepl("alpha", param)) %>%
  filter(!grepl("mu", param), !grepl("sigma", param)) %>%
  mutate(
    cov_num = substring(param, regexpr("[[]", param) + 1, regexpr(",", param) - 1),
    species = substring(param, regexpr(",", param) + 2, regexpr("[]]", param) - 1)
  ) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
  ) %>%
  mutate(
    scenario = factor(scenario),
    seed = factor(seed)
  ) %>%
  mutate(
    cov_clean = case_when(
      cov_num == 1 ~ "Intercept",
      cov_num == 2 ~ "Nightly air temp",
      cov_num == 3 ~ "Nightly precip",
      cov_num == 4 ~ "Lunar illumination"
    )
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
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = factor(scenario)),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  filter(cov_clean != "Intercept") %>%
  left_join(
    ., 
    rhat_tbl, 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = 1,
    position = position_dodge(.75)
  ) +
  facet_grid(cov_clean ~ species_char, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Vetting scenario",
    y = "Posterior distribution",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "red", high = "blue") +
  coord_flip()
ggsave("simulation 2a/plots/sim2a_relact_all.png", p, width = 2800, units = "px")

# theta plots
p <- sim_models %>%
  filter(grepl("theta", param)) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
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
  )  %>%
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = scenario),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  left_join(
    ., 
    rhat_tbl %>% 
      mutate(scenario = as.numeric(as.character(scenario))) %>%
      mutate(seed = as.numeric(as.character(seed))), 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  filter(
    species_true_char %in% c("LACI", "MYLU", "MYEV", "MYCI"),
    species_id_char %in% c("LACI", "MYLU", "MYEV", "MYCI")
  ) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = .75,
    position = position_dodge(.75)
  ) +
  facet_grid(species_true_char ~ species_id_char, scales = "free") +
  theme_bw() +
  labs(
    x = "Vetting scenario",
    y = "Posterior distribution",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient(low = "red", high = "blue") +
  coord_flip()
ggsave("simulation 2a/plots/sim2a_theta.png", p, width = 2800, units = "px")

p1 <- sim_models %>%
  filter(grepl("theta", param)) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
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
  )  %>%
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = scenario),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  left_join(
    ., 
    rhat_tbl %>% 
      mutate(scenario = as.numeric(as.character(scenario))) %>%
      mutate(seed = as.numeric(as.character(seed))), 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  filter(
    # species_true_char %in% c("LACI", "MYLU", "MYEV", "MYCI"),
    species_true_char %in% c("LACI"),
    species_id_char %in% c("LACI", "MYLU", "MYEV", "MYCI")
  ) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = .75,
    position = position_dodge(.75)
  ) +
  facet_grid(species_true_char ~ species_id_char, scales = "free") +
  theme_bw() +
  labs(
    x = "",
    y = "",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "red", high = "blue") +
  coord_flip() 

p2 <- sim_models %>%
  filter(grepl("theta", param)) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
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
  )  %>%
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = scenario),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  left_join(
    ., 
    rhat_tbl %>% 
      mutate(scenario = as.numeric(as.character(scenario))) %>%
      mutate(seed = as.numeric(as.character(seed))), 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  filter(
    # species_true_char %in% c("LACI", "MYLU", "MYEV", "MYCI"),
    species_true_char %in% c("MYCI"),
    species_id_char %in% c("LACI", "MYLU", "MYEV", "MYCI")
  ) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = .75,
    position = position_dodge(.75)
  ) +
  facet_grid(species_true_char ~ species_id_char, scales = "free") +
  theme_bw() +
  labs(
    x = "",
    y = "",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "red", high = "blue") +
  coord_flip()


p3 <- sim_models %>%
  filter(grepl("theta", param)) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
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
  )  %>%
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = scenario),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  left_join(
    ., 
    rhat_tbl %>% 
      mutate(scenario = as.numeric(as.character(scenario))) %>%
      mutate(seed = as.numeric(as.character(seed))), 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  filter(
    # species_true_char %in% c("LACI", "MYLU", "MYEV", "MYCI"),
    species_true_char %in% c("MYEV"),
    species_id_char %in% c("LACI", "MYLU", "MYEV", "MYCI")
  ) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = .75,
    position = position_dodge(.75)
  ) +
  facet_grid(species_true_char ~ species_id_char, scales = "free") +
  theme_bw() +
  labs(
    x = "",
    y = "",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "red", high = "blue") +
  coord_flip() 

p4 <- sim_models %>%
  filter(grepl("theta", param)) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
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
  )  %>%
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = scenario),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  left_join(
    ., 
    rhat_tbl %>% 
      mutate(scenario = as.numeric(as.character(scenario))) %>%
      mutate(seed = as.numeric(as.character(seed))), 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  filter(
    # species_true_char %in% c("LACI", "MYLU", "MYEV", "MYCI"),
    species_true_char %in% c("MYLU"),
    species_id_char %in% c("LACI", "MYLU", "MYEV", "MYCI")
  ) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = .75,
    position = position_dodge(.75)
  ) +
  facet_grid(species_true_char ~ species_id_char, scales = "free") +
  theme_bw() +
  labs(
    x = "",
    y = "",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "red", high = "blue") +
  coord_flip()

p5 <- gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 4)
ggsave("simulation 2a/plots/sim2a_theta.png", p5, width = 2800, units = "px")

p6 <- gridExtra::grid.arrange(
  p1 +
    geom_point(
      data = tibble(
        species_true_char = "LACI", 
        species_id_char = rep(c("LACI", "MYCI", "MYEV", "MYLU"), each = 2),
        mean = c(.85, 1, 0, .05, 0, .05, 0, .05),
        scenario_clean = "HH"
      ), aes(x = scenario_clean, y = mean), alpha = 0
    ),
  p3 +
    geom_point(
      data = tibble(
        species_true_char = "MYEV", 
        species_id_char = rep(c("LACI", "MYCI", "MYEV", "MYLU"), each = 2),
        mean = c(0, .05, 0, .05, .85, 1, 0, .05),
        scenario_clean = "HH"
      ), aes(x = scenario_clean, y = mean), alpha = 0
    ),
  nrow = 2
)
ggsave("simulation 2a/plots/sim2a_theta_fixed.png", p6, width = 2800, units = "px")

p <- sim_models %>%
  filter(grepl("theta", param)) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
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
  )  %>%
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = scenario),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  left_join(
    ., 
    rhat_tbl %>% 
      mutate(scenario = as.numeric(as.character(scenario))) %>%
      mutate(seed = as.numeric(as.character(seed))), 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = .75,
    position = position_dodge(.75)
  ) +
  facet_grid(species_true_char ~ species_id_char, scales = "free") +
  theme_bw() +
  labs(
    x = "Vetting scenario",
    y = "Posterior distribution",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_color_gradient(low = "red", high = "blue") +
  coord_flip()
ggsave("simulation 2a/plots/sim2a_theta_all.png", p, width = 2800, units = "px")

p <- sim_models %>%
  filter(grepl("theta", param)) %>%
  left_join(
    ., 
    sim_vetted %>% group_by(scenario, seed) %>% summarize(total_vetted = sum(total_calls)),
    by = c("scenario", "seed")
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
  )  %>%
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) %>%
  left_join(
    ., 
    sim_vetted %>% 
      group_by(scenario, seed) %>% 
      summarize(total_vetted = sum(total_calls)) %>% 
      group_by(scenario) %>% 
      summarize(mean_effort = mean(total_vetted)) %>%
      mutate(scenario = scenario),
    by = "scenario"
  ) %>%
  mutate(
    scenario_clean_ = forcats::fct_reorder(scenario_clean, mean_effort),
    scenario_clean = factor(scenario_clean, levels = levels(scenario_clean_))
  ) %>%
  left_join(
    ., 
    rhat_tbl %>% 
      mutate(scenario = as.numeric(as.character(scenario))) %>%
      mutate(seed = as.numeric(as.character(seed))), 
    by = c("scenario", "seed")
  ) %>%
  filter(max_rhat <= 1.5) %>%
  filter(
    species_true_char %in% c("LACI", "MYLU", "MYEV", "MYCI"),
    species_id_char %in% c("LACI", "MYLU", "MYEV", "MYCI")
  ) %>%
  ggplot() +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = lwr, ymax = upr, group = seed),
    position = position_dodge(.75)
  ) +
  geom_linerange(
    aes(x = scenario_clean, y = mean, ymin = q1, ymax = q3, color = total_vetted, group = seed),
    lwd = 1.5,
    position = position_dodge(.75)
  ) +
  geom_point(
    aes(x = scenario_clean, y = mean, group = seed), 
    size = .75,
    position = position_dodge(.75)
  ) +
  facet_grid(species_true_char ~ species_id_char, scales = "free") +
  theme_bw() +
  labs(
    x = "Vetting scenario",
    y = "Posterior distribution",
    col = "Total vetted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient(low = "red", high = "blue") +
  coord_flip()
ggsave("simulation 2a/plots/sim2a_theta.png", p, width = 2800, units = "px")

# look at empirical classification probs for each data set
# out <- list()
# for(i in 1:27){
#   data <- readRDS(paste0("simulation 2a/data/data_", i, ".rds"))
#   out[[i]] <- data %>% filter(type == "unambiguous") %>%
#     group_by(true_spp_char, id_spp_char) %>%
#     summarize(total = sum(count)) %>%
#     group_by(true_spp_char) %>%
#     mutate(total_calls = sum(total)) %>%
#     mutate(theta = total / total_calls) %>%
#     mutate(
#       scenario = ceiling(i/3),
#       seed = case_when(
#         i %% 3 == 1 ~ 1000,
#         i %% 3 == 2 ~ 2000,
#         i %% 3 == 0 ~ 3000
#       )
#     )
# }
# emp_thetas <- do.call("rbind", out)
# saveRDS(emp_thetas, file = "simulation 2a/emp_thetas.rds")
emp_thetas <- readRDS("simulation 2a/emp_thetas.rds") %>%
  mutate(
    scenario_clean = case_when(
      scenario == 1 ~ "HH",
      scenario == 2 ~ "HM",
      scenario == 3 ~ "HL",
      scenario == 4 ~ "MH",
      scenario == 5 ~ "MM",
      scenario == 6 ~ "ML",
      scenario == 7 ~ "LH",
      scenario == 8 ~ "LM",
      scenario == 9 ~ "LL"
    )
  ) 

emp_thetas %>%
  filter(
    true_spp_char %in% c("LACI", "MYLU", "MYEV", "MYCI"),
    id_spp_char %in% c("LACI", "MYLU", "MYEV", "MYCI")
  ) %>%
  ggplot() +
  geom_point(aes(x = scenario_clean, y = theta, group = seed), position = position_dodge(.75)) +
  facet_grid(true_spp_char ~ id_spp_char) +
  theme_bw() +
  coord_flip()


emp_thetas %>%
  filter(
    true_spp_char %in% c("LACI", "MYLU", "MYEV", "MYCI"),
    id_spp_char %in% c("LACI", "MYLU", "MYEV", "MYCI")
  ) %>% 
  filter(true_spp_char == id_spp_char) %>%
  summarize(min_diag_theta = min(theta))

left_join(
  sim_vetted %>% mutate(scenario = factor(scenario), seed = factor(seed)),
  rhat_tbl, 
  by = c("scenario", "seed")
) %>%
  filter(max_rhat >= 1.5) %>%
  group_by(scenario, seed) %>%
  summarize(min_calls = min(total_calls)) %>%
  left_join(
    ., 
    sim_vetted %>% mutate(scenario = factor(scenario), seed = factor(seed), min_calls = total_calls)
  )

# model summaries
data <- readRDS("simulation 2a/data/data_1.rds")
model1 <- readRDS("simulation 2a/models/fit_4.rds")
# nimble_summary(model0)

# traceplots 
# plot
nchains <- length(model1)
niter <- nrow(model1[[1]])
tmp <- do.call("rbind", model1)
plot_tbl <- tibble(
  trace = c(tmp),
  param = rep(colnames(tmp), each = nchains*niter),
  chain = factor(rep(rep(1:nchains, each = niter), ncol(tmp))),
  iteration = rep(rep(1:niter, nchains), ncol(tmp))
)

sum_tbl <- plot_tbl %>%
  filter(!grepl("mu", param), !grepl("sigma", param)) %>%
  filter(iteration >= 2500) %>%
  group_by(param) %>%
  summarize(
    mean = mean(trace),
    lwr = quantile(trace, .025),
    upr = quantile(trace, .975),
    q1 = quantile(trace, .25),
    q3 = quantile(trace, .75)
  )

p <- sum_tbl %>%
  filter(grepl("beta", param)) %>%
  mutate(
    cov_num = substring(param, regexpr("[[]", param) + 1, regexpr(",", param) - 1),
    species = substring(param, regexpr(",", param) + 2, regexpr("[]]", param) - 1)
  ) %>%
  left_join(
    .,
    tibble(
      species_char = as.character(sort(unique(data$id_spp_char)))
    ) %>% mutate(species = 1:n() %>% as.character)
  ) %>%
  left_join(
    .,
    tibble(
      cov = model.matrix(
        ~ year + mean_elevation + annual_precip + annual_mean_temp,
        data = readRDS("simulation 2a/data/data_1.rds") %>%
          ungroup %>%
          mutate(year = factor(year)) %>%
          select(
            site, year, mean_elevation, annual_precip, annual_mean_temp
          ) %>%
          distinct()
      ) %>% colnames()
    ) %>% mutate(cov_num = 1:n() %>% as.character)
  ) %>%
  mutate(
    cov_clean = case_when(
      cov == "(Intercept)" ~ "Intercept",
      cov == "year2017" ~ "Year2017",
      cov == "year2018" ~ "Year2018",
      cov == "year2019" ~ "Year2019",
      cov == "year2020" ~ "Year2020",
      cov == "mean_elevation" ~ "Elevation",
      cov == "annual_precip" ~ "Precipitation",
      cov == "annual_mean_temp" ~ "Temperature",
      TRUE ~ cov
    )
  ) %>%
  mutate(cov_clean = factor(cov_clean, levels = c("Intercept", "Year2017", "Year2018", "Year2019", "Year2020", "Elevation", "Precipitation", "Temperature"))) %>%
  ggplot() +
  # geom_boxplot(
  #   aes(x = cov_clean, lower = q1, upper = q3, middle = mean, ymin = lwr, ymax = upr, fill = species_char, width = .5),
  #   stat = "identity"
  # ) +
  geom_linerange(
    aes(x = species_char, y = mean, ymin = lwr, ymax = upr)
  ) +
  geom_linerange(
    aes(x = species_char, y = mean, ymin = q1, ymax = q3),
    lwd = 1.5
  ) +
  geom_point(
    aes(x = species_char, y = mean), size = 2
  ) +
  facet_grid(rows = vars(cov_clean)) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Species",
    y = "Posterior distribution"
  ) +
  theme(legend.position = "none") +
  geom_hline(aes(yintercept = 0), linetype = "dotdash")
ggsave("simulation 2a/plots/model_occ_covpanel.png", p, width = 2800, units = "px")

p <- sum_tbl %>%
  filter(grepl("alpha", param)) %>%
  mutate(
    cov_num = substring(param, regexpr("[[]", param) + 1, regexpr(",", param) - 1),
    species = substring(param, regexpr(",", param) + 2, regexpr("[]]", param) - 1)
  ) %>%
  left_join(
    ., 
    tibble(
      species_char = as.character(sort(unique(data$id_spp_char)))
    ) %>% mutate(species = 1:n() %>% as.character)
  ) %>%
  left_join(
    .,
    tibble(
      cov = model.matrix(
        ~ min_air_temp + precip + lunar_illum,
        data = readRDS("simulation 2a/data/data_1.rds") %>%
          ungroup %>%
          mutate(year = factor(year)) %>%
          select(
            min_air_temp, precip, lunar_illum
          ) %>%
          distinct()
      ) %>% colnames()
    ) %>% mutate(cov_num = 1:n() %>% as.character)
  ) %>%
  mutate(
    cov_clean = case_when(
      cov == "(Intercept)" ~ "Intercept",
      cov == "min_air_temp" ~ "Temperature",
      cov == "precip" ~ "Precipitation",
      cov == "lunar_illum" ~ "Illumination",
      TRUE ~ cov
    )
  ) %>%
  mutate(cov_clean = factor(cov_clean, levels = rev(c("Intercept", "Temperature", "Precipitation", "Illumination")))) %>%
  ggplot() +
  # geom_boxplot(
  #   aes(x = cov_clean, lower = q1, upper = q3, middle = mean, ymin = lwr, ymax = upr, fill = species_char, width = .5),
  #   stat = "identity"
  # ) +
  geom_linerange(
    aes(x = species_char, y = mean, ymin = lwr, ymax = upr)
  ) +
  geom_linerange(
    aes(x = species_char, y = mean, ymin = q1, ymax = q3), #, color = cov_clean
    lwd = 1.5
  ) +
  geom_point(
    aes(x = species_char, y = mean), size = 1
  ) +
  # facet_wrap(~ species_char, scales = "free_x", ncol = 1, strip.position = "right") +
  facet_grid(rows = vars(cov_clean)) + 
  coord_flip() +
  theme_bw() +
  labs(
    x = "Species",
    y = "Posterior distribution"
  ) +
  theme(legend.position = "none") +
  geom_hline(aes(yintercept = 0), linetype = "dotdash")
ggsave("simulation 2a/plots/model_relact_covpanel.png", p, width = 2800, units = "px")

p <- sum_tbl %>%
  filter(grepl("theta", param)) %>%
  mutate(
    true_spp = substring(param, regexpr("[[]", param) + 1, regexpr(",", param) - 1),
    id_spp = substring(param, regexpr(",", param) + 2, regexpr("[]]", param) - 1)
  ) %>% 
  left_join(
    .,
    tibble(
      true_spp_chr = sort(unique(data$true_spp_char))
    ) %>% mutate(true_spp = 1:n() %>% as.character)
  ) %>%
  left_join(
    .,
    tibble(
      id_spp_chr = sort(unique(data$id_spp_char))
    ) %>% mutate(id_spp = 1:n() %>% as.character)
  ) %>%
  mutate(
    correct = case_when(
      true_spp_chr == id_spp_chr ~ "yes",
      TRUE ~ "no"
    )
  ) %>%
  ggplot() + 
  geom_bar(aes(x = true_spp_chr, y = mean, fill = id_spp_chr, col = correct), size = 1.1, position = "fill", stat = "identity") +
  theme_bw() +
  labs(
    x = "manual",
    fill = "auto",
    y = "proportion"
  ) +
  scale_color_manual(values = c(NA, "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE, size = FALSE) +
  labs(
    y = "Classification probability",
    x = "Manually identified spp",
    fill = "Identified spp"
  )
ggsave("simulation 2a/plots/model_theta.png", p, width = 2800, units = "px")

sum_tbl %>%
  filter(grepl("theta", param)) %>%
  mutate(
    true_spp = substring(param, regexpr("[[]", param) + 1, regexpr(",", param) - 1),
    id_spp = substring(param, regexpr(",", param) + 2, regexpr("[]]", param) - 1)
  ) %>% 
  left_join(
    .,
    tibble(
      true_spp_chr = sort(unique(data$true_spp_char))
    ) %>% mutate(true_spp = 1:n() %>% as.character)
  ) %>%
  left_join(
    .,
    tibble(
      id_spp_chr = sort(unique(data$id_spp_char))
    ) %>% mutate(id_spp = 1:n() %>% as.character)
  ) %>%
  mutate(
    correct = case_when(
      true_spp_chr == id_spp_chr ~ "yes",
      TRUE ~ "no"
    )
  ) %>%
  filter(true_spp_chr == id_spp_chr)

# naive theta
data %>%
  filter(type == "unambiguous") %>%
  select(true_spp_char, id_spp_char, count) %>%
  group_by(true_spp_char, id_spp_char) %>%
  summarize(total = sum(count)) %>%
  mutate(correct = ifelse(id_spp_char == true_spp_char, T, F)) %>%
  ggplot() + 
  geom_bar(aes(x = true_spp_char, y = total, fill = id_spp_char, col = correct), size = 1.1, position = "fill", stat = "identity") +
  theme_bw() +
  labs(
    x = "manual",
    fill = "auto",
    y = "proportion"
  ) +
  scale_color_manual(values = c(NA, "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE, size = FALSE) +
  labs(
    y = "Classification probability",
    x = "Manually identified spp",
    fill = "Identified spp"
  )
