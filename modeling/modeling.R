# clear workspace
rm(list = ls())

# packages
packs <- c("dplyr", "nimble", "htmltools", "ggplot2", "sf", "Rcpp", "RcppArmadillo", "inline", "mvtnorm", "readr", "parallel", "xtable", "rstan", "coda", "vegan", "tidyverse", "lubridate", "bcmaps", "knitr", "kableExtra", "ggalluvial")
sapply(packs, require, character.only = T)
rm(packs)

# convenience
`%notin%` <- Negate("%in%")
options(mc.cores = parallel::detectCores())

###################################
### Naive model - no covariates ###
###################################
# from data prep
bc_final <- readRDS("bc_final.rds")

# format data to pass to model
cd_format <- function(data, seed = NULL){
  # function to format data for CD model
  # each row represents a single call
  # not super generalized right now - can work on 
  
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
  
  # insert NAs for manual call for some visits - unique to BC data
  tmp4 <- tmp3 %>% 
    ungroup() %>%
    select(grts, year, quad, sample_night) %>%
    distinct %>%
    group_by(grts, year) %>%
    sample_frac(.5) %>%
    mutate(man = NA) %>%
    ungroup()
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
test <- cd_format(bc_final, seed = 08052021)

data <- test %>%
  group_by(grts, year) %>%
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

cd_model_nophi <- function(cd_data, save = T, filename = "test.rds", parallel = F, niter = 5000, nchains = 3, nburnin){
  
  fit_model <- function(seed = 1, code, data, constants, inits, niter, nchains, thin = 1, nburnin = 0){
    library(nimble)
    
    # R model
    model <- nimbleModel(code, constants, data)
    
    # C model
    model_c <- compileNimble(model)
    
    # R mcmc
    model_conf <- configureMCMC(model, enableWAIC = TRUE)
    model_conf$addMonitors(c("z", "occ_resid"))
    
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
      nburnin = nburnin, 
      WAIC = TRUE
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
        
        occ_resid[site, k] <- z[site, k] - psi[k]
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
  alpha0 <- matrix(1, nspecies, nspecies)
  # diag(alpha0) <- 1
  init_func <- function(seed = 1, nspecies = 8, nrows = nrow(site_df), ncols = nspecies){
    set.seed(seed)
    
    out <- list(
      psi = runif(nspecies),
      lambda = abs(rnorm(nspecies, 0, 100)),
      theta = t(apply(alpha0, 1, function(x) rdirch(1, x))),
      z = matrix(1, nrow = nrows, ncol = ncols)
    )
    
    return(out)
  }
  
  
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
      niter = niter,
      nchains = 1,
      thin = 1,
      nburnin = nburnin
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
      inits = lapply(1:nchains, init_func, nspecies = length(unique(cd_data$id_spp))),
      niter = niter,
      nchains = nchains,
      thin = 1,
      nburnin = nburnin
    )
  }
  
  if(save) saveRDS(fit, file = filename)
  
  out <- fit
  return(out)
}
naive_fit <- cd_model_nophi(
  cd_data = data, save = T, filename = "naive_fit_alpha01.rds", parallel = F,
  nchains = 3, niter = 5000, nburnin = 2500
)

# summarize
nimble_summary <- function(fit, warmup = 0, thin = 1){
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
model0 <- readRDS("naive_fit_alpha01.rds")$samples
# nimble_summary(model0)

# traceplots 
# plot
nchains <- length(model0)
niter <- nrow(model0[[1]])
tmp <- do.call("rbind", model0)
plot_tbl <- tibble(
  trace = c(tmp),
  param = rep(colnames(tmp), each = nchains*niter),
  chain = factor(rep(rep(1:nchains, each = niter), ncol(tmp))),
  iteration = rep(rep(1:niter, nchains), ncol(tmp))
)

plot_tbl %>%
  filter(grepl("psi", param)) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = trace, col = chain)) +
  # geom_hline(data = truth_tbl %>% filter(grepl("psi", param)), aes(yintercept = truth), col = "red") +
  facet_wrap(~ param) +
  theme_bw()

plot_tbl %>%
  filter(grepl("lambda", param)) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = trace, col = chain)) +
  # geom_hline(data = truth_tbl %>% filter(grepl("psi", param)), aes(yintercept = truth), col = "red") +
  facet_wrap(~ param) +
  # ylim(0, 30) +
  theme_bw()

plot_tbl %>%
  filter(grepl("theta", param)) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = trace, col = chain)) +
  # geom_hline(data = truth_tbl %>% filter(grepl("psi", param)), aes(yintercept = truth), col = "red") +
  facet_wrap(~ param) +
  ylim(0, 1) +
  theme_bw() 

sum_tbl <- plot_tbl %>%
  filter(grepl("lambda", param) | grepl("psi", param) | grepl("theta", param)) %>%
  group_by(param) %>%
  summarize(
    mean = mean(trace),
    lwr = quantile(trace, .025),
    upr = quantile(trace, .975)
  )

sum_tbl %>%
  filter(grepl("psi", param)) %>%
  mutate(
    param = factor(param, levels = paste0("psi[", 1:10, "]")),
    species = substring(param, regexpr("[[]", param) + 1, regexpr("[]]", param) - 1)
  ) %>%
  left_join(
    ., 
    tibble(
      species_chr = sort(unique(data$id_spp_char))
    ) %>% mutate(species = 1:n() %>% as.character)
  ) %>%
  ggplot() +
  geom_pointrange(aes(x = param, y = mean, ymin = `lwr`, ymax = `upr`, col = species_chr)) +
  coord_flip() +
  theme_bw() +
  ylim(0, 1)

sum_tbl %>%
  filter(grepl("lambda", param)) %>%
  mutate(
    param = factor(param, levels = paste0("lambda[", 1:10, "]")),
    species = substring(param, regexpr("[[]", param) + 1, regexpr("[]]", param) - 1)
  ) %>%
  left_join(
    ., 
    tibble(
      species_chr = sort(unique(data$id_spp_char))
    ) %>% mutate(species = 1:n() %>% as.character)
  ) %>%
  ggplot() +
  geom_pointrange(aes(x = param, y = mean, ymin = `lwr`, ymax = `upr`, col = species_chr)) +
  coord_flip() +
  theme_bw() 

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
  guides(color = FALSE, size = FALSE)

sum <- nimble_summary(model0)

# parameter estimates for simulation
lambda_true <- sum[which(grepl("lambda", rownames(sum))),1]
names(lambda_true) <- levels(data$id_spp_char)

psi_true <- sum[which(grepl("psi", rownames(sum))),1]
names(psi_true) <- levels(data$id_spp_char)

theta_true <- sum[which(grepl("theta", rownames(sum))),1] %>% matrix(., nrow = 10, ncol = 10)
colnames(theta_true) <- rownames(theta) <- levels(data$id_spp_char)

save(lambda_true, psi_true, theta_true, file = "parameter_estimates.rdata")

############################
### Model 1 - covariates ###
############################
test <- cd_format(bc_final, seed = 08052021)
data <- test %>%
  group_by(grts, year) %>%
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

# detection_model = ~ min_air_temp + precip + lunar_illum

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
  initial_site_df <- data %>% 
    group_by(site, id_spp) %>%
    summarize(
      total_count = sum(count)
    ) %>%
    mutate(naive_occ = ifelse(total_count == 0, 0, 1)) %>%
    ungroup %>%
    left_join(
      ., 
      data %>% ungroup %>% select(site, site_covs),
      by = "site"
    ) %>%
    distinct() %>%
    mutate(year = factor(year)) %>%
    mutate_at(.vars = site_covs_, .funs = ~(scale(.) %>% as.vector))
  
  beta_inits <- matrix(NA, nrow = ncol(X), ncol = nspecies)
  for(i in 1:nspecies){
    df <- initial_site_df %>% filter(id_spp == i)
    beta_inits[,i] <- coef(glm(
      as.formula(paste0("naive_occ ~ ", paste0(site_covs, collapse = " + "))), 
      family = "binomial", 
      data = df
    ))
  }
  
  # initialize alpha
  initial_act_df <- data %>% 
    group_by(site, visit, id_spp) %>%
    summarize(
      total_count = sum(count)
    ) %>%
    mutate(naive_occ = ifelse(total_count == 0, 0, 1)) %>%
    ungroup %>%
    left_join(
      ., 
      data %>% ungroup %>% select(site, visit, act_covs),
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

Sys.time()
model1 <- cd_model_nophi_covs(
  cd_data = data,
  save = T, 
  filename = "model1.rds", 
  parallel = T,
  niter = 5000,
  nchains = 3, 
  nburnin = 0,
  hierarchical = T,
  site_model = ~ year + mean_elevation + annual_precip + annual_mean_temp,
  detection_model = ~ min_air_temp + precip + lunar_illum 
)
Sys.time()

model1 <- readRDS("model1.rds")
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

plot_tbl %>%
  filter(iteration >= 2500) %>%
  filter(grepl("beta", param)) %>%
  filter(!grepl("mu", param), !grepl("sigma", param)) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = trace, col = chain)) +
  # geom_hline(data = truth_tbl %>% filter(grepl("psi", param)), aes(yintercept = truth), col = "red") +
  facet_wrap(~ param) +
  theme_bw()

plot_tbl %>%
  filter(iteration >= 2500) %>%
  filter(grepl("alpha", param)) %>%
  filter(!grepl("mu", param), !grepl("sigma", param)) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = trace, col = chain)) +
  # geom_hline(data = truth_tbl %>% filter(grepl("psi", param)), aes(yintercept = truth), col = "red") +
  facet_wrap(~ param) +
  # ylim(0, 30) +
  theme_bw()

plot_tbl %>%
  filter(iteration >= 2500) %>%
  filter(grepl("theta", param)) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = trace, col = chain)) +
  # geom_hline(data = truth_tbl %>% filter(grepl("psi", param)), aes(yintercept = truth), col = "red") +
  facet_wrap(~ param) +
  ylim(0, 1) +
  theme_bw() 

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

sum_tbl %>%
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
      cov = model.matrix( ~ year + mean_elevation + annual_precip + annual_mean_temp, data %>% mutate(year = factor(year))) %>% colnames(.)
    ) %>% mutate(cov_num = 1:n() %>% as.character)
  ) %>%
  mutate(
    cov_clean = case_when(
      cov == "(Intercept)" ~ "Intercept",
      cov == "year2016" ~ "Year2016",
      cov == "year2017" ~ "Year2017",
      cov == "year2019" ~ "Year2019",
      cov == "year2020" ~ "Year2020",
      cov == "mean_elevation" ~ "Elevation",
      cov == "annual_precip" ~ "Precipitation",
      cov == "annual_mean_temp" ~ "Temperature",
      TRUE ~ cov
    )
  ) %>%
  mutate(cov_clean = factor(cov_clean, levels = rev(c("Year2016", "Year2017", "Intercept", "Year2019", "Year2020", "Elevation", "Precipitation", "Temperature")))) %>%
  ggplot() +
  # geom_boxplot(
  #   aes(x = cov_clean, lower = q1, upper = q3, middle = mean, ymin = lwr, ymax = upr, fill = species_char, width = .5),
  #   stat = "identity"
  # ) +
  geom_linerange(
    aes(x = cov_clean, y = mean, ymin = lwr, ymax = upr)
  ) +
  geom_linerange(
    aes(x = cov_clean, y = mean, ymin = q1, ymax = q3, color = species_char),
    lwd = 1.5
  ) +
  geom_point(
    aes(x = cov_clean, y = mean), size = 2
  ) +
  facet_grid(rows = vars(species_char)) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Occupancy probability coefficients",
    y = "Posterior distribution"
  ) +
  theme(legend.position = "none") +
  geom_hline(aes(yintercept = 0), linetype = "dotdash")

sum_tbl %>%
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
      cov = model.matrix( ~ year + mean_elevation + annual_precip + annual_mean_temp, data %>% mutate(year = factor(year))) %>% colnames(.)
    ) %>% mutate(cov_num = 1:n() %>% as.character)
  ) %>%
  mutate(
    cov_clean = case_when(
      cov == "(Intercept)" ~ "Intercept",
      cov == "year2016" ~ "Year2016",
      cov == "year2017" ~ "Year2017",
      cov == "year2019" ~ "Year2019",
      cov == "year2020" ~ "Year2020",
      cov == "mean_elevation" ~ "Elevation",
      cov == "annual_precip" ~ "Precipitation",
      cov == "annual_mean_temp" ~ "Temperature",
      TRUE ~ cov
    )
  ) %>%
  mutate(cov_clean = factor(cov_clean, levels = rev(c("Year2016", "Year2017", "Intercept", "Year2019", "Year2020", "Elevation", "Precipitation", "Temperature")))) %>%
  ggplot() +
  # geom_boxplot(
  #   aes(x = cov_clean, lower = q1, upper = q3, middle = mean, ymin = lwr, ymax = upr, fill = species_char, width = .5),
  #   stat = "identity"
  # ) +
  geom_linerange(
    aes(x = species_char, y = mean, ymin = lwr, ymax = upr)
  ) +
  geom_linerange(
    aes(x = species_char, y = mean, ymin = q1, ymax = q3, color = cov_clean),
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

sum_tbl %>%
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
      cov = model.matrix( ~ min_air_temp + precip + lunar_illum , data)  %>% colnames(.)
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
    aes(x = cov_clean, y = mean, ymin = lwr, ymax = upr)
  ) +
  geom_linerange(
    aes(x = cov_clean, y = mean, ymin = q1, ymax = q3, color = species_char),
    lwd = 1.5
  ) +
  geom_point(
    aes(x = cov_clean, y = mean), size = 1
  ) +
  # facet_wrap(~ species_char, scales = "free_x", ncol = 1, strip.position = "right") +
  facet_grid(rows = vars(species_char)) + 
  coord_flip() +
  theme_bw() +
  labs(
    x = "Relative activity coefficients",
    y = "Posterior distribution"
  ) +
  theme(legend.position = "none") +
  geom_hline(aes(yintercept = 0), linetype = "dotdash")

sum_tbl %>%
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
      cov = model.matrix( ~ min_air_temp + precip + lunar_illum , data)  %>% colnames(.)
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
    aes(x = species_char, y = mean, ymin = q1, ymax = q3, color = cov_clean),
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
    x = "Relative activity coefficients",
    y = "Posterior distribution"
  ) +
  theme(legend.position = "none") +
  geom_hline(aes(yintercept = 0), linetype = "dotdash")


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



