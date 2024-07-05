

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(runjags)
library(RWiener)
library(here)
library(multitool)
library(glue)

nobs <- 720
ncores <- 4 # Number of cores used for parallel processing

set.seed(2865905)

source('scripts/functions_DDM.R')
source('scripts/DBDA2E-utilities.R')

# I take the following six steps for the simulation analyses:

# 1. Based on observed mean RTs, accuracy rates, and correlations between task conditions, simulate 10,000 trials per participant
# 2. Estimate DDM parameters based on simulated trials (which will be highly reliable because of the large number of trials).
# 3. For each DDM parameter, calculate the mean and SD, and use these to simulate *ground-truth* DDM parameters for individual participants
# 4. Simulate RTs/accuracies based on *ground-truth* DDM parameters, with the same number of trials as the original study
# 5. Apply Hierarchical Bayesian DDM to the simulated RTs/accuracies
# 6. Compare estimated DDM parameters with simulated *ground-truth* DDM parameters, under different conditions (e.g., different model specifications)

# Flanker settings

ntrials_fm_con <- 29 # Nr of trials in congruent condition
ntrials_fm_inc <- 16 # Nr of trials in incongruent condition

ffm_tcon_m       <- 0.822
ffm_tcon_sd      <- 0.160
ffm_tinc_m       <- 0.521
ffm_tinc_sd      <- 0.155

ffm_tcon_rta_m   <- 0.820
ffm_tcon_rta_sd  <- 0.135
ffm_tinc_rta_m   <- 0.904
ffm_tinc_rta_sd  <- 0.172

ffm_rta_cor      <- 0.78
ffm_acc_cor      <- 0.28

sfm_tcon_m       <- 0.853
sfm_tcon_sd      <- 0.142
sfm_tinc_m       <- 0.581
sfm_tinc_sd      <- 0.173

sfm_tcon_rta_m   <- 0.815
sfm_tcon_rta_sd  <- 0.127
sfm_tinc_rta_m   <- 0.906
sfm_tinc_rta_sd  <- 0.163

# 1. Simulate 1,000 trials per participant based on observed RTs --------

sim_fm_raw <-
  tibble(
    id   = 1:nobs
  ) %>%
  # Simulate mean/SD RTs per subject, mean accuracy, as well as rate and shape settings for the Gamma distribution (to simulate right-skewed RTs)
  mutate(
    faux::rnorm_multi(
      n = nobs,
      varnames = c("rt_con_m", "rt_con_sd", "rt_inc_m", "rt_inc_sd"),
      mu = c(ffm_tcon_rta_m,  ffm_tcon_rta_sd, ffm_tinc_rta_m, ffm_tinc_rta_sd),
      sd = c(ffm_tcon_rta_sd, ffm_tcon_rta_sd/2, ffm_tinc_rta_sd, ffm_tinc_rta_sd/2),
      r = c(   1,         .80,    ffm_rta_cor,   .80,
               .80,           1,    .80,           .80,
               ffm_rta_cor, .80,      1,           .80,
               .80,         .80,    .80,             1)
    ),
    across(c(rt_inc_sd, rt_con_sd), ~ifelse(. < 0, . * -1, .)),

    acc_con_m = rnorm(nobs, mean = ffm_tcon_m, sd = ffm_tcon_sd),
    acc_inc_m = rnorm(nobs, mean = ffm_tinc_m, sd = ffm_tinc_sd),

    rate_con  = rt_con_m / rt_con_sd,
    shape_con = rt_con_m * rate_con,

    rate_inc  = rt_inc_m / rt_inc_sd,
    shape_inc = rt_inc_m * rate_inc,

    across(c(acc_con_m, acc_inc_m), ~ifelse(. > 1, 1, .)),
    across(c(acc_con_m, acc_inc_m), ~ifelse(. < 0, 0, .))
  ) |>
  mutate(
    trial_rts = pmap(list(id, rt_con_m, rt_con_sd, rt_inc_m, rt_inc_sd, acc_con_m, acc_inc_m, rate_con, shape_con, rate_inc, shape_inc), function(id, rt_con_m, rt_con_sd, rt_inc_m, rt_inc_sd, acc_con_m, acc_inc_m, rate_con, shape_con, rate_inc, shape_inc) {

      bind_rows(
        con_trials <- rgamma(n = 1000, shape = shape_con, rate = rate_con) |>
          as_tibble() |>
          rename(rt = value) |>
          filter(rt > 0.2) |>
          mutate(
            id = id,
            condition = "congruent",
            acc = sample(x = c(0,1), size = n(), prob = c(1-acc_con_m, acc_con_m), replace = TRUE),
            rt_con_m = rt_con_m,
            rt_con_sd = rt_con_sd
          ),
        inc_trials <- rgamma(n = 1000, shape = shape_inc, rate = rate_inc) |>
          as_tibble() |>
          rename(rt = value) |>
          filter(rt > 0.2) |>
          mutate(
            id = id,
            condition = "incongruent",
            acc = sample(x = c(0,1), size = n(), prob = c(1-acc_inc_m, acc_inc_m), replace = TRUE),
            rt_inc_m = rt_inc_m,
            rt_inc_sd = rt_inc_sd
          )
      )
    })
  ) |>
  select(id, trial_rts)



# 2. Estimate DDM parameters based on simulated trials --------------------


sim_fm_ddm <-
  sim_fm_raw |>
  select(trial_rts) |>
  unnest(trial_rts) |>
  select(id, rt, correct = acc, condition)

write_DDM_files(data = sim_fm_ddm, path = "data/flanker", vars = c("rt", "correct", "condition"), task = "flanker")


fast_dm_settings(task = "flanker",
                 path = "data/flanker",
                 model_version = "_mod2",
                 method = "ml",
                 depend = c("depends v condition", "depends a condition", "depends t0 condition"),
                 format = "TIME RESPONSE condition")

# Compute DDM parameters
execute_fast_dm(task = "flanker", path = "data/flanker", model_version = "_mod2")



# Read DDM results
sim_fm_ddm <- read_DDM(task = "flanker", path = "data/flanker", model_version = "_mod2") |>
  summarise(
    v_con_m   = mean(v_congruent),
    v_con_sd  = sd(v_congruent),
    v_inc_m   = mean(v_incongruent),
    v_inc_sd  = sd(v_incongruent),
    a_con_m   = mean(a_congruent),
    a_con_sd  = sd(a_congruent),
    a_inc_m   = mean(a_incongruent),
    a_inc_sd  = sd(a_incongruent),
    t0_con_m  = mean(t0_congruent),
    t0_con_sd = sd(t0_congruent),
    t0_inc_m  = mean(t0_incongruent),
    t0_inc_sd = sd(t0_incongruent),
  ) |>
  as.list()


# 3. Simulate DDM ground-truth --------------------------------------------

fm_ddm_gt <-
  tibble(
    id          = 1:2000,
    ntrials_con = ntrials_fm_con,
    ntrials_inc = ntrials_fm_inc,
    v_con_gt    = rnorm(2000, mean = sim_fm_ddm$v_con_m, sd = sim_fm_ddm$v_con_sd),
    v_inc_gt    = rnorm(2000, mean = sim_fm_ddm$v_inc_m, sd = sim_fm_ddm$v_inc_sd),
    a_con_gt    = rnorm(2000, mean = sim_fm_ddm$a_con_m, sd = sim_fm_ddm$a_con_sd),
    a_inc_gt    = rnorm(2000, mean = sim_fm_ddm$a_inc_m, sd = sim_fm_ddm$a_inc_sd),
    t0_con_gt   = rnorm(2000, mean = sim_fm_ddm$t0_con_m, sd = sim_fm_ddm$t0_con_sd),
    t0_inc_gt   = rnorm(2000, mean = sim_fm_ddm$t0_inc_m, sd = sim_fm_ddm$t0_inc_sd),
  ) |>
  # Remove very small non-decision times because the cause issues with fitting the model (and are biologically not very plausible)
  filter(a_con_gt > 0.2 & a_inc_gt > 0.2 & t0_con_gt > 0.01 & t0_inc_gt > 0.01) |>
  mutate(id = 1:n()) |>
  filter(id %in% sample(1:n(), size = nobs)) |>
  mutate(id = 1:nobs)



# 4. Simulate RTs/accuracies based on ground truth ------------------------

future::plan(future::multisession, workers = ncores)

fm_ddm_rt <- fm_ddm_gt %>%
  mutate(
    responses = furrr::future_pmap(., function(id, ntrials_con, ntrials_inc, v_con_gt, v_inc_gt, a_con_gt, a_inc_gt, t0_con_gt, t0_inc_gt) {

      bind_rows(
        # Simulate RTs/accuracy for congruent condition
        RWiener::rwiener(n=ntrials_con, alpha = a_con_gt, tau = t0_con_gt, beta = 0.5, delta = v_con_gt) |>
          as_tibble() |>
          mutate(
            con = 'congruent'
          ),
        # Simulate RTs/accuracy for incongruent condition
        RWiener::rwiener(n = ntrials_inc, alpha = a_inc_gt, tau = t0_inc_gt, beta = 0.5, delta = v_inc_gt) |>
          as_tibble() |>
          mutate(
            con = 'incongruent'
          )
      )
    }, .options = furrr::furrr_options(seed = TRUE))
  ) |>
  unnest(responses) |>
  select(subject = id, choice = resp, RT = q, condition = con, ntrials_con, ntrials_inc) |>
  mutate(choice = ifelse(choice == 'upper', 1, 0))

future::plan(future::sequential)


# 5. Fit HDDM to simulated RTs/accuracies ------------------------------------

# Drift rate and boundary separation separately for congruent and incongruent condition
# One fixed boundary separation across conditions
# Starting point fixed to 0.5

## Model Specification ----

model_fm <- "model {
  #likelihood function
  for (t in 1:nTrials) {
    y[t] ~ dwiener(alpha[condition[t], subject[t]],
                   tau[condition[t], subject[t]],
                   0.5,
                   delta[condition[t], subject[t]])
  }

  for (s in 1:nSubjects) {
    for (c in 1:nCon) {
      tau[c, s]  ~ dnorm(muTau[c], precTau) T(.0001, 1)
      delta[c, s] ~ dnorm(muDelta[c] , precDelta) T(-10, 10)
      alpha[c, s]  ~ dnorm(muAlpha[c], precAlpha) T(.05, 10)
    }

  }

  #priors
  for (c in 1:nCon){
    muTau[c] ~ dunif(.0001, 1)
    muDelta[c] ~ dunif(-10, 10)
    muAlpha[c] ~ dunif(.05, 10)
  }


  precAlpha  ~ dgamma(.001, .001)
  precTau ~ dgamma(.001, .001)
  precDelta ~ dgamma(.001, .001)
}"


## Prepare data ----

fm_ddm_rt_hddm <- fm_ddm_rt |>
  mutate(
    condition = ifelse(condition == "congruent", 1, 2),
    subject = rep(1:nobs, each = ntrials_fm_con + ntrials_fm_inc))

#Change error response times to negative for JAGS weiner module
y_fm <- round(ifelse(fm_ddm_rt_hddm$choice == 0, (fm_ddm_rt_hddm$RT*-1), fm_ddm_rt_hddm$RT),3)

condition_fm <- as.numeric(fm_ddm_rt_hddm$condition)



#Create numbers for looping purposes
nTrials_fm <- nrow(fm_ddm_rt_hddm)
nSubjects_fm <- length(unique(fm_ddm_rt_hddm$subject))
nCondition_fm <- max(condition_fm)



#Create a list of the data; this gets sent to JAGS
datalist_fm <- list(y = y_fm, condition = condition_fm,
                    subject = fm_ddm_rt_hddm$subject, nTrials = nTrials_fm,
                    nSubjects = nSubjects_fm, nCon = nCondition_fm)

## JAGS Specifications ----

#Need to tell JAGS where to start the samplers
#This function choses initial values randomly
initfunction <- function(chain){
  return(list(
    muAlpha = runif(2, .06, 9.95),
    muTau = runif(2, .01, .05),
    muDelta = runif(2, -9.9, 9.9),
    precAlpha = runif(1, .01, 100),
    precTau = runif(1, .01, 100),
    precDelta = runif(1, .01, 100),
    y = rep(NA, length(y_fm)),
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample.int(1e10, 1, replace = F)))
}

#Create list of parameters to be monitored
parameters <- c("alpha", "tau", "delta", "muAlpha",
                "muTau", "muDelta", "precAlpha", "precTau", "precDelta",
                "deviance")

nUseSteps = 1000 # Specify number of steps to run
nChains = 3# Specify number of chains to run (one per processor)

#Run the model in runjags
ddm_fm_mod <- run.jags(method = "parallel",
                       model = model_fm,
                       monitor = parameters,
                       data = datalist_fm,
                       inits = initfunction,
                       n.chains = nChains,
                       adapt = 1000, #how long the samplers "tune"
                       burnin = 2000, #how long of a burn in
                       sample = 2000,
                       thin = 1, #thin if high autocorrelation to avoid huge files
                       modules = c("wiener", "lecuyer"),
                       summarise = F,
                       plots = F)


# 6. Compare estimated and ground-truth DDM parameters

## Unpack Results ----

#Convert the runjags object to a coda format
codaSamples_fm <- as.mcmc.list(ddm_fm_mod)
mcmc_fm <- as.matrix(codaSamples_fm, chains = F) |>
  as_tibble()


# Traces for convergence checks
ddm_fm_traces_mod2 <- mcmc_fm |>
  select(starts_with("mu")) |>
  mutate(
    n = rep(1:2000, 3),
    chains = rep(1:3, each = 2000))


# Combine simulated and estimated DDM parameters
ddm_fm_data_mod2 <- mcmc_fm |>
  pivot_longer(everything(), names_to = "parameter", values_to = "estimated") |>
  group_by(parameter) |>
  summarise(estimated = mean(estimated, na.rm = T)) |>
  filter(str_detect(parameter, pattern = 'deviance|^mu|^prec', negate = T)) |>
  separate(col = parameter, into = c('parameter', 'id'), sep = "\\[") |>
  mutate(
    id = str_remove(id, pattern = "\\]$"),
    id = ifelse(parameter %in% c('delta', 'tau', 'alpha'),
                str_replace_all(id, "([0-9]*),([0-9]*)", "\\2,\\1"),
                id
    )
  ) |>
  separate(id, into = c('id', 'condition')) |>
  mutate(
    id = as.numeric(id),
    parameter = case_when(
      parameter == 'alpha' ~ 'a',
      parameter == 'tau' ~ 't0',
      parameter == 'delta' ~ 'v'
    ),
    condition = case_when(
      condition == 1 ~ "con",
      condition == 2 ~ "inc",
      is.na(condition) ~ NA_character_
    )
  ) |>
  left_join(
    fm_ddm_gt |>
      select(id, ends_with('_gt')) |>
      pivot_longer(-id, names_to = 'parameter', values_to = 'groundtruth') |>
      mutate(parameter = str_remove_all(parameter, "_gt$")) |>
      separate(parameter, into = c('parameter', 'condition'), sep = "_")
  )

ddm_fm_cor_mod2 <- ddm_fm_data_mod2 |>
  group_by(parameter) |>
  summarise(r = cor(estimated, groundtruth))


save(ddm_fm_traces_mod2, ddm_fm_data_mod2, ddm_fm_cor_mod2, file = "analysis_objects/ddm_fm_mod2_results.RData")


remove_DDM_files(path = "data/flanker")
