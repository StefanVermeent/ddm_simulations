

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(runjags)
library(RWiener)
library(here)
library(multitool)
library(glue)

nobs <- 720
ncores <- 4 # Number of cores used for parallel processing

set.seed(9474)

source('scripts/functions_DDM.R')

# I take the following six steps for the simulation analyses:

# 1. Based on observed mean RTs, accuracy rates, and correlations between task conditions, simulate 10,000 trials per participant
# 2. Estimate DDM parameters based on simulated trials (which will be highly reliable because of the large number of trials).
# 3. For each DDM parameter, calculate the mean and SD, and use these to simulate *ground-truth* DDM parameters for individual participants
# 4. Simulate RTs/accuracies based on *ground-truth* DDM parameters, with the same number of trials as the original study
# 5. Apply Hierarchical Bayesian DDM to the simulated RTs/accuracies
# 6. Compare estimated DDM parameters with simulated *ground-truth* DDM parameters, under different conditions (e.g., different model specifications)


# MSIT

ntrials_msit_con <- 24
ntrials_msit_inc <- 24

fmsit_tcon_m      <- 0.982
fmsit_tcon_sd     <- 0.041
fmsit_tinc_m      <- 0.820
fmsit_tinc_sd     <- 0.192

fmsit_tcon_rta_m  <- 0.832
fmsit_tcon_rta_sd <- 0.149
fmsit_tinc_rta_m  <- 1.430
fmsit_tinc_rta_sd <- 0.227


# 1. Simulate 1,000 trials per participant based on observed RTs --------

sim_msit_raw <-
  tibble(
    id   = 1:nobs
  ) %>%
  # Simulate mean/SD RTs per subject, mean accuracy, as well as rate and shape settings for the Gamma distribution (to simulate right-skewed RTs)
  mutate(
    faux::rnorm_multi(
      n = nobs,
      varnames = c("rt_con_m", "rt_con_sd", "rt_inc_m", "rt_inc_sd"),
      mu = c(fmsit_tcon_rta_m,  fmsit_tcon_rta_sd, fmsit_tinc_rta_m, fmsit_tinc_rta_sd),
      sd = c(fmsit_tcon_rta_sd, fmsit_tcon_rta_sd/2, fmsit_tinc_rta_sd, fmsit_tinc_rta_sd/2),
      r = c(   1,         .80,      .61,          .80,
               .80,         1,      .80,          .80,
               0.61,      .80,      1,            .80,
               .80,         .80,    .80,           1)
    ),
    across(c(rt_inc_sd, rt_con_sd), ~ifelse(. < 0, . * -1, .)),

    acc_con_m = rnorm(nobs, mean = fmsit_tcon_m, sd = fmsit_tcon_sd),
    acc_inc_m = rnorm(nobs, mean = fmsit_tinc_m, sd = fmsit_tinc_sd),

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


sim_msit_ddm <-
  sim_msit_raw |>
  select(trial_rts) |>
  unnest(trial_rts) |>
  select(id, rt, correct = acc, condition)

write_DDM_files(data = sim_msit_ddm, path = "data/msit", vars = c("rt", "correct", "condition"), task = "msit")


fast_dm_settings(task = "msit",
                 path = "data/msit",
                 model_version = "_mod2",
                 method = "ml",
                 depend = c("depends v condition", "depends a condition", "depends t0 condition"),
                 format = "TIME RESPONSE condition")

# Compute DDM parameters
execute_fast_dm(task = "msit", path = "data/msit", model_version = "_mod2")



# Read DDM results
sim_msit_ddm <- read_DDM(task = "msit", path = "data/msit", model_version = "_mod2") |>
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

msit_ddm_gt <-
  tibble(
    id          = 1:2000,
    ntrials_con = ntrials_msit_con,
    ntrials_inc = ntrials_msit_inc,
    v_con_gt    = rnorm(2000, mean = sim_msit_ddm$v_con_m,  sd = sim_msit_ddm$v_con_sd),
    v_inc_gt    = rnorm(2000, mean = sim_msit_ddm$v_inc_m,  sd = sim_msit_ddm$v_inc_sd),
    a_con_gt    = rnorm(2000, mean = sim_msit_ddm$a_con_m,  sd = sim_msit_ddm$a_con_sd),
    a_inc_gt    = rnorm(2000, mean = sim_msit_ddm$a_inc_m,  sd = sim_msit_ddm$a_inc_sd),
    t0_con_gt   = rnorm(2000, mean = sim_msit_ddm$t0_con_m, sd = sim_msit_ddm$t0_con_sd),
    t0_inc_gt   = rnorm(2000, mean = sim_msit_ddm$t0_inc_m, sd = sim_msit_ddm$t0_inc_sd),
  ) |>
  # Remove very small non-decision times and boundary separations because they cause issues with fitting the model (and are biologically not very plausible)
  filter(a_con_gt > 0.2 & a_inc_gt > 0.2 & t0_con_gt > 0.01 & t0_inc_gt > 0.01) |>
  mutate(id = 1:n()) |>
  filter(id %in% sample(1:n(), size = nobs))|>
  mutate(id = 1:nobs)



# 4. Simulate RTs/accuracies based on ground truth ------------------------

future::plan(future::multisession, workers = ncores)

msit_ddm_rt <- msit_ddm_gt %>%
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

model_msit <- "model {
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

msit_ddm_rt_hddm <- msit_ddm_rt |>
  mutate(
    condition = ifelse(condition == "congruent", 1, 2),
    subject = rep(1:nobs, each = ntrials_msit_con + ntrials_msit_inc))

#Change error response times to negative for JAGS weiner module
y_msit <- round(ifelse(msit_ddm_rt_hddm$choice == 0, (msit_ddm_rt_hddm$RT*-1), msit_ddm_rt_hddm$RT),3)

condition_msit <- as.numeric(msit_ddm_rt_hddm$condition)



#Create numbers for looping purposes
nTrials_msit <- nrow(msit_ddm_rt_hddm)
nSubjects_msit <- length(unique(msit_ddm_rt_hddm$subject))
nCondition_msit <- max(condition_msit)



#Create a list of the data; this gets sent to JAGS
datalist_msit <- list(y = y_msit, condition = condition_msit,
                      subject = msit_ddm_rt_hddm$subject, nTrials = nTrials_msit,
                      nSubjects = nSubjects_msit, nCon = nCondition_msit)

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
    y = rep(NA, length(y_msit)),
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
ddm_msit_mod <- run.jags(method = "parallel",
                         model = model_msit,
                         monitor = parameters,
                         data = datalist_msit,
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
codaSamples_msit <- as.mcmc.list(ddm_msit_mod)
mcmc_msit <- as.matrix(codaSamples_msit, chains = F) |>
  as_tibble()


# Traces for convergence checks
ddm_msit_traces <- mcmc_msit |>
  select(starts_with("mu")) |>
  mutate(
    n = rep(1:2000, 3),
    chains = rep(1:3, each = 2000))


# Combine simulated and estimated DDM parameters
ddm_msit_data <- mcmc_msit |>
  pivot_longer(everything(), names_to = "parameter", values_to = "estimated") |>
  group_by(parameter) |>
  summarise(estimated = mean(estimated, na.rm = T)) |>
  filter(str_detect(parameter, pattern = 'deviance|^mu|^prec', negate = T)) |>
  separate(col = parameter, into = c('parameter', 'id'), sep = "\\[") |>
  mutate(
    id = str_remove(id, pattern = "\\]$"),
    id = ifelse(parameter %in% c('delta', 'tau'),
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
    msit_ddm_gt |>
      select(id, ends_with('_gt')) |>
      pivot_longer(-id, names_to = 'parameter', values_to = 'groundtruth') |>
      mutate(parameter = str_remove_all(parameter, "_gt$")) |>
      separate(parameter, into = c('parameter', 'condition'), sep = "_")
  )

ddm_msit_cor <- ddm_msit_data |>
  group_by(parameter) |>
  summarise(r = cor(estimated, groundtruth))


save(ddm_msit_traces, ddm_msit_data, ddm_msit_cor, file = "analysis_objects/ddm_msit_results.RData")
