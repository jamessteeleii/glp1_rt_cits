# Functions for targets ----

##### Functions for pre-registration power simulations ----

### Read in and prepare data to use for determining parameters in simulations ----

  # Read and prepare Karakakis et al (10.1016/j.metabol.2024.156113) data (baseline lean mass) to use for GLP1-RA SMD
  read_prep_karakakis_data <- function(file) {
    karakakis_baseline <- read_csv(file) |>
      slice_tail(n= 56) |>
      fill(`Author, year`) |>
      group_by(`Author, year`) |>
      mutate(
        study = cur_group_id()
      ) |>
      group_by(`Author, year`, Intervention) |>
      mutate(
        arm = cur_group_id()
      ) |>
      ungroup() |>
      select(study, arm, `No. of patients`, `Baseline lean mass`) |>
      rename(n = `No. of patients`,
             lean_mass = `Baseline lean mass`) |>
      separate(lean_mass, into = c("mean", "sd"), sep = c(" ")) |>
      mutate(
        mean = parse_number(mean),
        sd = parse_number(sd)
      )
    
    karakakis_baseline
    
  }
  
  # Read and prepare Benito et al. (10.3390/ijerph17041285) data for effects of RT "ONLY" on LMM
  read_prep_benito_data <- function(file) {
    benito_data <- read_csv(file) |>
      separate(study_arm, into = c("study", "arm"), sep = "_")|>
      group_by(study) |>
      mutate(
        study = cur_group_id()
      ) |>
      group_by(study, arm) |>
      mutate(
        arm = cur_group_id()
      ) |>
      filter(outcome == "ffm") |>
      mutate(
        duration_centre = duration/12
      ) |>
      filter(!is.na(duration))
    
    benito_data
    
  }
  
  # Read and prepare Murphy and Koehler (10.1111/sms.14075) data for effects of energy restriction during RT on LMM
  read_prep_murphy_koehler_data <- function(file) {
    murphy_koehler_data <- read_csv(file) |>
      group_by(study) |>
      mutate(
        study = cur_group_id()
      ) |>
      group_by(study, energy_deficit) |>
      mutate(
        arm = cur_group_id()
      ) |>
      ungroup() |>
      filter(outcome == "LM") |>
      mutate(
        vi = ((ci.ub - yi)/1.96)^2
      ) 
    
  }
  
  # Prepare data from Krupa et al. (10.3945/ajcn.116.137232) for effects of energy restriction alone
  # Note energy restriction in intervention group ~300 kcal/day 
  # Estimated roughly via simulation values of median split from summary data for high and low restriction groups
  prep_krupa_data <- function() {
    energy_deficit_12 <- tibble(
      energy_deficit = rnorm(1000, -378.1, 17.86 * sqrt(143))
    ) |>
      mutate(
        arm = case_when(
          energy_deficit < median(energy_deficit) ~ "high_cr",
          energy_deficit > median(energy_deficit) ~ "low_cr",
          energy_deficit == median(energy_deficit) ~ "median"
        )
      ) |>
      group_by(arm) |>
      summarise(
        energy_deficit = mean(energy_deficit)
      )
    
    energy_deficit_24 <-  tibble(
      energy_deficit = rnorm(1000,-298.7, 18.13 * sqrt(143))
    ) |>
      mutate(
        arm = case_when(
          energy_deficit < median(energy_deficit) ~ "high_cr",
          energy_deficit > median(energy_deficit) ~ "low_cr",
          energy_deficit == median(energy_deficit) ~ "median"
        )
      ) |>
      group_by(arm) |>
      summarise(
        energy_deficit = mean(energy_deficit)
      )
    
    krupa_data <- tibble(
      arm = c("con", "con", "low_cr", "low_cr", "high_cr", "high_cr"),
      energy_deficit = c(-38.3,-26.9,
                         energy_deficit_12$energy_deficit[2],
                         energy_deficit_24$energy_deficit[2],
                         energy_deficit_12$energy_deficit[1],
                         energy_deficit_24$energy_deficit[1]),
      time = rep(c(12,24), 3),
      n = c(75,75,57,57,58,58),
      pre_se = c(0.99,0.99,1.24,1.24,1.25,1.25),
      pre_sd = pre_se * sqrt(n),
      delta = c(-0.1,0,-1.76,-1.44,-2.32,-2.53),
      m2i = 0
    )
    
    krupa_data
  }
  
  # Read and prepare Beaver et al. (10.1002/oby.24172) data for examining assumption of interindividual response of GLP1-RA on lean mass
  read_prep_beaver_data <- function(file) {
    beaver_data <- read_csv(file) |>
      janitor::clean_names() |>
      fill(study) |>
      mutate(abs_mean_change = abs(mean_group_total_wt)) |>
      group_by(study, group) |>
      mutate(
        arm = cur_group_id()
      ) |>
      ungroup()|>
      select(study, arm, group, n_treat, sd_group_total_wt, abs_mean_change) |>
      pivot_wider(id_cols = study,
                  names_from = group,
                  values_from = c(n_treat, sd_group_total_wt, abs_mean_change))
    
  }

### Prepare parameters for simulations and check assumptions ----

  # Overall estimate for raw effect of GLP1-RAs vs placebo/control on lean mass from 10.1016/j.metabol.2024.156113
  get_raw_glp1_effect <- function() {
    
    raw_glp1_effect <- -0.86 
    
  }
  
  # Overall estimate for SMD effect of GLP1-RAs vs placebo/control on lean mass from 10.1016/j.metabol.2024.156113 
  # First meta-analyse baseline standard deviations for lean mass to then standardise overall raw effect estimate
  get_SMD_glp1_effect <- function(data, raw_glp1_effect) {
    
    karakakis_baseline <- data
    
    karakakis_baseline <- escalc(
      measure = "SDLN",
      mi = mean,
      sdi = sd,
      ni = n,
      data = karakakis_baseline
    )
    
    karakakis_baseline_meta_sd <- rma.mv(yi, vi,
                                         random = list(~ 1 | study, ~ 1 | arm),
                                         data = karakakis_baseline,
                                         method="REML", test="t"
    )
    
    baseline_sd <- exp(karakakis_baseline_meta_sd$beta[1])
    
    SMD_glp1_effect <- raw_glp1_effect / baseline_sd
    
  }
  
  # Overall estimate for raw effect of RT on lean mass from 10.3390/ijerph17041285
  # Note, median length is 8 weeks, though the majority of muscle gain occurs in first ~8-12 weeks with negligible diff at 24wks
  get_raw_RT_effect <- function(data) {
    
    benito_data <- data
    
    benito_data <- escalc(
      measure = "MC",
      m1i = post_mean,
      m2i = pre_mean,
      sd1i = post_sd,
      sd2i = pre_sd,
      ni = post_n,
      ri = 0.9,
      data = benito_data
    )
    
    benito_meta_raw <- rma.mv(yi, vi,
                              random = list(~ 1 | study, ~ 1 | arm),
                              data = benito_data,
                              method="REML", test="t"
    )
    
    raw_RT_effect <- benito_meta_raw$beta[1] 
    
  }
  
  # Overall estimate for SMD effect of RT on lean mass from 10.3390/ijerph17041285
  # Note, median length is 8 weeks, though the majority of muscle gain occurs in first ~8-12 weeks with negligible diff at 24wks
  get_SMD_RT_effect <- function(data) {
    
    benito_data <- data
    
    benito_data <- escalc(
      measure = "SMCR",
      m1i = post_mean,
      m2i = pre_mean,
      sd1i = post_sd,
      sd2i = pre_sd,
      ni = post_n,
      ri = 0.9,
      data = benito_data
    )
    
    benito_meta_SMD <- rma.mv(yi, vi,
                              random = list(~ 1 | study, ~ 1 | arm),
                              data = benito_data,
                              method="REML", test="t"
    )
    
    SMD_RT_effect <- benito_meta_SMD$beta[1] 
    
  }
  
  # Overall estimate for SMD effect of energy restriction during RT on lean mass from 10.1111/sms.14075
  get_SMD_RT_energy_restriction_effect <- function(data) {
    
    murphy_koehler_data <- data
    
    murphy_koehler_meta_smd <- rma.mv(yi, vi,
                                      random = list(~ energy_restriction | study, ~ 1 | arm),
                                      mods = ~ energy_restriction,
                                      data = murphy_koehler_data,
                                      method="REML", test="t"
    )
    
    SMD_RT_energy_restriction_effect <- murphy_koehler_meta_smd$beta[2]
    
  }
  
  # Overall estimate for SMD effect of energy restriction ONLY for low and high restriction on lean mass from 10.3945/ajcn.116.137232
  get_SMD_energy_restriction_effect <- function(data) {
    
    krupa_data <- data
    
    krupa_data <- escalc(
      measure = "SMCR",
      m1i = delta,
      m2i = m2i,
      sd1i = pre_sd,
      ni = n,
      ri = 0.7,
      data = krupa_data
    ) |>
      mutate(
        ci.lb = yi - sqrt(vi) * 1.96,
        ci.ub = yi + sqrt(vi) * 1.96
      ) 
    
  }
  
  # Check assumption of interindividual response variation to GLP1-RA on lean mass using 10.1002/oby.24172
  check_glp1_irv_assumption <- function(data){
    
    beaver_data <- data
    
    beaver_data <- escalc(
      measure = "VR",
      m1i = abs_mean_change_GLP,
      m2i = abs_mean_change_Placebo,
      sd1i = sd_group_total_wt_GLP,
      sd2i = sd_group_total_wt_Placebo,
      n1i = n_treat_GLP,
      n2i = n_treat_Placebo,
      data = beaver_data
    )
    
    beaver_meta_irv_VR <- rma(yi, vi,
                              data = beaver_data,
                              method="REML", test="t"
    )
    
    beaver_data <- escalc(
      measure = "CVR",
      m1i = abs_mean_change_GLP,
      m2i = abs_mean_change_Placebo,
      sd1i = sd_group_total_wt_GLP,
      sd2i = sd_group_total_wt_Placebo,
      n1i = n_treat_GLP,
      n2i = n_treat_Placebo,
      data = beaver_data
    )
    
    beaver_meta_irv_CVR <- rma(yi, vi,
                               data = beaver_data,
                               method="REML", test="t"
    )
    
    irv_glp1_effect <- tibble(
      effect = c("lnVR", "lnCVR"),
      estimate = c(beaver_meta_irv_VR$b, beaver_meta_irv_CVR$b),
      ci.lb = c(beaver_meta_irv_VR$ci.lb, beaver_meta_irv_CVR$ci.lb),
      ci.ub = c(beaver_meta_irv_VR$ci.ub, beaver_meta_irv_CVR$ci.ub)
    )
    
  }
  
  # Set intercepts i.e., baseline lean mass values
  # Kieser members ~55yrs old so values taken for lean mass for 50 yr olds from 10.1002/jcsm.12712
  # Note, assuming median = mean and a normal distribution to determine sigma
  set_intercepts_lean_mass <- function() {
    
    intercepts_lean_mass <- tibble(
      male_lean_mass = 66,
      female_lean_mass = 45.4,
      male_lean_mass_sigma = (50.8 - 66) / qnorm(0.03),
      female_lean_mass_sigma = (34.8 - 45.4) / qnorm(0.03),
    )
    
  }
  
  # Set error term based on approx measurement error taken from 
  set_measurement_error <- function() {
    # sigma <- 0.35 10.1123/ijsnem.2018-0283
    sigma <- 2.5 # from above but conservative just based on triangulating agreement between devices
    # sigma <- 0.6 # https://www.frontiersin.org/journals/nutrition/articles/10.3389/fnut.2024.1491931/full 
    
  }
  
### Simulate data for an interrupted time series with control ----
  
  # Define all simulation parameters
  # Note approx 60% male users and 40% female - from 10.1007/s00228-023-03539-8 - so weighted intercept and sigma
  define_simulation_parameters <- function(intercepts_lean_mass,
                                           raw_RT_effect,
                                           raw_glp1_effect,
                                           measurement_error) {
    
    simulation_parameters <- crossing(
      rep = 1:1000, # number of replicates
      participant_n = seq(10,200, by = 10), # range of participant N
      measurement_n = seq(5, 11, by = 2), # number of measurements total (before and after intervention introduction)
      b_intercept = weighted.mean(c(intercepts_lean_mass$male_lean_mass,intercepts_lean_mass$female_lean_mass), c(0.6,0.4)), # centred baseline intercept
      b_sex = intercepts_lean_mass$male_lean_mass - intercepts_lean_mass$female_lean_mass, # male-female weight diff in 
      b_rt = c(raw_RT_effect/12, (raw_glp1_effect/12)*-1), # effect of rt by week
      b_glp1 = raw_glp1_effect/12, # effect of glp1 by week
      b_period = 1, 
      b_intervention = 1,
      rho_AR = c(0,0.5,1),
      u_intercept = weighted.mean(c(intercepts_lean_mass$male_lean_mass_sigma,intercepts_lean_mass$female_lean_mass_sigma), c(0.6,0.4)), # weighted random intercept as roughly 60:40 split of males and females using glp1
      sigma = measurement_error
    )
    
  }

  
  
  simulate_power <- function(simulation_parameters) {
    
    sim <- function(participant_n = as.double(), measurement_n = as.double(),
                    b_intercept = as.double(), b_sex = as.double(), b_rt = as.double(), b_glp1 = as.double(), # fixed effects
                    b_period = 1, b_intervention = 1, # fixed effects 
                    u_intercept = as.double(), # random intercept
                    sigma = as.double(), # measurement error,
                    rho_AR = as.double(), # autocorrelation coefficient,
                    prop_missing = as.double(), # proportion of missing data
                    ... # helps the function work with pmap() below
    ) {
      
      weeks <- measurement_n*12
      
      # set up data structure
      data <- add_random(participant = participant_n) %>%
        # add within participant time
        add_within("participant", time = seq(from=0,to=weeks, by=12)) |>
        mutate(period = case_when(
          time >= weeks/2 ~ 1,
          .default = 0
        )) |>
        # add and code categorical variables
        add_between("participant", intervention = c(0, 1)) %>% # con = 0, glp1 = 1
        add_between("participant", sex = c(-0.5, 0.5), .prob = c(0.4, 0.6)) %>% # female = -0.5, male  = 0.5
        # add random effects 
        add_ranef("participant", u_intercept = u_intercept) %>%
        add_ranef(sigma = sigma) %>%
        # calculate DV
        mutate(ffm = b_intercept + u_intercept + 
                 (b_sex*sex) + 
                 (b_rt*time) + (b_glp1*time*period*intervention) +
                 (b_period*period*intervention) +
                 sigma) |>
        group_by(participant) |>
        mutate(
          sigma_AR = case_when(
            time == 0 ~ sigma,
            .default = sigma + rho_AR*lag(sigma)
          ),
          ffm_AR = b_intercept + u_intercept + 
            (b_sex*sex) + 
            (b_rt*time) + (b_glp1*time*period*intervention) +
            (b_period*period*intervention) +
            sigma_AR
        )
      
      
      model <- lme(ffm_AR ~ sex + time + intervention + time:intervention + period + period:time + period:intervention + period:intervention:time,
                   random = ~ 1 | participant, data = data, correlation = corARMA(form = ~ 1 | participant, p = 1, q = 1))
      
      test_glp1_slope <- avg_slopes(
        model,
        by = c("intervention","period"),
        variables = "time",
        equivalence = c(-0.07, 0.13)
      ) |>
        slice_tail()
      
      test_glp1_slope <- as.data.frame(test_glp1_slope) |>
        bind_cols(
          tibble(
            glp1_effect_estimate = intervals(summary(model))$fixed[9,2],
            glp1_effect_ci.lb = intervals(summary(model))$fixed[9,1],
            glp1_effect_ci.ub = intervals(summary(model))$fixed[9,3],
            glp1_effect_p.value = summary(model)$tTable[9,5],
            RT_effect_estimate = intervals(summary(model))$fixed[3,2],
            RT_effect_ci.lb = intervals(summary(model))$fixed[3,1],
            RT_effect_ci.ub = intervals(summary(model))$fixed[3,3],
            RT_effect_p.value = summary(model)$tTable[3,5]
            
          )
        )
      
      rm(data)
      rm(model)
      
      gc()
      
      return(test_glp1_slope)
      
    }
    
    safe_sim <- safely(sim)
    
    # handlers(global = TRUE)  # allow progress bars
    # 
    # with_progress({
    #   
    #   p <- progressor(steps = nrow(simulation_parameters))
    #   
    #   safe_simulations <- simulation_parameters %>%
    #     mutate(analysis = pmap(., function(...) {
    #       p()
    #       safe_sim(...)
    #     }))
    # })
    
    safe_simulations <- simulation_parameters %>%
      mutate(analysis = pmap(., safe_sim))
    
    
    safe_simulations_clean <- safe_simulations |>
      filter(map_lgl(analysis, ~ !is.null(.x$result) && nrow(.x$result) > 0)) %>%
      mutate(analysis = map(analysis, "result")) %>%
      unnest(analysis)

    return(safe_simulations_clean)
  }

    
  
  