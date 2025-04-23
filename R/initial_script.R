library(tidyverse)
library(patchwork)
library(metafor)
library(emmeans)
library(faux)

# cHECKLIST https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-021-01235-8
# MISSINGNESS https://pubmed.ncbi.nlm.nih.gov/34326669/


##### Setting parameters for simulations ----

# Baseline lean mass data - taken from studies included in 10.1016/j.metabol.2024.156113 
karakakis_baseline <- read_csv("data/karakakis_baseline_study_data.csv") |>
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

karakakis_baseline <- escalc(
  measure = "MN",
  mi = mean,
  sdi = sd,
  ni = n,
  data = karakakis_baseline
)

karakakis_baseline_meta_mean <- rma.mv(yi, vi,
                                  random = list(~ 1 | study, ~ 1 | arm),
                                  data = karakakis_baseline,
                                  method="REML", test="t"
                                  )

baseline_mean <- karakakis_baseline_meta_mean$beta[1]

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

# Mean effects  - overall estimate for lean mass from GLP1-RAs vs placebo/control from 10.1016/j.metabol.2024.156113
# Converted to Mean effects  using meta-analysed baseline standard deviation
# Duration of intervention ~ 24 weeks

mean_effect_glp1 <- -0.86

standardised_mean_effect_glp1 <- mean_effect_glp1 / baseline_sd

# Mean effects of RT ONLY on lean mass - taken from studies included in 10.3390/ijerph17041285
# Duration of intervention ~ 10 weeks

benito_data <- read_csv("data/benito_study_data.csv") |>
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

benito_data <- escalc(
  measure = "MC",
  m1i = post_mean,
  m2i = pre_mean,
  sd1i = post_sd,
  sd2i = pre_sd,
  ni = post_n,
  ri = 0.7,
  data = benito_data
)

benito_meta_raw <- rma.mv(yi, vi,
                                       random = list(~ duration_centre | study, ~ duration_centre | arm),
                          mods = ~ duration_centre,
                                       data = benito_data,
                                       method="REML", test="t"
)

mean_effect_rt <- benito_meta_raw$beta[1]

benito_data <- escalc(
  measure = "SMCR",
  m1i = post_mean,
  m2i = pre_mean,
  sd1i = post_sd,
  sd2i = pre_sd,
  ni = post_n,
  ri = 0.7,
  data = benito_data
)

benito_meta_smd<- rma.mv(yi, vi,
                          random = list(~ 1 | study, ~ 1 | arm),
                          data = benito_data,
                          method="REML", test="t"
)

standardised_mean_effect_rt <- benito_meta_smd$beta[1]


# Mean effects  - overall estimate for lean mass for energy restriction during RT from 10.1111/sms.14075
# Mean energy restriction ~600 kcal/day
# Duration of intervention ~13 weeks

murphy_koehler_data <- read_csv("data/murphy_koehler_study_data.csv") |>
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

murphy_koehler_meta_smd <- rma.mv(yi, vi,
                         random = list(~ energy_restriction | study, ~ 1 | arm),
                         mods = ~ energy_restriction,
                         data = murphy_koehler_data,
                         method="REML", test="t"
)

standardised_mean_effect_energy_restriction_rt <- murphy_koehler_meta_smd$beta[2]


# Mean effects  - overall estimate for lean mass for energy restriction ONLY from 10.3945/ajcn.116.137232
# Energy restriction ~300 kcal/day - estimate roughly via simulation values of median split from summary data
# Duration of intervention 52 and 104 weeks
# to compare with RT plus energy restriction to see if any interction or can assume additive effects

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

# Note we don't assume any interindividual slope variation as little evidence for IRV from 10.1002/oby.24172 data

beaver_data <- read_csv("data/beaver_study_data.csv") |>
  janitor::clean_names() |>
  fill(study) |>
  mutate(abs_mean_change = abs(mean_group_total_wt)) |>
  group_by(study, group) |>
  mutate(
    arm = cur_group_id()
  ) |>
  ungroup()

beaver_data <- escalc(
  measure = "SDLN",
  sdi = sd_group_total_wt,
  ni = n_treat,
  data = beaver_data
)

beaver_meta_irv <- rma.mv(yi, vi,
                          random = list(~ 1 | study, ~ 1 | arm),
                          mods = ~ log(abs_mean_change) + group,
                          data = beaver_data,
                          method="REML", test="t"
                          )


# Simulate data for an interrupted time series with control ----

# approx 60% male users and 40% female - from 10.1007/s00228-023-03539-8

# based on our client weight data and estimates of FFM% for ~55yr olds from 10.1002/jcsm.12712
male_ffm <- 86.8 * 0.66
female_ffm <- 72 * 0.45
male_ffm_sigma <- 17 * 0.66
female_ffm_sigma <- 17.6 * 0.45

# we'll assume that all the intervention effect estimates above are similar over a ~12wk period even if estimated from longer
# for example the CALERIE study was 12/24 months, but showed similar degree of weight loss in total to that seen in 10.1002/fsn3.4442 over ~3 months on similar % deficit
# We'll assume linear trends and treat all intervention slopes as being by week i.e., effects assumed to be over 12wks so slope is effect divided by 12

# define parameters
b_intercept <- mean(c(male_ffm,female_ffm))
b_sex <- male_ffm - female_ffm # male-female weight diff in 
b_rt <- mean_effect_rt/12 # effect of rt by week
b_glp1 <- mean_effect_glp1/12 # effect of glp1 by week
b_period <- 1 
b_intervention <- 1
rho_AR <- 0.5
u_intercept <- (0.6 * male_ffm_sigma) + (0.4 * female_ffm_sigma) # weighted random intercept as roughly 60:40 split of males and females using glp1
sigma <- 0.35 # rough estimate of measurement error taken from 10.1123/ijsnem.2018-0283

# set up data structure
data <- add_random(subj = 100) %>%
  # add within participant time
  add_within("subj", time = c(0,6,12,18,24,30)) |>
  mutate(period = case_when(
    time > 12 ~ 1,
    .default = 0
  )) |>
  # add and code categorical variables
  add_between("subj", intervention = c(0, 1)) %>% # con = 0, glp1 = 1
  add_between("subj", sex = c(-0.5, 0.5), .prob = c(0.4, 0.6)) %>% # female = -0.5, male  = 0.5
  # add random effects 
  add_ranef("subj", u_intercept = u_intercept) %>%
  add_ranef(sigma = sigma) %>%
  # calculate DV
  mutate(ffm = b_intercept + u_intercept + 
           (b_sex*sex) + 
           (b_rt*time) + (b_glp1*time*period*intervention) +
           (b_period*period*intervention) +
           sigma) |>
  group_by(subj) |>
  mutate(
    sigma = case_when(
      time == 0 ~ sigma,
      .default = sigma + rho_AR*lag(sigma)
    ),
    ffm_AR = b_intercept + u_intercept + 
      (b_sex*sex) + 
      (b_rt*time) + (b_glp1*time*period*intervention) +
      (b_period*period*intervention) +
      sigma
  )
         


data |>
  ggplot(aes(x=factor(time), y=ffm,
             color = factor(intervention),
             group = interaction(factor(intervention), factor(period)))) +
  # geom_point() +
  stat_smooth(data = data |> filter(period == 0),
              method = "lm", fullrange = TRUE, linetype = "dashed", se = FALSE) +
  geom_smooth(method = "lm")

library(nlme)

model <- lme(ffm_AR ~ sex + time + intervention + time:intervention + period + period:time + period:intervention + period:intervention:time,
                   random = ~ 1 | subj, data = data, correlation = corARMA(form = ~ 1 | subj, p = 1, q = 1))


summary(model)


plot(nlme::ACF(model, resType="normalized"), alpha = 0.05)

library(marginaleffects)

avg_slopes(
  model,
  # newdata = datagrid(
  #   sex = 0,
  #   time = c(18,24,30),
  #   intervention = 1,
  #   period = 1,
  #   # subj = NA
  # ),
  by = c("intervention","period"),
  variables = "time",
  equivalence = c(-0.07, 0.13), 
  re.form = NA
)

# TO VARY IN SIMULATIONS
# NO. REVIEW TIME POINTS (EVERY ~12 WEEKS, VARY AROUND THAT EXACT TIME E.G. BY ~2 WEEKS) - ADD MISSINGNESS E.G., MISS REVIEW
# RANDOM SLOPE ASSUMPTION FOR TIME (I.E,. RT) AND ALSO GLP1 - NOT LIKELY, BUT BE CONSERVATIVE
# the review numbers determine duration of intervention too - check Benito data based on duration of 24 weeks

