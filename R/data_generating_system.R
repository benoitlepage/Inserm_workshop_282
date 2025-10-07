# 
# sex         L0.1 = sex 
# low_par_edu L0.2 = low parent educational level
# edu         A = exposure of interest
# phys        L1.1 = physical activity
# occupation  L1.2 = occupation
# smoking     M = mediator of interest
# death       Y binary
# score       Y continuous

rm(list = ls())

# ---------------------------------------------------------------------------- #
# 1) data generating mechanism: ----
# ---------------------------------------------------------------------------- #
### function to create data sets and estimate "true" parameters via Monte Carlo simulations
### when {edu and smoking} are specified

GenerateData.CDE <- function(N, 
                             inter = rep(1, N), # presence of A*M interaction
                             psi = FALSE) { # FALSE = simulate data.fame only
  ### rexpit function
  rexpit <- function (x) rbinom(length(x), 1, plogis(x))
  
  ### baseline confounders L0
  sex <- rbinom(N, size = 1, prob = 0.45) # 0 = women; 1 = men
  low_par_edu <- rexpit(qlogis(0.7) + log(1.5) * sex ) # low parent education
  
  ### exposure A: low educational level = 1
  edu <- rexpit(qlogis(0.5) + log(0.8) * sex  + log(2) * low_par_edu)
  edu0 <- rep(0, N)
  edu1 <- rep(1, N)
  
  ### intermediate counfounders L1
  # physical activity: 1 = yes ; 0 = no
  phys <- rexpit(qlogis(0.6) + log(1.5) * sex  + log(0.8) * low_par_edu + 
                   log(0.7) * edu)
  phys0 <- rexpit(qlogis(0.6) + log(1.5) * sex  + log(0.8) * low_par_edu + 
                    log(0.7) * edu0)
  phys1 <- rexpit(qlogis(0.6) + log(1.5) * sex  + log(0.8) * low_par_edu + 
                    log(0.7) * edu1)
  
  # occupation: 1 = manual; 0 = non-manual
  occupation <- rexpit(qlogis(0.6) + log(1.5) * sex  + log(1.5) * low_par_edu + 
                   log(3) * edu + log(2) * phys)
  occupation0 <- rexpit(qlogis(0.6) + log(1.5) * sex  + log(1.5) * low_par_edu + 
                          log(3) * edu0 + log(2) * phys0)  
  occupation1 <- rexpit(qlogis(0.6) + log(1.5) * sex  + log(1.5) * low_par_edu + 
                          log(3) * edu1 + log(2) * phys1) 
  
  ### mediator
  smoking <- rexpit(qlogis(0.3) + log(1.8) * sex  + log(1.5) * low_par_edu + 
                      log(2) * edu + log(0.7) * phys + log(1.8) * occupation)
  smoking0 <- rep(0, N)
  smoking1 <- rep(1, N)
  smoking_tot0 <- rexpit(qlogis(0.3) + log(1.8) * sex  + log(1.5) * low_par_edu + 
                           log(2) * edu0 + log(0.7) * phys0 + log(1.8) * occupation0)
  smoking_tot1 <- rexpit(qlogis(0.3) + log(1.8) * sex  + log(1.5) * low_par_edu + 
                           log(2) * edu1 + log(0.7) * phys1 + log(1.8) * occupation1)
                     
  
  ### outcomes
  death <- rexpit(qlogis(0.05) + log(1.5) * sex  + log(1.6) * low_par_edu + 
                    log(1.7) * edu + log(0.8) * phys + log(1.5) * occupation + 
                    log(2.5) * smoking + log(1.5) * edu * smoking * inter)
  death00 <- rexpit(qlogis(0.05) + log(1.5) * sex  + log(1.6) * low_par_edu + 
                      log(1.7) * edu0 + log(0.8) * phys0 + log(1.5) * occupation0 + 
                      log(2.5) * smoking0 + log(1.5) * edu0 * smoking0 * inter)
  death01 <- rexpit(qlogis(0.05) + log(1.5) * sex  + log(1.6) * low_par_edu + 
                      log(1.7) * edu0 + log(0.8) * phys0 + log(1.5) * occupation0 + 
                      log(2.5) * smoking1 + log(1.5) * edu0 * smoking1 * inter)
  death10 <- rexpit(qlogis(0.05) + log(1.5) * sex  + log(1.6) * low_par_edu + 
                      log(1.7) * edu1 + log(0.8) * phys1 + log(1.5) * occupation1 + 
                      log(2.5) * smoking0 + log(1.5) * edu1 * smoking0 * inter)
  death11 <- rexpit(qlogis(0.05) + log(1.5) * sex  + log(1.6) * low_par_edu + 
                      log(1.7) * edu1 + log(0.8) * phys1 + log(1.5) * occupation1 + 
                      log(2.5) * smoking1 + log(1.5) * edu1 * smoking1 * inter)
  death_tot0 <- rexpit(qlogis(0.05) + log(1.5) * sex  + log(1.6) * low_par_edu + 
                         log(1.7) * edu0 + log(0.8) * phys0 + log(1.5) * occupation0 + 
                         log(2.5) * smoking_tot0 + log(1.5) * edu0 * smoking_tot0 * inter)
  death_tot1 <- rexpit(qlogis(0.05) + log(1.5) * sex  + log(1.6) * low_par_edu + 
                         log(1.7) * edu1 + log(0.8) * phys1 + log(1.5) * occupation1 + 
                         log(2.5) * smoking_tot1 + log(1.5) * edu1 * smoking_tot1 * inter)

  score <- rnorm(N, mean = 50 + 5 * sex  -5 * low_par_edu + 
                   -10 * edu + 8 * phys -7 * occupation + 
                   -15 * smoking + -8 * edu * smoking * inter, 
                 sd = 15)
  score00 <- rnorm(N, mean = 50 + 5 * sex  -5 * low_par_edu + 
                     -10 * edu0 + 8 * phys0 -7 * occupation0 + 
                     -15 * smoking0 + -8 * edu0 * smoking0 * inter,
                   sd = 15)
  score01 <- rnorm(N, mean = 50 + 5 * sex  -5 * low_par_edu + 
                     -10 * edu0 + 8 * phys0 -7 * occupation0 + 
                     -15 * smoking1 + -8 * edu0 * smoking1 * inter,
                   sd = 15)  
  score10 <- rnorm(N, mean = 50 + 5 * sex  -5 * low_par_edu + 
                     -10 * edu1 + 8 * phys1 -7 * occupation1 + 
                     -15 * smoking0 + -8 * edu1 * smoking0 * inter,
                   sd = 15)   
  score11 <- rnorm(N, mean = 50 + 5 * sex  -5 * low_par_edu + 
                     -10 * edu1 + 8 * phys1 -7 * occupation1 + 
                     -15 * smoking1 + -8 * edu1 * smoking1 * inter, 
                   sd = 15)   
  score_tot0 <- rnorm(N, mean = 50 + 5 * sex  -5 * low_par_edu + 
                        -10 * edu0 + 8 * phys0 -7 * occupation0 + 
                        -15 * smoking_tot0 + -8 * edu0 * smoking_tot0 * inter, 
                      sd = 15)
  score_tot1 <- rnorm(N, mean = 50 + 5 * sex  -5 * low_par_edu + 
                        -10 * edu1 + 8 * phys1 -7 * occupation1 + 
                        -15 * smoking_tot1 + -8 * edu1 * smoking_tot1 * inter, 
                      sd = 15)
                   
  
  if (psi == FALSE) {
    return(data.sim = data.frame(subjid = 1:N,
                                 sex = sex, 
                                 low_par_edu = low_par_edu,
                                 edu = edu,
                                 phys = phys,
                                 occupation = occupation,
                                 smoking = smoking,
                                 death = death,
                                 score = score))
  } else {
    return(Psi = list(EY00_death = mean(death00),
                      EY01_death = mean(death01),
                      EY10_death = mean(death10),
                      EY11_death = mean(death11),
                      EY0_death = mean(death_tot0), 
                      EY1_death = mean(death_tot1),
                      ATE_death = mean(death_tot1) - mean(death_tot0),
                      CDE0_death = mean(death10) - mean(death00),
                      CDE1_death = mean(death11) - mean(death01),
                      EY00_score = mean(score00),
                      EY01_score = mean(score01),
                      EY10_score = mean(score10),
                      EY11_score = mean(score11),
                      CDE0_score = mean(score10) - mean(score00),
                      CDE1_score = mean(score11) - mean(score01),
                      EY0_score = mean(score_tot0),
                      EY1_score = mean(score_tot1),
                      ATE_score = mean(score_tot1) - mean(score_tot0)))
  }
}

# ---------------------------------------------------------------------------- #
# 2) Simulate a data set of size N = 10000 ----
# ---------------------------------------------------------------------------- #
set.seed(1234)
df <- GenerateData.CDE(N = 10000, inter = rep(1, 10000), psi = FALSE)
summary(df)

# ---------------------------------------------------------------------------- #
# 3) Calculate "true" estimands  ----
# ---------------------------------------------------------------------------- #
set.seed(1234)
true <- GenerateData.CDE(N = 1e6, inter = rep(1, 1e6), psi = TRUE)
true
# $EY00_death
# [1] 0.10123
# 
# $EY01_death
# [1] 0.218588
# 
# $EY10_death
# [1] 0.168742
# 
# $EY11_death
# [1] 0.424843
# 
# $EY0_death
# [1] 0.161479
# 
# $EY1_death
# [1] 0.344455
# 
# $ATE_death
# [1] 0.182976
# 
# $CDE0_death
# [1] 0.067512
# 
# $CDE1_death
# [1] 0.206255
# 
# $EY00_score
# [1] 47.94446
# 
# $EY01_score
# [1] 32.94178
# 
# $EY10_score
# [1] 36.35606
# 
# $EY11_score
# [1] 13.36906
# 
# $CDE0_score
# [1] -11.5884
# 
# $CDE1_score
# [1] -19.57272
# 
# $EY0_score
# [1] 40.60899
# 
# $EY1_score
# [1] 20.84092
# 
# $ATE_score
# [1] -19.76807
