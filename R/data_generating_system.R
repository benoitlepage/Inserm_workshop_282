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
  low_par_edu <- rexpit(qlogis(0.7) + log(1.5) * sex) # low parent education
  
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
  occupation <- rexpit(qlogis(0.5) + log(1.3) * sex  + log(1.2) * low_par_edu + 
                   log(2.5) * edu + log(2) * phys) 
  occupation0 <- rexpit(qlogis(0.5) + log(1.3) * sex  + log(1.2) * low_par_edu + 
                          log(2.5) * edu0 + log(2) * phys0)  
  occupation1 <- rexpit(qlogis(0.5) + log(1.3) * sex  + log(1.2) * low_par_edu + 
                          log(2.5) * edu1 + log(2) * phys1) 
  
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
                      EY0_score = mean(score_tot0),
                      EY1_score = mean(score_tot1),
                      ATE_score = mean(score_tot1) - mean(score_tot0),
                      CDE0_score = mean(score10) - mean(score00),
                      CDE1_score = mean(score11) - mean(score01)))
  }
}

# ---------------------------------------------------------------------------- #
# 2) Simulate a data set of size N = 10000 ----
# ---------------------------------------------------------------------------- #
set.seed(1234)
df <- GenerateData.CDE(N = 10000, inter = rep(1, 10000), psi = FALSE)
summary(df)
write.csv2(df, "data/df.csv")
rm(df)
df <- read.csv2("data/df.csv")

# ---------------------------------------------------------------------------- #
# 3) Calculate "true" estimands  ----
# ---------------------------------------------------------------------------- #
set.seed(1234)
true <- GenerateData.CDE(N = 1e6, inter = rep(1, 1e6), psi = TRUE)
true
# $EY00_death
# [1] 0.097281
# 
# $EY01_death
# [1] 0.210761
# 
# $EY10_death
# [1] 0.164146
# 
# $EY11_death
# [1] 0.416408
# 
# $EY0_death
# [1] 0.153827
# 
# $EY1_death
# [1] 0.334302
# 
# $ATE_death      -> Average total effect death
# [1] 0.180475
# 
# $CDE0_death     -> CDE(m=0) death
# [1] 0.066865
# 
# $CDE1_death
# [1] 0.205647
# 
# $EY00_score
# [1] 48.78932
# 
# $EY01_score
# [1] 33.78664
# 
# $EY10_score
# [1] 36.99206
# 
# $EY11_score
# [1] 14.00505
# 
# $EY0_score
# [1] 41.70859
# 
# $EY1_score
# [1] 21.75026
# 
# $ATE_score
# [1] -19.95833   -> Average total effect score
# 
# $CDE0_score
# [1] -11.79726   -> CDE(m=0) score
# 
# $CDE1_score
# [1] -19.78158


# ---------------------------------------------------------------------------- #
# 4) test CMAverse  ----
# ---------------------------------------------------------------------------- #
library(CMAverse)

# ---------------------------------------------------------------------------- #
## 4.1) estimation par g-computation on the OR scale ----
# ---------------------------------------------------------------------------- #
set.seed(1234)
gformula_death_OR <- cmest(data = df, 
                           model = "gformula", # for parametric g-computation
                           outcome = "death", # outcome variable
                           exposure = "edu", # exposure variable
                           mediator = "smoking", # mediator
                           basec = c("sex",     # confounders
                                     "low_par_edu"), 
                           postc = c("phys", 
                                     "occupation"), # intermediate confounder (post-exposure)
                           EMint = TRUE, # exposures*mediator interaction
                           mreg = list("logistic"), # g(M=1|L1,A,L0)
                           yreg = "logistic",# Qbar.L2 = P(Y=1|M,L1,A,L0) # linear to get risk differences # loglinear for rate ratios
                           postcreg = list("logistic", "logistic"), # Qbar.L1 = P(L1=1|A,L0)
                           astar = 0,
                           a = 1,
                           mval = list(0), # do(M=0) to estimate CDE_m
                           estimation = "imputation", # if model= gformula
                           inference = "bootstrap",
                           boot.ci.type = "per", # for percentile, other option: "bca"
                           nboot = 2) # we should use a large number of bootstrap samples
summary(gformula_death_OR)
# Outcome regression:
# Call:
#   glm(formula = death ~ edu + smoking + edu * smoking + sex + low_par_edu + 
#         phys + occupation, family = binomial(), data = getCall(x$reg.output$yreg)$data, 
#       weights = getCall(x$reg.output$yreg)$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -3.02998    0.10428 -29.057  < 2e-16 ***
#   edu          0.65285    0.10332   6.319 2.63e-10 ***
#   smoking      1.03577    0.09909  10.453  < 2e-16 ***
#   sex          0.36821    0.04993   7.374 1.65e-13 ***
#   low_par_edu  0.44093    0.06151   7.169 7.58e-13 ***
#   phys        -0.10026    0.05037  -1.990  0.04654 *  
#   occupation   0.27883    0.06311   4.418 9.96e-06 ***
#   edu:smoking  0.32733    0.12155   2.693  0.00708 **    => interaction A*M
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 11536  on 9999  degrees of freedom
# Residual deviance: 10274  on 9992  degrees of freedom
# AIC: 10290
# Number of Fisher Scoring iterations: 5
# 
# 
# # Mediator regressions: 
# Call:
#   glm(formula = smoking ~ edu + sex + low_par_edu + phys + occupation, 
#       family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data, 
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -0.85996    0.06143 -13.998   <2e-16 ***
#   edu          0.70047    0.04402  15.911   <2e-16 ***
#   sex          0.62499    0.04365  14.318   <2e-16 ***
#   low_par_edu  0.41552    0.04772   8.707   <2e-16 ***
#   phys        -0.38703    0.04418  -8.761   <2e-16 ***
#   occupation   0.59325    0.04946  11.993   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 13464  on 9999  degrees of freedom
# Residual deviance: 12593  on 9994  degrees of freedom
# AIC: 12605
# Number of Fisher Scoring iterations: 4
# 
# 
# # Regressions for mediator-outcome confounders affected by the exposure: 
# Call:
#   glm(formula = phys ~ edu + sex + low_par_edu, family = binomial(), 
#       data = getCall(x$reg.output$postcreg[[1L]])$data, weights = getCall(x$reg.output$postcreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)  0.38815    0.04765   8.145 3.78e-16 ***
#   edu         -0.34545    0.04205  -8.215  < 2e-16 ***
#   sex          0.43714    0.04124  10.599  < 2e-16 ***
#   low_par_edu -0.19185    0.04689  -4.092 4.28e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 13729  on 9999  degrees of freedom
# Residual deviance: 13519  on 9996  degrees of freedom
# AIC: 13527
# Number of Fisher Scoring iterations: 4
# 
# Call:
#   glm(formula = occupation ~ edu + sex + low_par_edu, family = binomial(), 
#       data = getCall(x$reg.output$postcreg[[2L]])$data, weights = getCall(x$reg.output$postcreg[[2L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)  0.37498    0.05083   7.378 1.61e-13 ***
#   edu          0.85924    0.04734  18.150  < 2e-16 ***
#   sex          0.36173    0.04783   7.562 3.97e-14 ***
#   low_par_edu  0.09333    0.05217   1.789   0.0736 .          the model assumes independence with phys ?
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 8578.1  on 9999  degrees of freedom
# Residual deviance: 8037.8  on 9996  degrees of freedom
# AIC: 8045.8
# Number of Fisher Scoring iterations: 5
# 
# 
# # Effect decomposition on the odds ratio scale via the g-formula approach
# Direct counterfactual imputation estimation with 
# bootstrap standard errors, percentile confidence intervals and p-values 
#                     Estimate Std.error  95% CIL 95% CIU  P.val    
#   Rcde            2.0871979 0.2413521 2.0294568   2.354 <2e-16 *** CDE(m0), OR = 2.09
#   rRpnde          2.4201624 0.2522129 2.2122230   2.551 <2e-16 ***
#   rRtnde          2.5809151 0.2715732 2.2940236   2.659 <2e-16 ***
#   rRpnie          1.2171029 0.0134554 1.2002146   1.218 <2e-16 ***
#   rRtnie          1.2979457 0.0187438 1.2445946   1.270 <2e-16 ***
#   Rte             3.1412393 0.3617331 2.7533207   3.239 <2e-16 *** ATE, OR = 3.14
#   ERcde           0.5446765 0.1263706 0.4826422   0.652 <2e-16 ***
#   rERintref       0.8754859 0.1258423 0.7302174   0.899 <2e-16 ***
#   rERintmed       0.5039740 0.0960648 0.3412806   0.470 <2e-16 ***
#   rERpnie         0.2171029 0.0134554 0.2002183   0.218 <2e-16 ***
#   ERcde(prop)     0.2543744 0.0120060 0.2749897   0.291 <2e-16 ***
#   rERintref(prop) 0.4088688 0.0110529 0.4014920   0.416 <2e-16 ***
#   rERintmed(prop) 0.2353655 0.0114882 0.1944175   0.210 <2e-16 ***
#   rERpnie(prop)   0.1013913 0.0124414 0.0975362   0.114 <2e-16 ***
#   rpm             0.3367568 0.0009532 0.3073882   0.309 <2e-16 ***
#   rint            0.6442343 0.0004353 0.6107591   0.611 <2e-16 ***
#   rpe             0.7456256 0.0120060 0.7088802   0.725 <2e-16 ***


# ---------------------------------------------------------------------------- #
## 4.2) estimation par g-computation on the RD scale ----
# ---------------------------------------------------------------------------- #
set.seed(1234)
gformula_death_RD <- cmest(data = df, 
                           model = "gformula", # for parametric g-computation
                           outcome = "death", # outcome variable
                           exposure = "edu", # exposure variable
                           mediator = "smoking", # mediator
                           basec = c("sex",     # confounders
                                     "low_par_edu"), 
                           postc = c("phys", 
                                     "occupation"), # intermediate confounder (post-exposure)
                           EMint = TRUE, # exposures*mediator interaction
                           mreg = list("logistic"), # g(M=1|L1,A,L0)
                           yreg = "linear",# Qbar.L2 = P(Y=1|M,L1,A,L0) # linear for RD
                           postcreg = list("logistic", "logistic"), # Qbar.L1 = P(L1=1|A,L0)
                           astar = 0,
                           a = 1,
                           mval = list(0), # do(M=0) to estimate CDE_m
                           estimation = "imputation", # if model= gformula
                           inference = "bootstrap",
                           boot.ci.type = "per", # for percentile, other option: "bca"
                           nboot = 2) # we should use a large number of bootstrap samples
summary(gformula_death_RD)
#                  Estimate Std.error   95% CIL 95% CIU  P.val    
#   cde           0.0715898 0.0102149 0.0648080   0.079 <2e-16 ***  CDE(m=0) = 0.0716
#   rpnde         0.1415395 0.0124769 0.1284949   0.145 <2e-16 ***
#   rtnde         0.1735670 0.0166478 0.1519051   0.174 <2e-16 ***
#   rpnie         0.0242051 0.0002929 0.0222740   0.023 <2e-16 ***
#   rtnie         0.0562325 0.0038779 0.0460777   0.051 <2e-16 ***
#   te            0.1977720 0.0163549 0.1745726   0.197 <2e-16 ***  ATE = 0.1978
#   rintref       0.0699497 0.0022620 0.0636869   0.067 <2e-16 ***
#   rintmed       0.0320275 0.0041709 0.0234102   0.029 <2e-16 ***
#   cde(prop)     0.3619816 0.0210889 0.3711416   0.399 <2e-16 ***
#   rintref(prop) 0.3536883 0.0188551 0.3395706   0.365 <2e-16 ***
#   rintmed(prop) 0.1619413 0.0100660 0.1340541   0.148 <2e-16 ***
#   rpnie(prop)   0.1223887 0.0122998 0.1133772   0.130 <2e-16 ***
#   rpm           0.2843301 0.0022338 0.2609549   0.264 <2e-16 ***
#   rint          0.5156296 0.0087890 0.4871483   0.499 <2e-16 ***
#   rpe           0.6380184 0.0210889 0.6005255   0.629 <2e-16 ***

set.seed(1234)
gformula_score <- cmest(data = df, 
                        model = "gformula", # for parametric g-computation
                        outcome = "score", # outcome variable
                        exposure = "edu", # exposure variable
                        mediator = "smoking", # mediator
                        basec = c("sex",     # confounders
                                  "low_par_edu"), 
                        postc = c("phys", 
                                  "occupation"), # intermediate confounder (post-exposure)
                        EMint = TRUE, # exposures*mediator interaction
                        mreg = list("logistic"), # g(M=1|L1,A,L0)
                        yreg = "linear",# Qbar.L2 = P(Y=1|M,L1,A,L0) # linear to get risk differences # loglinear for rate ratios
                        postcreg = list("logistic", "logistic"), # Qbar.L1 = P(L1=1|A,L0)
                        astar = 0,
                        a = 1,
                        mval = list(0), # do(M=0) to estimate CDE_m
                        estimation = "imputation", # if model= gformula
                        inference = "bootstrap",
                        boot.ci.type = "per", # for percentile, other option: "bca"
                        nboot = 2) # we should use a large number of bootstrap samples
                           
summary(gformula_score)
# Direct counterfactual imputation estimation with 
# bootstrap standard errors, percentile confidence intervals and p-values 
#                   Estimate  Std.error    95% CIL 95% CIU  P.val    
#   cde           -14.585773   0.103593 -11.784720 -11.646 <2e-16 *** CDE(m=0) = -14.59
#   rpnde         -17.954120   0.161499 -15.354122 -15.137 <2e-16 ***
#   rtnde         -19.496366   0.091134 -16.716362 -16.594 <2e-16 ***
#   rpnie          -3.360983   0.468785  -3.268967  -2.639 <2e-16 ***
#   rtnie          -4.903230   0.539151  -4.725743  -4.001 <2e-16 ***
#   te            -22.857350   0.377652 -19.862890 -19.356 <2e-16 *** ATE = -22.86
#   rintref        -3.368346   0.265092  -3.708579  -3.352 <2e-16 ***
#   rintmed        -1.542246   0.070366  -1.456776  -1.362 <2e-16 ***
#   cde(prop)       0.638122   0.006224   0.593309   0.602 <2e-16 ***
#   rintref(prop)   0.147364   0.016989   0.168794   0.192 <2e-16 ***
#   rintmed(prop)   0.067473   0.002204   0.070378   0.073 <2e-16 ***
#   rpnie(prop)     0.147042   0.021009   0.136331   0.165 <2e-16 ***
#   rpm             0.214514   0.023213   0.206709   0.238 <2e-16 ***
#   rint            0.214836   0.014785   0.242134   0.262 <2e-16 ***
#   rpe             0.361878   0.006224   0.398329   0.407 <2e-16 ***

# ---------------------------------------------------------------------------- #
## 4.3) estimation par IPTW on the RD scale ----
# ---------------------------------------------------------------------------- #
set.seed(1234)
iptw_death_RD <- cmest(data = df, 
                       model = "msm", # using MSM estimated by IPTW
                       outcome = "death", # outcome variable
                       exposure = "edu", # exposure variable
                       mediator = "smoking", # mediator
                       basec = c("sex",     # confounders
                                 "low_par_edu"), 
                       postc = c("phys", 
                                 "occupation"), # intermediate confounder (post-exposure)
                       EMint = TRUE, # exposures*mediator interaction
                       ereg = "logistic", # exposure regression model g(A=1|L(0))
                       mreg = list("logistic"), # g(M=1|L1,A,L0)
                       yreg = "linear",# Qbar.L2 = P(Y=1|M,L1,A,L0) # linear to get risk differences # loglinear for rate ratios
                       postcreg = list("logistic", "logistic"), # Qbar.L1 = P(L1=1|A,L0)
                       wmnomreg = list("logistic"), #g(M=1|A) wgt nominator
                       wmdenomreg = list("logistic"), # g(M=1|L1,A,L(0)) wgt denominator
                       astar = 0,
                       a = 1,
                       mval = list(0), # do(M=0) to estimate CDE_m
                       estimation = "imputation", # if model= gformula
                       inference = "bootstrap",
                       boot.ci.type = "per", # for percentile, other option: "bca"
                       nboot = 2) # we should use a large number of bootstrap samples

summary(iptw_death_RD)
# # Outcome regression:
# Call:
#   glm(formula = death ~ edu + smoking + edu * smoking, family = gaussian(), 
#       data = getCall(x$reg.output$yreg)$data, weights = getCall(x$reg.output$yreg)$weights)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) 0.082933   0.008865   9.355  < 2e-16 ***
#   edu         0.070447   0.012775   5.514 3.59e-08 ***
#   smoking     0.119740   0.012960   9.239  < 2e-16 ***
#   edu:smoking 0.142064   0.017184   8.267  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for gaussian family taken to be 0.1690611)
# Null deviance: 1880.9  on 9999  degrees of freedom
# Residual deviance: 1689.9  on 9996  degrees of freedom
# AIC: 10969
# Number of Fisher Scoring iterations: 2
# 
# 
# # Mediator regressions: 
# Call:
#   glm(formula = smoking ~ edu, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data, 
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -0.11692    0.03149  -3.713 0.000205 ***
#   edu          0.78869    0.04174  18.895  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 13576  on 9999  degrees of freedom
# Residual deviance: 13214  on 9998  degrees of freedom
# AIC: 13197
# Number of Fisher Scoring iterations: 4
# 
# 
# # Mediator regressions for weighting (denominator): 
# Call:
#   glm(formula = smoking ~ edu + sex + low_par_edu + phys + occupation, 
#       family = binomial(), data = getCall(x$reg.output$wmdenomreg[[1L]])$data, 
#       weights = getCall(x$reg.output$wmdenomreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -0.85996    0.06143 -13.998   <2e-16 ***
#   edu          0.70047    0.04402  15.911   <2e-16 ***
#   sex          0.62499    0.04365  14.318   <2e-16 ***
#   low_par_edu  0.41552    0.04772   8.707   <2e-16 ***
#   phys        -0.38703    0.04418  -8.761   <2e-16 ***
#   occupation   0.59325    0.04946  11.993   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 13578  on 9999  degrees of freedom
# Residual deviance: 12685  on 9994  degrees of freedom
# AIC: 12697
# Number of Fisher Scoring iterations: 4
# 
# 
# # Mediator regressions for weighting (nominator): 
# Call:
#   glm(formula = smoking ~ edu, family = binomial(), data = getCall(x$reg.output$wmnomreg[[1L]])$data, 
#       weights = getCall(x$reg.output$wmnomreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -0.13313    0.03151  -4.225 2.39e-05 ***
#   edu          0.81446    0.04178  19.493  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 13578  on 9999  degrees of freedom
# Residual deviance: 13192  on 9998  degrees of freedom
# AIC: 13196
# Number of Fisher Scoring iterations: 4
# 
# 
# # Exposure regression for weighting: 
# Call:
#   glm(formula = edu ~ sex + low_par_edu, family = binomial(), data = getCall(x$reg.output$ereg)$data, 
#       weights = getCall(x$reg.output$ereg)$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)  0.02910    0.04186   0.695    0.487    
#   sex         -0.23178    0.04154  -5.580 2.41e-08 ***
#   low_par_edu  0.63793    0.04599  13.871  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 13497  on 9999  degrees of freedom
# Residual deviance: 13285  on 9997  degrees of freedom
# AIC: 13291
# Number of Fisher Scoring iterations: 4
# 
# 
# # Effect decomposition on the mean difference scale via the marginal structural model
# Direct counterfactual imputation estimation with 
# bootstrap standard errors, percentile confidence intervals and p-values 
# 
#                 Estimate Std.error  95% CIL 95% CIU  P.val    
#   cde           0.070447  0.007889 0.067239   0.078 <2e-16 *** CDE(m=0) = 0.0704
#   rpnde         0.137886  0.005815 0.137662   0.145 <2e-16 ***
#   rtnde         0.164690  0.003025 0.168127   0.172 <2e-16 ***
#   rpnie         0.022592  0.002277 0.020573   0.024 <2e-16 ***
#   rtnie         0.049395  0.000512 0.050350   0.051 <2e-16 ***
#   te            0.187281  0.005303 0.188700   0.196 <2e-16 *** ATE = 0.187
#   rintref       0.067440  0.002075 0.067636   0.070 <2e-16 ***
#   rintmed       0.026803  0.002789 0.026718   0.030 <2e-16 ***
#   cde(prop)     0.376155  0.030641 0.356284   0.397 <2e-16 ***
#   rintref(prop) 0.360098  0.020701 0.345418   0.373 <2e-16 ***
#   rintmed(prop) 0.143118  0.018617 0.136462   0.161 <2e-16 ***
#   rpnie(prop)   0.120629  0.008677 0.109013   0.121 <2e-16 ***
#   rpm           0.263747  0.009939 0.257133   0.270 <2e-16 ***
#   rint          0.503216  0.039318 0.481879   0.535 <2e-16 ***
#   rpe           0.623845  0.030641 0.602550   0.644 <2e-16 ***

# ---------------------------------------------------------------------------- #
# 5) test ltmle  ----
# ---------------------------------------------------------------------------- #
library(ltmle)
library(SuperLearner)
library(hal9001)

# ---------------------------------------------------------------------------- #
## 5.1) ATE death ----
# ---------------------------------------------------------------------------- #
Qform <- c(death = "Q.kplus1 ~ sex + low_par_edu + edu")
gform <- c("edu ~ sex + low_par_edu")
set.seed(1234)
ATE_tmle_death <- ltmle(data = subset(df, 
                                      select = c(sex, low_par_edu,
                                                 edu,
                                                 death)),
                        Anodes = "edu",
                        Ynodes = "death",
                        Qform = Qform,
                        gform = gform,
                        gbounds = c(0.01, 1),
                        abar = list(1,0), # vector of the counterfactual exposure 
                        SL.library = "glm", # glm without superlearner
                        variance.method = "tmle") # or ic
                                         
summary(ATE_tmle_death, estimator = "tmle")
# Additive Treatment Effect:
#   Parameter Estimate:  0.18649 
#    Estimated Std Err:  0.0082388 
#              p-value:  <2e-16 
#    95% Conf Interval: (0.17034, 0.20264)

# Odds Ratio:
#   Parameter Estimate:  2.9326
#  Est Std Err log(OR):  0.052901 
#              p-value:  <2e-16 
#    95% Conf Interval: (2.6437, 3.2529)

summary(ATE_tmle_death, estimator = "iptw")
# Additive Treatment Effect:
#   Parameter Estimate:  0.18657 
#    Estimated Std Err:  0.0083215 
#              p-value:  <2e-16 
#    95% Conf Interval: (0.17026, 0.20288) 

# Odds Ratio:
#   Parameter Estimate:  2.9343 
#  Est Std Err log(OR):  0.053393 
#              p-value:  <2e-16 
#    95% Conf Interval: (2.6428, 3.258) 

set.seed(1234)
ATE_gcomp_death <- ltmle(data = subset(df, 
                                       select = c(sex, low_par_edu,
                                                  edu,
                                                  death)),
                         Anodes = "edu",
                         Ynodes = "death",
                         Qform = Qform,
                         gform = gform,
                         gbounds = c(0.01, 1),
                         abar = list(1,0), # vector of the counterfactual exposure 
                         SL.library = "glm", # glm without superlearner
                         gcomp = TRUE, # to get gcomp estimation
                         variance.method = "ic") # not trustful with gcomp
summary(ATE_gcomp_death)                      
# Additive Treatment Effect:
#   Parameter Estimate:  0.18672 
#    Estimated Std Err:  0.0082388 
#              p-value:  <2e-16 
#    95% Conf Interval: (0.17057, 0.20286) 
# 
# Odds Ratio:
#   Parameter Estimate:  2.9379 
#  Est Std Err log(OR):  0.052951 
#              p-value:  <2e-16 
#    95% Conf Interval: (2.6483, 3.2592) 

# ---------------------------------------------------------------------------- #
## 5.2) ATE score ----
# ---------------------------------------------------------------------------- #
Qform <- c(score = "Q.kplus1 ~ sex + low_par_edu + edu")
gform <- c("edu ~ sex + low_par_edu")
SL_library <- list(Q=c("SL.mean","SL.glm","SL.lm"),
                   g=c("SL.mean","SL.glm","SL.lm"))
set.seed(1234)
ATE_tmle_score <- ltmle(data = subset(df, 
                                      select = c(sex, low_par_edu,
                                                 edu,
                                                 score)),
                        Anodes = "edu",
                        Ynodes = "score",
                        Qform = Qform,
                        gform = gform,
                        gbounds = c(0.01, 1),
                        abar = list(1,0), # vector of the counterfactual exposure 
                        SL.library = SL_library, # superlearner
                        variance.method = "ic") # faster than tmle
ATE_tmle_score$fit$g 
ATE_tmle_score$fit$Q 

summary(ATE_tmle_score, estimator = "tmle")
# Additive Treatment Effect:
# Parameter Estimate:  -19.76 
#  Estimated Std Err:  0.37737 
#            p-value:  <2e-16 
#  95% Conf Interval: (-20.5, -19.021) 

summary(ATE_tmle_score, estimator = "iptw")
# Additive Treatment Effect:
#   Parameter Estimate:  -19.768 
#    Estimated Std Err:  0.38332 
#              p-value:  <2e-16 
#    95% Conf Interval: (-20.52, -19.017) 



# ---------------------------------------------------------------------------- #
## 5.3) CDE death ----
# ---------------------------------------------------------------------------- #
Qform <- c(phys="Q.kplus1 ~ sex + low_par_edu + edu",
           death="Q.kplus1 ~ sex + low_par_edu + edu + phys + occupation + 
                  smoking + edu*smoking")
gform <- c("edu ~ sex + low_par_edu",
           "smoking ~ sex + low_par_edu + edu + phys + occupation")
set.seed(1234)
CDE_tmle_death <- ltmle(data = subset(df, 
                                      select = c(sex, low_par_edu,
                                                 edu,
                                                 phys, occupation,
                                                 smoking,
                                                 death)),
                        Anodes = c("edu", "smoking"),
                        Lnodes = c("phys", "occupation"),
                        Ynodes = "death",
                        Qform = Qform,
                        gform = gform,
                        gbounds = c(0.01, 1),
                        abar = list(c(1,0), # do(A=1, M=0)
                                    c(0,0)), # do(A=0, M=0)
                        SL.library = "glm", # glm without superlearner
                        variance.method = "ic") # faster than tmle

# check distribution of propensity scores
summary(CDE_tmle_death$cum.g.unbounded[,,1])
summary(CDE_tmle_death$cum.g.unbounded[,,2])

boxplot(CDE_tmle_death$cum.g.unbounded[,,1])
boxplot(CDE_tmle_death$cum.g.unbounded[,,2])

boxplot(CDE_tmle_death$cum.g[,,1])
boxplot(CDE_tmle_death$cum.g[,,2])
# no need for truncation in this example

summary(CDE_tmle_death, estimator = "tmle")
# Additive Treatment Effect:
# Parameter Estimate:  0.063685 
#  Estimated Std Err:  0.011188 
#            p-value:  1.2525e-08 
#  95% Conf Interval: (0.041758, 0.085613) 

# Odds Ratio:
#  Parameter Estimate:  1.7979 
# Est Std Err log(OR):  0.10418 
#             p-value:  1.7948e-08 
#   95% Conf Interval: (1.4658, 2.2051) 

summary(CDE_tmle_death, estimator = "iptw")
# Additive Treatment Effect:
# Parameter Estimate:  0.063619 
#  Estimated Std Err:  0.011248 
#            p-value:  1.5511e-08 
#  95% Conf Interval: (0.041572, 0.085665) 

# Odds Ratio:
#  Parameter Estimate:  1.7979 
# Est Std Err log(OR):  0.1048 
#             p-value:  2.1767e-08 
#   95% Conf Interval: (1.4641, 2.2079)

set.seed(1234)
CDE_gcomp_death <- ltmle(data = subset(df, 
                                       select = c(sex, low_par_edu,
                                                  edu,
                                                  phys, occupation,
                                                  smoking,
                                                  death)),
                         Anodes = c("edu", "smoking"),
                         Lnodes = c("phys", "occupation"),
                         Ynodes = "death",
                         Qform = Qform,
                         gform = gform,
                         gbounds = c(0.01, 1),
                         abar = list(c(1,0), # do(A=1, M=0)
                                     c(0,0)), # do(A=0, M=0)
                         SL.library = "glm", # glm without superlearner
                         gcomp = TRUE,
                         variance.method = "ic") # faster than tmle
                                      
summary(CDE_gcomp_death)
# Additive Treatment Effect:
#   Parameter Estimate:  0.066487 
#    Estimated Std Err:  0.011188 
#              p-value:  2.8052e-09 
#    95% Conf Interval: (0.044558, 0.088415) 

# Odds Ratio:
#   Parameter Estimate:  1.8397 
#  Est Std Err log(OR):  0.104 
#              p-value:  4.5891e-09 
#    95% Conf Interval: (1.5004, 2.2556) 

