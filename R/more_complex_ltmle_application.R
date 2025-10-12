### Application of the ltmle package on: 
###  - a longitudinal structure, with 2 waves
###  - exposure and mediators are 3-level catetorical variables
###  - survival with possible informative censoring


rm(list = ls())

# ---------------------------------------------------------------------------- #
# 1) data generating mechanism: ----
# ---------------------------------------------------------------------------- #

GenerateData.CDE <- function(N) { 
  # rexpit function
  rexpit <- function (x) rbinom(length(x), 1, plogis(x))
  
  # baseline confounders L0
  L0_1 <- rbinom(N, size = 1, prob = 0.45) 
  L0_2 <- rexpit(qlogis(0.6) + log(1.5) * L0_1) 
  
  # exposure A: treatment
  A1_1 <- rexpit(qlogis(0.2) + log(0.9) * L0_1  + log(1.5) * L0_2)
  A1_2 <- ifelse(A1_1 == 1,
                 0,
                 rexpit(qlogis(0.5) + log(0.8) * L0_1  + log(2) * L0_2))
  
  # censoring C1
  C1 <- rexpit(qlogis(0.03) + log(1.4) * L0_1  + log(1.4) * L0_2 + 
                 log(1.1) * A1_1 + log(1.4) * A1_2)
  
  # death Y1
  Y1 <- rexpit(qlogis(0.05) + log(1.2) * L0_1  + log(1.5) * L0_2 + 
                 log(2) * A1_1 + log(3) * A1_2)  
  
  ### intermediate counfounders L1
  L1 <- ifelse(Y1 == 1,
               NA,
               rnorm(N, mean = 50 + (5 * L0_1) + 
                       (-3 * L0_2) + (4 * A1_1) + (10 * A1_2),
                     sd = 15))

  # mediator M2: continue treatment
  M2_1 <- ifelse(Y1 == 1, 
                 NA, 
                 rexpit(qlogis(0.4) + log(0.9) * L0_1[Y1 == 0]  + log(1.5) * L0_2[Y1 == 0] + 
                          log(1.3) * A1_1[Y1 == 0] + log(1.6) * A1_2[Y1 == 0] + 
                          log(1.02) * L1[Y1 == 0]))
  M2_2 <- rep(NA, N)
  M2_2[Y1 == 1] <- NA
  M2_2[Y1 == 0 & M2_1 == 1] <- 0
  M2_2[Y1 == 0 & M2_1 == 0] <- rexpit(qlogis(0.4) + log(0.9) * L0_1[Y1 == 0 & M2_1 == 0]  + log(1.5) * L0_2[Y1 == 0 & M2_1 == 0] + 
                                        log(1.5) * A1_1[Y1 == 0 & M2_1 == 0] + log(2) * A1_2[Y1 == 0 & M2_1 == 0] + 
                                        log(1.03) * L1[Y1 == 0 & M2_1 == 0])
  
  # censoring C2
  C2 <- rep(NA, N)
  C2[Y1 == 1] <- NA
  C2[Y1 == 0 & C1 == 1] <- 1
  C2[Y1 == 0 & C1 == 0] <- rexpit(qlogis(0.03) + log(1.4) * L0_1[Y1 == 0 & C1 == 0]  + log(1.4) * L0_2[Y1 == 0 & C1 == 0] + 
                                    log(1.05) * A1_1[Y1 == 0 & C1 == 0] + log(1.2) * A1_2[Y1 == 0 & C1 == 0] + 
                                    log(1.01) * L1[Y1 == 0 & C1 == 0] + 
                                    log(1.1) * M2_1[Y1 == 0 & C1 == 0] + log(1.4) * M2_2[Y1 == 0 & C1 == 0])

  # death Y2
  Y2 <- rep(NA, N)
  Y2[Y1 == 1] <- 1
  Y2[Y1 == 0] <- rexpit(qlogis(0.05) + log(1.2) * L0_1[Y1 == 0]  + log(1.5) * L0_2[Y1 == 0] + 
                          log(1.5) * A1_1[Y1 == 0] + log(2) * A1_2[Y1 == 0] + 
                          log(1.01) * L1[Y1 == 0] + 
                          log(1.5) * M2_1[Y1 == 0] + log(2) * M2_2[Y1 == 0])  
                               
  df <- data.frame(L0_1 = L0_1, L0_2 = L0_2,
                   A1_1 = A1_1, A1_2 = A1_2, C1,
                   Y1 = Y1, L1 = L1,
                   M2_1 = M2_1, M2_2 = M2_2, C2,
                   Y2 = Y2)
  
  df$Y1 <- ifelse(df$C1 == 1, NA, df$Y1)
  df$L1 <- ifelse(df$C1 == 1, NA, df$L1)
  df$M2_1 <- ifelse(df$C1 == 1, NA, df$M2_1)
  df$M2_2 <- ifelse(df$C1 == 1, NA, df$M2_2)
  df$C2 <- ifelse(df$C1 == 1, NA, df$C2)
  df$Y2 <- ifelse(df$C2 == 1, NA, df$Y2)

  return(df)
}

# ---------------------------------------------------------------------------- #
# 2) generate 1 data frame ----
# ---------------------------------------------------------------------------- #
set.seed(1234)
df <- GenerateData.CDE(N = 10000)

View(df)
df |> 
  with(table(A1_1, A1_2))

df |> 
  subset(subset = (C1 == 1)) |>
  summary() # 519 censored after C1

df |> 
  subset(subset = (C1 == 0)) |>
  summary() # uncensored after C1

df |> 
  subset(subset = (C1 == 0 & Y1 == 1)) |>
  summary() # among uncensored at C1, 619 death at Y1

df |> 
  subset(subset = (C1 == 0 & Y1 == 0)) |>
  summary() # uncensored at C1 and alive at Y1

df |> 
  subset(subset = (C1 == 0 & Y1 == 0 & C2 == 1)) |>
  summary() # 737 censored at C2

df |> 
  subset(subset = (C1 == 0 & Y1 == 0)) |>
  with(table(M2_1, M2_2)) 

# OK


# ---------------------------------------------------------------------------- #
# 3) analyse using ltmle ----
# ---------------------------------------------------------------------------- #
library(ltmle)
df_ltmle <- df

# format censoring variables
df_ltmle$C1 <- BinaryToCensoring(is.censored = df_ltmle$C1)
df_ltmle$C2 <- BinaryToCensoring(is.censored = df_ltmle$C2)
# once a Ynode jumps to 1 (e.g. death), all subsequent Ynode values should be 1
df_ltmle$Y2 <- ifelse(df_ltmle$Y1 == 1, 1, df_ltmle$Y2) 

head(df_ltmle, 13)
#    L0_1 L0_2 A1_1 A1_2         C1 Y1       L1 M2_1 M2_2         C2 Y2
# 1     0    1    0    0 uncensored  0 84.37728    1    0 uncensored  0
# 2     1    1    1    0 uncensored  0 56.79832    0    1 uncensored  0
# 3     1    1    0    1 uncensored  0 68.84374    1    0 uncensored  0
# 4     1    0    0    0 uncensored  0 78.65583    1    0 uncensored  0
# 5     1    1    0    0 uncensored  0 61.33529    1    0 uncensored  0
# 6     1    1    0    0 uncensored  0 69.81963    0    1 uncensored  1
# 7     0    1    0    0 uncensored  0 42.79730    1    0 uncensored  1
# 8     0    1    0    0 uncensored  0 26.72749    1    0 uncensored  0
# 9     1    1    1    0 uncensored  0 51.65862    1    0 uncensored  0
# 10    0    1    0    0 uncensored  0 29.31751    1    0 uncensored  0
# 11    1    0    0    0 uncensored  1       NA   NA   NA       <NA>  1
# 12    0    0    0    0 uncensored  0 54.69811    0    1 uncensored  0
# 13    0    0    0    1 uncensored  0 36.30073    1    0 uncensored  0

SL.library <- c("SL.glm", "SL.mean")

# deal with multicategorical exposures and mediators
# we can incorporate deterministic knowledge
det.g <- function(data, current.node, nodes) {
  if (names(data)[current.node] != "A1_2" &  names(data)[current.node] != "M2_2") {
    return(NULL) # for other variables
  } else if (names(data)[current.node] == "A1_2") {
    is.deterministic <- data$A1_1 == 1 #if we're regressing A1_2, then: if A1_1=1 then P(A1_2 = 1) = 0
    return(list(is.deterministic = is.deterministic, prob1 = 0))
  } else if (names(data)[current.node] == "M2_2") {
    is.deterministic <- data$M2_1 == 1 #if we're regressing M2_2, then: if M2_1=1 then P(M2_2 = 1) = 0
    return(list(is.deterministic = is.deterministic, prob1 = 0))
  }
    else {
    stop("something went wrong!") #defensive programming
  }
}

CDE_A2vA0_M0 <- ltmle(data = df_ltmle,
                Anodes = c("A1_1", "A1_2", "M2_1", "M2_2"), # exposure and mediator
                Cnodes = c("C1", "C2"),
                Lnodes = c("L1"),
                Ynodes = c("Y1", "Y2"),
                survivalOutcome = TRUE,
                abar = list(c(0,1,0,0), # A1 = 2 # EY(A=2,M=0)
                            c(0,0,0,0)), # M2 = 0 # EY(A=0,M=0)
                deterministic.g.function = det.g,
                SL.library = SL.library,
                gcomp = FALSE,
                variance.method = "ic")
summary(CDE_A2vA0_M0)
# Treatment Estimate:
#   Parameter Estimate:  0.25043 
#    Estimated Std Err:  0.026775 
#              p-value:  <2e-16 
#    95% Conf Interval: (0.19795, 0.3029) 
# 
# Control Estimate:
#   Parameter Estimate:  0.21722 
#    Estimated Std Err:  0.027116 
#              p-value:  1.1395e-15 
#    95% Conf Interval: (0.16408, 0.27037) 
# 
# Additive Treatment Effect:
#   Parameter Estimate:  0.033203 
#    Estimated Std Err:  0.038102 
#              p-value:  0.38353 
#    95% Conf Interval: (-0.041476, 0.10788) 

CDE_A2vA0_M0$fit$g

CDE_A2vA0_M0$fit$Q


## Average total effect
## we do not need the mediators, but we need the censoring variables
## (which are considered as exposures)
ATE_A2vA0 <- ltmle(data = subset(df_ltmle, select = -c(M2_1, M2_2)),
                      Anodes = c("A1_1", "A1_2"), # exposure and mediator
                      Cnodes = c("C1", "C2"),
                      Lnodes = c("L1"),
                      Ynodes = c("Y1", "Y2"),
                      survivalOutcome = TRUE,
                      abar = list(c(0,1), # A1 = 2 # EY(A=2)
                                  c(0,0)), # M2 = 0 # EY(A=0)
                      deterministic.g.function = det.g,
                      SL.library = SL.library,
                      gcomp = FALSE,
                      variance.method = "ic")
summary(ATE_A2vA0)
# Treatment Estimate:
#   Parameter Estimate:  0.17078 
# Estimated Std Err:  0.025236 
# p-value:  1.3114e-11 
# 95% Conf Interval: (0.12132, 0.22025) 
# 
# Control Estimate:
#   Parameter Estimate:  0.1577 
# Estimated Std Err:  0.021452 
# p-value:  1.9593e-13 
# 95% Conf Interval: (0.11566, 0.19975) 
# 
# Additive Treatment Effect:
#   Parameter Estimate:  0.013082 
# Estimated Std Err:  0.033124 
# p-value:  0.69289 
# 95% Conf Interval: (-0.051839, 0.078003) 