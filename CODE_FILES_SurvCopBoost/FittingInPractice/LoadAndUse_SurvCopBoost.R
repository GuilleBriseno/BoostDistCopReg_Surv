setwd("../../BoostingCopulas_Papers/SurvCopBoost_CODE_FILES/")

library("VineCopula")
library("VGAM")
library("gamboostLSS")
library("mvtnorm")
library("pROC")
library("scoringRules")



## Source marginal distributions
source("FittingInPractice/AFT_RightCensoring_Margins.R")

## Source copula functions
source("FittingInPractice/Bivariate_AFT_RightCensoring_Copulas.R")

## Source wrapper function SurvCopBoost: 
source("FittingInPractice/SurvCopBoost.R")


#### Use SurvCopBoost...
# #### Some example:
# Fit <- SurvCopBoost(form_list, marings = c("WEI", "LOGLOGISTIC"), copula = c("GUMBEL"), 
#                     response_1 = resp1, response_2 = resp2, data = dat, 
#                     mstops = c(1000, 1000, 1000),
#                     oobag_weights = boost_weights, s_step = 0.1, stabilization = "L2")



# where form_list is a list of length 3 containing formulas. (or K, where K is the number of parameters in the joint survival function)
# 3 formulas, assuming that all parameters of the first and second margins have the same formula, as well as the copula parameter.

# resp1, resp2 must be like this: resp1 <- cbind(time1, cens1), resp2 <- cbind(time2, cens2), where
# time1 / time2 is the name of the follow-up time of margin 1 / 2 and cens1 / cens2 is a binary variable that corresponds to the 
# right-censoring indiator. 

# dat is a data frame that has ONLY the covariates. It must NOT contain the variables time1, time2, cens1, cens2. 

# boost_weights is a vector of length nrow(dat) that contains binary entries. The out-of-bag risk is computed using the observations
# with weight equal to zero. 

