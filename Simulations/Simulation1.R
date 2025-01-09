
library("VineCopula")
library("VGAM")
library("gamboostLSS")
library("mvtnorm")
library("pROC")
library("scoringRules")

setwd("../../BoostingCopulas_Papers/SurvCopBoost_CODE_FILES/")


## Source marginal distributions
source("FittingInPractice/AFT_RightCensoring_Margins.R")

## Source copula functions
source("FittingInPractice/Bivariate_AFT_RightCensoring_Copulas.R")

## Source simulation functions: 
source("Simulations/SemiCompetingRisks_RightCensoring_SimulationFunction.R")
source("Copulas/BivariateAFT/Bivariate_AFT_RightCensoring_SimulationFunction_TOEPLITZ.R")



#### careful with your number of cores!
the_cores <- 1


n.train <- 1000
n.mstop <- 1000
n.test <- 1000

# Set different number of candidate covariates in the model!
#p <- 1000
#p <- 500
p <- 10



CHOSENCOPULA <- 3 # Clayton copula without rotation


# Censoring rate
censoring_rate <- "mild"


# step length
boost.nu.steplength <- 0.1


TryMSTOP_Margins <- 400 # for p = 1000...   wei x loglogistic NONLINEAR
TryMSTOP_Dependence <- 750 ## for p = 1000.     wei x loglogistic NONLINEAR




## number of replications
the_seeds <- 1:500


results <- mclapply(the_seeds,
                    #sim_TwoStage_Estimation_TOEP,                # Use this line for LINEAR DGP
                    sim_TwoStage_Estimation_NONLINEARDGP_TOEP,    # Use this line for NON-LINEAR DGP
                    mc.cores = the_cores,
                    COPFAM = CHOSENCOPULA,
                    n.train = n.train,
                    n.mstop = n.mstop,
                    n.test = n.test,
                    censoring = censoring_rate,
                    TryMSTOP_Margins = TryMSTOP_Margins,
                    TryMSTOP_Dependence = TryMSTOP_Dependence,
                    boost.nu.steplength = boost.nu.steplength,
                    p  = p,
                    mc.preschedule = FALSE)



#save(results, file = "")

