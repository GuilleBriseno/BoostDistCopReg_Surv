library("VineCopula")
library("VGAM")
library("gamboostLSS")
library("mvtnorm")
library("pROC")
library("scoringRules")


### load copula functions and marginal distributions
source("FittingInPractice/AFT_RightCensoring_Margins.R")
source("FittingInPractice/Bivariate_AFT_RightCensoring_Copulas.R")


#### Load the data and weights:
load("Application/SCR_APPLICATION_DATA_FILES.RData")



### This function computes the log-score for a univariate right-censored response: 
LogScore_RightCens <- function(PDF, SURV, time, cens){
  
  loglik <- cens * log(PDF) + (1 - cens) * log(SURV)
  
  loglik <- sum(loglik)
  
  return(-loglik)
  
}

#### This function computes the log-score for a bivariate right-censored response: 
LogScore_Biv <- function(PDF1, PDF2, SURV1, SURV2, COPPAR, cens1, cens2, COPFAM){
  
  ### Copula CDF
  COPCDF <- VineCopula::BiCopCDF(u1 = SURV1, u2 = SURV2, family = COPFAM, par = COPPAR)
  
  ### Copula density
  COPDENSITY <- VineCopula::BiCopPDF(u1 = SURV1, u2 = SURV2, family = COPFAM, par = COPPAR)
  
  ### Copula h-function w.r.t. M1
  HWRT_M1 <- VineCopula::BiCopHfunc1(u1 = SURV1, u2 = SURV2, family = COPFAM, par = COPPAR)
  
  ### Copula h-function w.r.t. M2
  HWRT_M2 <- VineCopula::BiCopHfunc1(u2 = SURV1, u1 = SURV2, family = COPFAM, par = COPPAR)
  
  
  loglik <- cens1 * cens2 * ( log(COPDENSITY) + log(PDF1) + log(PDF2) ) + 
    
    (1 - cens1) * cens2 * ( log(HWRT_M2) + log(PDF2) ) +
    
    cens1 * (1 - cens2) * ( log(HWRT_M1) + log(PDF1) ) +
    
    (1 - cens1) * (1 - cens2) * ( log(COPCDF) )
  
  
  loglik <- sum(loglik)
  
  return(-loglik)
}


dim(SCR_DATA_FITTING)
dim(SCR_DATA_EXTERNAL)


## n total
nrow(SCR_DATA_FITTING) + nrow(SCR_DATA_EXTERNAL)

## n_train
nrow(SCR_DATA_FITTING) - sum(weights_mstop)

## n mstop
sum(weights_mstop)

## n external
nrow(SCR_DATA_EXTERNAL)



table(c(SCR_DATA_EXTERNAL$nonterminal_status, SCR_DATA_FITTING$nonterminal_status), 
      c(SCR_DATA_EXTERNAL$terminal_status, SCR_DATA_FITTING$terminal_status))



dim(SCR_DATA_FITTING)
length(weights_mstop)

rownames(SCR_DATA_FITTING) <- NULL
rownames(SCR_DATA_EXTERNAL) <- NULL


table(c(SCR_DATA_EXTERNAL$nonterminal_status, SCR_DATA_FITTING$nonterminal_status), 
      c(SCR_DATA_EXTERNAL$terminal_status, SCR_DATA_FITTING$terminal_status))


############################################################################# DEFINE HYPERPARAMETERS: 
boost_nu <- 0.005


the_stabilization <- "L2"


###### These are for initialisation: 
TryMstop_NONTERMINAL <- 4000

TryMstop_NONTERMINAL_AGE <- 5000


### Weights have been created already: 
table(weights_mstop) / nrow(SCR_DATA_FITTING)


#########################################################################################################
#########################################################################################################
#########################################################################################################     
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
###############################+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ FIT MODEL WITH SOME CLINICAL COVARIATES AND SNPs
#########################################################################################################
#########################################################################################################
#########################################################################################################     
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
observations_FITTING  <- 1:nrow(SCR_DATA_FITTING)
observations_EXTERNAL <- (nrow(SCR_DATA_FITTING) + 1):(nrow(SCR_DATA_FITTING) + nrow(SCR_DATA_EXTERNAL))

ncol(SCR_DATA_FITTING)
head(colnames(SCR_DATA_FITTING)) 
tail(colnames(SCR_DATA_FITTING)) 

non_SNPs_columns <- c(1, 2, 3, 4, 5, 11762:11766)

JUST_SNPS <- rbind(SCR_DATA_FITTING[,-non_SNPs_columns], 
                   SCR_DATA_EXTERNAL[,-non_SNPs_columns])


all.equal(colnames(SCR_DATA_FITTING[,-non_SNPs_columns]), colnames(SCR_DATA_EXTERNAL[,-non_SNPs_columns]))

### re-scale 
JUST_SNPS <- scale(as.matrix( JUST_SNPS ), center = TRUE, scale = TRUE)

dim(JUST_SNPS)



SCR_Nonterminal_Data           <- data.frame(nonterminal_event = SCR_DATA_FITTING$nonterminal_event, 
                                             nonterminal_status = SCR_DATA_FITTING$nonterminal_status,
                                             residual_tumor_size = SCR_DATA_FITTING$residual_tumor_size, 
                                             TumorStage_2 = SCR_DATA_FITTING$TumorStage_2,
                                             TumorStage_3 = SCR_DATA_FITTING$TumorStage_3,
                                             TumorStage_4 = SCR_DATA_FITTING$TumorStage_4,
                                             #
                                             JUST_SNPS[observations_FITTING, ] ) 


SCR_Nonterminal_Data_External  <- data.frame(nonterminal_event = SCR_DATA_EXTERNAL$nonterminal_event, 
                                             nonterminal_status = SCR_DATA_EXTERNAL$nonterminal_status,
                                             residual_tumor_size = SCR_DATA_EXTERNAL$residual_tumor_size, 
                                             TumorStage_2 = SCR_DATA_EXTERNAL$TumorStage_2,
                                             TumorStage_3 = SCR_DATA_EXTERNAL$TumorStage_3,
                                             TumorStage_4 = SCR_DATA_EXTERNAL$TumorStage_4,
                                             #
                                             JUST_SNPS[observations_EXTERNAL, ] ) 


head(colnames(SCR_Nonterminal_Data), 10)
tail(colnames(SCR_Nonterminal_Data), 10)


SCR_Nonterminal_Data$residual_tumor_size <- factor( SCR_Nonterminal_Data$residual_tumor_size )

SCR_Nonterminal_Data$TumorStage_2        <- factor( SCR_Nonterminal_Data$TumorStage_2 )
SCR_Nonterminal_Data$TumorStage_3        <- factor( SCR_Nonterminal_Data$TumorStage_3 )
SCR_Nonterminal_Data$TumorStage_4        <- factor( SCR_Nonterminal_Data$TumorStage_4 )


SCR_Nonterminal_Data_External$residual_tumor_size <- factor( SCR_Nonterminal_Data_External$residual_tumor_size )

SCR_Nonterminal_Data_External$TumorStage_2 <- factor( SCR_Nonterminal_Data_External$TumorStage_2 )
SCR_Nonterminal_Data_External$TumorStage_3 <- factor( SCR_Nonterminal_Data_External$TumorStage_3 )
SCR_Nonterminal_Data_External$TumorStage_4 <- factor( SCR_Nonterminal_Data_External$TumorStage_4 )

dim(SCR_Nonterminal_Data)
dim(SCR_Nonterminal_Data_External)

######################################################################################################### MODEL FORMULA
nonterminal_formula <- formula(cbind(nonterminal_event, nonterminal_status) ~ . )


#########################################################################################################
#########################################################################################################
#########################################################################################################     WEIBULL
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

NON_TERMINAL_EVENT_WEIBULL_MODEL <- glmboostLSS(formula = nonterminal_formula, 
                                                data = SCR_Nonterminal_Data,
                                                families = Custom_WeibullFamily(stabilization = the_stabilization),
                                                control = boost_control(mstop = TryMstop_NONTERMINAL, 
                                                                        nu = boost_nu, 
                                                                        risk = "oobag", 
                                                                        trace = TRUE), 
                                                weights = weights_mstop, 
                                                method = "noncyclic", 
                                                center = FALSE
)

plot(risk(NON_TERMINAL_EVENT_WEIBULL_MODEL, merge = T), type = "l")
which.min(risk(NON_TERMINAL_EVENT_WEIBULL_MODEL, merge = T))

WEIBULL_MSTOP_OPT_OOBAG <- which.min(risk(NON_TERMINAL_EVENT_WEIBULL_MODEL, merge = T))

NON_TERMINAL_EVENT_WEIBULL_MODEL <- glmboostLSS(formula = nonterminal_formula,
                                                data = SCR_Nonterminal_Data,
                                                families = Custom_WeibullFamily(stabilization = the_stabilization),
                                                control = boost_control(mstop = WEIBULL_MSTOP_OPT_OOBAG,
                                                                        nu = boost_nu,
                                                                        trace = TRUE),
                                                method = "noncyclic",
                                                center = FALSE
)


########################################################################################################## COMPUTE LOG SCORES OF TERMINAL EVENT: 
## WEIBULL DISTRIBUTION

### Compute parameters: 
NON_TERMINAL_MU <- exp( predict( NON_TERMINAL_EVENT_WEIBULL_MODEL$mu, type = "link", newdata = SCR_Nonterminal_Data_External) 
)

NON_TERMINAL_SIGMA <- exp( predict( NON_TERMINAL_EVENT_WEIBULL_MODEL$sigma, type = "link", newdata = SCR_Nonterminal_Data_External) 
)

### Compute PDF and SURV: WEIBULL DISTRIBUTION
PDF_NON_TERMINAL  <- dweibull(x = SCR_Nonterminal_Data_External$nonterminal_event, scale = NON_TERMINAL_MU, shape = NON_TERMINAL_SIGMA)
SURV_NON_TERMINAL <- pweibull(q = SCR_Nonterminal_Data_External$nonterminal_event, scale = NON_TERMINAL_MU, shape = NON_TERMINAL_SIGMA, lower.tail = FALSE)

### COMPUTE LOG SCORE on external data
LOGSCORE_WEIBULL_MODEL <- LogScore_RightCens(PDF = PDF_NON_TERMINAL, 
                                             SURV = SURV_NON_TERMINAL, 
                                             time = SCR_DATA_EXTERNAL$nonterminal_event, 
                                             cens = SCR_DATA_EXTERNAL$nonterminal_status)








# #########################################################################################################
# #########################################################################################################
# #########################################################################################################
# #########################################################################################################   LOG LOGISTIC
# #########################################################################################################
# #########################################################################################################
# #########################################################################################################
# #########################################################################################################
# #########################################################################################################

NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL <- glmboostLSS(formula = nonterminal_formula,
                                                    data = SCR_Nonterminal_Data,
                                                    families = Custom_LogLogisticFamily(stabilization = the_stabilization),
                                                    control = boost_control(mstop = TryMstop_NONTERMINAL,
                                                                            nu = boost_nu,
                                                                            risk = "oobag",
                                                                            trace = TRUE),
                                                    weights = weights_mstop,
                                                    method = "noncyclic",
                                                    center = FALSE
)


plot(risk(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL, merge = T), type = "l")
which.min(risk(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL, merge = T))

NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL <- NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL[5000]


LOGLOGISTIC_MSTOP_OPT_OOBAG <- which.min(risk(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL, merge = T))


NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL <- glmboostLSS(formula = nonterminal_formula,
                                                    data = SCR_Nonterminal_Data,
                                                    families = Custom_LogLogisticFamily(stabilization = the_stabilization),
                                                    control = boost_control(mstop = LOGLOGISTIC_MSTOP_OPT_OOBAG,
                                                                            nu = boost_nu,
                                                                            trace = TRUE),
                                                    method = "noncyclic",
                                                    center = FALSE
)



length(coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$mu))
length(coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$sigma))

mstop(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL)



########################################################################################################## COMPUTE LOG SCORES OF TERMINAL EVENT:

### LOG-LOGISTIC DISTRIBUTION  dfisk(x = y[,2], scale = mu2, shape1.a = sigma2)
### Compute parameters:
NON_TERMINAL_MU <-  exp( 
  predict( NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$mu, type = "link", newdata = SCR_Nonterminal_Data_External) 
)

NON_TERMINAL_SIGMA <- exp( 
  predict( NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$sigma, type = "link", newdata = SCR_Nonterminal_Data_External) 
)

### Compute PDF and SURV:
PDF_NON_TERMINAL  <- VGAM::dfisk(SCR_Nonterminal_Data_External$nonterminal_event, scale = NON_TERMINAL_MU, shape1.a = NON_TERMINAL_SIGMA)
SURV_NON_TERMINAL <- 1 - VGAM::pfisk(SCR_Nonterminal_Data_External$nonterminal_event, scale = NON_TERMINAL_MU, shape1.a = NON_TERMINAL_SIGMA)


### COMPUTE LOG SCORE on external data
LOGSCORE_LOGLOGISTIC_MODEL <- LogScore_RightCens(PDF = PDF_NON_TERMINAL,
                                                 SURV = SURV_NON_TERMINAL,
                                                 time = SCR_DATA_EXTERNAL$nonterminal_event,
                                                 cens = SCR_DATA_EXTERNAL$nonterminal_status)


# save(LOGSCORE_LOGLOGISTIC_MODEL, 
#      file = "BIVARIATE_AFT/SCR_Application_Results2.RData")

coeff_list_loglogistic <- list(COVNAMES_MU = names(coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$mu)),
                               COVNAMES_SIGMA = names(coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$sigma)),
                               coef_mu = coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$mu),
                               coef_sigma = coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$sigma),
                               offset_mu = attr(coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$mu), "offset"),
                               offset_sigma = attr(coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$sigma), "offset"),
                               MSTOP = mstop(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL)
                               )





#########################################################################################################
#########################################################################################################
#########################################################################################################
######################################################################################################### LOG NORMAL DISTRIBUTION
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

NON_TERMINAL_EVENT_LOGNORMAL_MODEL <- glmboostLSS(formula = nonterminal_formula, 
                                                  data = SCR_Nonterminal_Data,
                                                  families = Custom_LogNormalFamily(stabilization = the_stabilization),
                                                  control = boost_control(mstop = TryMstop_NONTERMINAL, 
                                                                          nu = boost_nu, 
                                                                          risk = "oobag", 
                                                                          trace = TRUE), 
                                                  weights = weights_mstop,
                                                  method = "noncyclic",
                                                  center = FALSE
)


plot(risk(NON_TERMINAL_EVENT_LOGNORMAL_MODEL, merge = T), type = "l")
which.min(risk(NON_TERMINAL_EVENT_LOGNORMAL_MODEL, merge = T))

NON_TERMINAL_EVENT_LOGNORMAL_MODEL <- NON_TERMINAL_EVENT_LOGNORMAL_MODEL[6000]

LOGNORMAL_MSTOP_OPT_OOBAG <- which.min(risk(NON_TERMINAL_EVENT_LOGNORMAL_MODEL, merge = T))

NON_TERMINAL_EVENT_LOGNORMAL_MODEL <- glmboostLSS(formula = nonterminal_formula, 
                                                  data = SCR_Nonterminal_Data,
                                                  families = Custom_LogNormalFamily(stabilization = the_stabilization),
                                                  control = boost_control(mstop = LOGNORMAL_MSTOP_OPT_OOBAG, 
                                                                          nu = boost_nu, 
                                                                          trace = TRUE), 
                                                  method = "noncyclic",
                                                  center = FALSE
)




########################################################################################################## COMPUTE LOG SCORES OF TERMINAL EVENT: 

### LOGNORMAL DISTRIBUTION
### Compute parameters: 
NON_TERMINAL_MU <- predict( NON_TERMINAL_EVENT_LOGNORMAL_MODEL$mu, type = "link", newdata = SCR_Nonterminal_Data_External)


NON_TERMINAL_SIGMA <- exp( predict( NON_TERMINAL_EVENT_LOGNORMAL_MODEL$sigma, type = "link", newdata = SCR_Nonterminal_Data_External) )

### Compute PDF and SURV: 
PDF_NON_TERMINAL  <- dlnorm(SCR_Nonterminal_Data_External$nonterminal_event, meanlog = NON_TERMINAL_MU, sdlog = NON_TERMINAL_SIGMA)
SURV_NON_TERMINAL <- 1 - plnorm(SCR_Nonterminal_Data_External$nonterminal_event, meanlog = NON_TERMINAL_MU, sdlog = NON_TERMINAL_SIGMA)


### COMPUTE LOG SCORE on external data
LOGSCORE_LOGNORMAL_MODEL <- LogScore_RightCens(PDF = PDF_NON_TERMINAL, 
                                               SURV = SURV_NON_TERMINAL, 
                                               time = SCR_DATA_EXTERNAL$nonterminal_event, 
                                               cens = SCR_DATA_EXTERNAL$nonterminal_status)





tail(risk(NON_TERMINAL_EVENT_WEIBULL_MODEL, merge = T), 1)
tail(risk(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL, merge = T), 1)
tail(risk(NON_TERMINAL_EVENT_LOGNORMAL_MODEL, merge = T), 1)

##### gather all the mstops:
MSTOPS_MODELS     <- list(WEIBULL = WEIBULL_MSTOP_OPT_OOBAG, 
                          LOGLOGISTIC = LOGLOGISTIC_MSTOP_OPT_OOBAG,
                          LOGNORMAL = LOGNORMAL_MSTOP_OPT_OOBAG)




#### PRINT THEM: (lower is better) # -------------------------------> Log LOGISTIC is the best
NON_TERMINAL_LOGSCORES <- list(WEIBULL = LOGSCORE_WEIBULL_MODEL, 
                               LOGLOGISTIC = LOGSCORE_LOGLOGISTIC_MODEL,
                               LOGNORMAL = LOGSCORE_LOGNORMAL_MODEL)

unlist(NON_TERMINAL_LOGSCORES)[order(unlist(NON_TERMINAL_LOGSCORES))]


#### Check out the model:
coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$mu)
coef(NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$sigma)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
#############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
#############################################################################################################################
##############################################################################################################################
##############################################################################################################################
#                                 COMPUTE PDF / SURVIVAL FUNCTIONS OF OPTIMAL DISTRIBUTIONS
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
#############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
#############################################################################################################################
##############################################################################################################################
##############################################################################################################################

####################### The best-fitting distribution  is: 
################ Terminal event ~ LOG-LOGISTIC

## compute eta hat from the best-fitting model with SNPs and other clinical covariates (AGE as P-Spline is built-in in the offset)
ETA_NON_TERMINAL_MU     <- predict( NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$mu, type = "link" )
ETA_NON_TERMINAL_SIGMA  <- predict( NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$sigma, type = "link" )


MU_NON_TERMINAL     <- exp( ETA_NON_TERMINAL_MU )

SIGMA_NON_TERMINAL  <- exp( ETA_NON_TERMINAL_SIGMA ) 



PDF_NON_TERMINAL     <- VGAM::dfisk(SCR_Nonterminal_Data$nonterminal_event, scale = MU_NON_TERMINAL, shape1.a = SIGMA_NON_TERMINAL)

SURV_NON_TERMINAL     <- 1 - VGAM::pfisk(SCR_Nonterminal_Data$nonterminal_event, scale = MU_NON_TERMINAL, shape1.a = SIGMA_NON_TERMINAL)



### On external data:
ETA_MU_EXTERNAL     <- predict( NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$mu, type = "link", newdata = SCR_Nonterminal_Data_External )

ETA_SIGMA_EXTERNAL  <- predict( NON_TERMINAL_EVENT_LOGLOGISTIC_MODEL$sigma, type = "link", newdata = SCR_Nonterminal_Data_External )

MU_EXTERNAL     <- exp( ETA_MU_EXTERNAL )

SIGMA_EXTERNAL  <- exp( ETA_SIGMA_EXTERNAL )


PDF_NON_TERMINAL_EXTERNAL       <- VGAM::dfisk(SCR_Nonterminal_Data_External$nonterminal_event, scale = MU_EXTERNAL, shape1.a = SIGMA_EXTERNAL)

SURV_NON_TERMINAL_EXTERNAL      <- 1 - VGAM::pfisk(SCR_Nonterminal_Data_External$nonterminal_event, scale = MU_EXTERNAL, shape1.a = SIGMA_EXTERNAL)


par(mfrow = c(2, 2))
plot(PDF_NON_TERMINAL)
plot(SURV_NON_TERMINAL)
plot(PDF_NON_TERMINAL_EXTERNAL)
plot(SURV_NON_TERMINAL_EXTERNAL)
par(mfrow = c(1, 1))





##############################################################################################################################
##############################################################################################################################
############### Compute quantities on the EXTERNAL validation data for computation of log-score of copula models!
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
################################################################################### SAVE SOME STUFF
##############################################################################################################################
save(SURV_NON_TERMINAL,
     PDF_NON_TERMINAL,
     PDF_NON_TERMINAL_EXTERNAL,
     SURV_NON_TERMINAL_EXTERNAL,
     LOGSCORE_LOGLOGISTIC_MODEL,
     coeff_list_loglogistic,
     boost_nu,
     the_stabilization, 
     file = "Applications/SCR_APPLICATION_NEWVERSION/NONTERMINAL_BEST_FITTING.RData"
)












