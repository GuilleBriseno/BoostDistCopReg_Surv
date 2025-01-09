library("VineCopula")
library("VGAM")
library("mvtnorm")
library("pROC")
library("scoringRules")
library("mboost")


### load copula functions and marginal distributions
source("FittingInPractice/AFT_RightCensoring_Margins.R")
source("FittingInPractice/Bivariate_AFT_RightCensoring_Copulas.R")



#### Load the data and weights:
load("Application/SCR_APPLICATION_DATA_FILES.RData")

#### Load quantities from the TERMINAL AND NON-TERMINAL MARGIN EVENTS
load("Application/SCR_APPLICATION_NEWVERSION/NONTERMINAL_BEST_FITTING.RData")
load("Application/SCR_APPLICATION_NEWVERSION/TERMINAL_BEST_FITTING.RData")

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


############################################################################# DEFINE HYPERPARAMETERS: 
boost_nu <- 0.005

the_stabilization <- "L2"

###### These are for initialisation: 
TryMSTOP_Dependence <- 3000

### Weights have been created already: 
table(weights_mstop) / nrow(SCR_DATA_FITTING)


######################################################################################################################################
######################################################################################################################################
######################################################## FINAL ADJUSTMENT TO  DATA
######################################################################################################################################
######################################################################################################################################
par(mfrow = c(2, 2))
plot(SURV_TERMINAL)
plot(SURV_NON_TERMINAL)
plot(PDF_TERMINAL)
plot(PDF_NON_TERMINAL)
par(mfrow = c(1,1))


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

sum(colMeans(JUST_SNPS))
sum(apply(JUST_SNPS, 2, sd)) == ncol(JUST_SNPS)


#################################################################### Design matrix:
residual_tumor_size_numeric           <- ifelse(as.numeric(SCR_DATA_FITTING$residual_tumor_size) == 2, 1, 0)
residual_tumor_size_numeric_external  <- ifelse(as.numeric(SCR_DATA_EXTERNAL$residual_tumor_size) == 2, 1, 0)

MatrixOfClinicalCovariates <- as.matrix(data.frame(ResidualTumorSize = residual_tumor_size_numeric, 
                                                   TumorStage2 = as.numeric(SCR_DATA_FITTING$TumorStage_2),
                                                   TumorStage3 = as.numeric(SCR_DATA_FITTING$TumorStage_3),
                                                   TumorStage4 = as.numeric(SCR_DATA_FITTING$TumorStage_4),
                                                   XInter = rep(1, nrow(SCR_DATA_FITTING))) 
)

MatrixOfClinicalCovariates_EXTERNAL <- as.matrix(data.frame(ResidualTumorSize = residual_tumor_size_numeric_external, 
                                                            TumorStage2 = as.numeric(SCR_DATA_EXTERNAL$TumorStage_2),
                                                            TumorStage3 = as.numeric(SCR_DATA_EXTERNAL$TumorStage_3),
                                                            TumorStage4 = as.numeric(SCR_DATA_EXTERNAL$TumorStage_4),
                                                            XInter = rep(1, nrow(SCR_DATA_EXTERNAL)))
)



DesignMatrix_DependenceModel <- cbind(MatrixOfClinicalCovariates, 
                                      JUST_SNPS[observations_FITTING,])

DesignMatrix_DependenceModel_External <- cbind(MatrixOfClinicalCovariates_EXTERNAL, 
                                               JUST_SNPS[observations_EXTERNAL,])

rownames(DesignMatrix_DependenceModel) <- NULL
rownames(DesignMatrix_DependenceModel_External) <- NULL

class(DesignMatrix_DependenceModel)
class(DesignMatrix_DependenceModel_External)

head(colnames(DesignMatrix_DependenceModel), 10)
tail(colnames(DesignMatrix_DependenceModel), 10)

#### (Survival function, pdf, status) of each margin:
ResponseVector_DependenceModel <- cbind(SURV_NON_TERMINAL, 
                                        PDF_NON_TERMINAL, 
                                        SCR_DATA_FITTING$nonterminal_status,
                                        SURV_TERMINAL, 
                                        PDF_TERMINAL, 
                                        SCR_DATA_FITTING$terminal_status)


################################################################################################################## 

### ALL COPULA FUNCTIONS: 

# BivAFT_RC_GaussCopula_RhoSoloFamily()
# BivAFT_RC_FrankCopula_RhoSoloFamily()
#
# BivAFT_RC_ClaytonCopula_RhoSoloFamily()
# BivAFT_RC_ClaytonCopula_90_RhoSoloFamily()
# BivAFT_RC_ClaytonCopula_180_RhoSoloFamily()
# BivAFT_RC_ClaytonCopula_270_RhoSoloFamily()
#
# BivAFT_RC_GumbelCopula_RhoSoloFamily()
# BivAFT_RC_GumbelCopula_90_RhoSoloFamily()
# BivAFT_RC_GumbelCopula_180_RhoSoloFamily()
# BivAFT_RC_GumbelCopula_270_RhoSoloFamily()
#
# BivAFT_RC_JoeCopula_RhoSoloFamily()
# BivAFT_RC_JoeCopula_90_RhoSoloFamily()
# BivAFT_RC_JoeCopula_180_RhoSoloFamily()
# BivAFT_RC_JoeCopula_270_RhoSoloFamily()


# rotated copulas by 90 and 270 degrees are not supported by the data
names_of_copulas <- c("GAUSS",
                      "C0", "C180", 
                      "G0", "G180", 
                      "J0", "J180", 
                      #
                      #
                      "C90", "G90", "J90", # 90 degrees
                      "C270", "G270", "J270", # 270 degrees
                      "FRANK")

List_Of_Copula_Functions <- list(BivAFT_RC_GaussCopula_RhoSoloFamily,
                                 #
                                 BivAFT_RC_ClaytonCopula_RhoSoloFamily,
                                 BivAFT_RC_ClaytonCopula_180_RhoSoloFamily,
                                 #
                                 BivAFT_RC_GumbelCopula_RhoSoloFamily,
                                 BivAFT_RC_GumbelCopula_180_RhoSoloFamily,
                                 #
                                 BivAFT_RC_JoeCopula_RhoSoloFamily,
                                 BivAFT_RC_JoeCopula_180_RhoSoloFamily,
                                 #
                                 ##### NOT SUPPORTED BY DATA
                                 #
                                 BivAFT_RC_ClaytonCopula_90_RhoSoloFamily, # 90 degrees
                                 BivAFT_RC_GumbelCopula_90_RhoSoloFamily,
                                 BivAFT_RC_JoeCopula_90_RhoSoloFamily,
                                 #
                                 #
                                 BivAFT_RC_ClaytonCopula_270_RhoSoloFamily, # 270 degrees
                                 BivAFT_RC_GumbelCopula_270_RhoSoloFamily,
                                 BivAFT_RC_JoeCopula_270_RhoSoloFamily,
                                 #
                                 # 
                                 BivAFT_RC_FrankCopula_RhoSoloFamily
)


List_Of_MSTOPS              <- vector(mode = "list", length = length(List_Of_Copula_Functions))

List_Of_RISKS               <- vector(mode = "list", length = length(List_Of_Copula_Functions))

List_Of_ETA_COPPAR_EXTERNAL <- vector(mode = "list", length = length(List_Of_Copula_Functions))

LIST_OF_COPULA_LOGSCORES    <- vector(mode = "list", length = length(List_Of_Copula_Functions))


names(List_Of_ETA_COPPAR_EXTERNAL)    <- names_of_copulas
names(List_Of_MSTOPS)                 <- names_of_copulas
names(List_Of_RISKS)                  <- names_of_copulas
names(LIST_OF_COPULA_LOGSCORES)       <- names_of_copulas


##### This loop fits a model to the implemented copula functions, use it up to 7 if you only want the copulas actually 
# supported by the data. 
# alternatively, let the loop run until length(List_Of_Copula_Functions) to try out ALL copulas. 
for(i in 1:length(List_Of_Copula_Functions)){
  
  
  ### Once we fit the copulas supported by the data, we increase the MSTOP just to check if the remaining
  ### copulas reach an optimal mstop: 
  if( i > 7 & i < 14){
    
    print("candidate MSTOP has been increased to 15000 iterations.")
    
    TryMSTOP_Dependence <- 15000
    
  }
  
  ## Only for FRANK copula
  if( i == 14){
    
    print("candidate MSTOP has been increased to 50000 iterations.")
    
    TryMSTOP_Dependence <- 50000
    
  }

  print(paste("Fitting the following copula:", names_of_copulas[i]))

  Copula_Family <- List_Of_Copula_Functions[[i]]


  COPULA_MODEL <- mboost::glmboost(y = ResponseVector_DependenceModel, 
                                   x = DesignMatrix_DependenceModel, 
                                   family = Copula_Family(),
                                   weights = weights_mstop, 
                                   control = mboost::boost_control(mstop = TryMSTOP_Dependence,
                                                                   nu = boost_nu,
                                                                   risk = "oobag", 
                                                                   trace = FALSE), 
                                   center = FALSE)

  COPULA_MODEL_MSTOP <- which.min(risk(COPULA_MODEL))

  List_Of_MSTOPS[[i]] <- COPULA_MODEL_MSTOP
  
  List_Of_RISKS[[i]]  <- risk(COPULA_MODEL)

  COPULA_MODEL <- mboost::glmboost(y = ResponseVector_DependenceModel, 
                                   x = DesignMatrix_DependenceModel, 
                                   family = Copula_Family(), 
                                   control = mboost::boost_control(mstop = COPULA_MODEL_MSTOP,
                                                                   nu = boost_nu,
                                                                   trace = FALSE), 
                                   center = FALSE)


  List_Of_ETA_COPPAR_EXTERNAL[[i]] <- as.numeric( predict( COPULA_MODEL, type = "link", newdata = DesignMatrix_DependenceModel_External ) )

  rm(COPULA_MODEL)
  
  
  
  # compute the log scores:
  if(i == 1){
    
    print("Computing log score of Gaussian copula")
    
    COPPAR_HAT_EXTERNAL <- tanh( List_Of_ETA_COPPAR_EXTERNAL[[i]] )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 1,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 2){
    
    print("Computing log score of Clayton 0 copula")
    
    COPPAR_HAT_EXTERNAL <- exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 3,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 3){
    
    print("Computing log score of Clayton 180 copula")
    
    COPPAR_HAT_EXTERNAL <- exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 13,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 4){
    
    print("Computing log score of Gumbel 0 copula")
    
    COPPAR_HAT_EXTERNAL <- exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] ) + 1
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 4,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 5){
    
    print("Computing log score of Gumbel 180 copula")
    
    COPPAR_HAT_EXTERNAL <- exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] ) + 1
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 14,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 6){
    
    print("Computing log score of Joe 0 copula")
    
    COPPAR_HAT_EXTERNAL <- exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] ) + 1
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 6,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 7){
    
    print("Computing log score of Joe 180 copula")
    
    COPPAR_HAT_EXTERNAL <- exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] ) + 1
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 16,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  ### Rotations (these are not supported by the data)
  # 90
  if(i == 8){
    
    print("Computing log score of Clayton 90 copula")
    
    COPPAR_HAT_EXTERNAL <-  - exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 23,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 9){
    
    print("Computing log score of Gumbel 90 copula")
    
    COPPAR_HAT_EXTERNAL <- - ( exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] ) + 1 )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 24,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 10){
    
    print("Computing log score of Joe 90 copula")
    
    COPPAR_HAT_EXTERNAL <- - ( exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] ) + 1 )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 26,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  ## 270
  if(i == 11){
    
    print("Computing log score of Clayton 270 copula")
    
    COPPAR_HAT_EXTERNAL <- - exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 33,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 12){
    
    print("Computing log score of Gumbel 270 copula")
    
    COPPAR_HAT_EXTERNAL <- - ( exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] ) + 1 )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 34,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  if(i == 13){
    
    print("Computing log score of Joe 270 copula")
    
    COPPAR_HAT_EXTERNAL <- - ( exp( List_Of_ETA_COPPAR_EXTERNAL[[i]] ) + 1 )
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 36,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
  ### FRANK Copula
  if(i == 14){
    
    print("Computing log score of Frank copula")
    
    COPPAR_HAT_EXTERNAL <- List_Of_ETA_COPPAR_EXTERNAL[[i]]
    
    ## Compute log-score
    LOGSCORE_COPULA <- LogScore_Biv(COPFAM = 5,
                                    PDF1 = PDF_NON_TERMINAL_EXTERNAL, 
                                    PDF2 = PDF_TERMINAL_EXTERNAL, 
                                    SURV1 = SURV_NON_TERMINAL_EXTERNAL, 
                                    SURV2 = SURV_TERMINAL_EXTERNAL, 
                                    COPPAR = COPPAR_HAT_EXTERNAL, 
                                    cens1 = SCR_DATA_EXTERNAL$nonterminal_status, 
                                    cens2 = SCR_DATA_EXTERNAL$terminal_status
    )
    
    LIST_OF_COPULA_LOGSCORES[[i]] <- LOGSCORE_COPULA
    
  }
  
}



List_Of_MSTOPS

##### LOGSCORE of the independence copula (no dependence between the margins):
LOGSCORE_INDEPENDENCE <- TERMINAL_LOGSCORES[[1]] + NON_TERMINAL_LOGSCORES[[2]]


##### COMPARE INDEPENDENCE VS THE BEST FITTING COPULA (GUMBEL) (lower is better)
unlist(LIST_OF_COPULA_LOGSCORES[order(unlist(LIST_OF_COPULA_LOGSCORES))])

LOGSCORE_INDEPENDENCE


#### RE-FIT Optimal model: 
OPTIMAL_COPULA_MODEL <- mboost::glmboost(y = ResponseVector_DependenceModel, 
                                         x = DesignMatrix_DependenceModel, 
                                         family = BivAFT_RC_GumbelCopula_RhoSoloFamily(),
                                         control = mboost::boost_control(mstop = List_Of_MSTOPS[[4]],
                                                                         nu = boost_nu,
                                                                         trace = FALSE), 
                                         center = FALSE)


summary(OPTIMAL_COPULA_MODEL)

COEFF_OPTIMAL_COPULA_MODEL <- coef(OPTIMAL_COPULA_MODEL)



##### 
save(List_Of_MSTOPS, 
     List_Of_RISKS,
     List_Of_ETA_COPPAR_EXTERNAL,
     LIST_OF_COPULA_LOGSCORES,
     LOGSCORE_INDEPENDENCE,
     COEFF_OPTIMAL_COPULA_MODEL,
     file = "Applications/SCR_APPLICATION/SCR_APPLICATION_DEPENDENCE_MODEL_OUTPUT.RData"
)






