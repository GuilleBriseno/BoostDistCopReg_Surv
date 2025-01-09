source("Simulations/SemiCompetingRisks_Copula_Gumbel_BivAFT.R")

### This function just computs the Brier score.
ComputeBrierScore <- function(S_hat, G_hat_atT, G_hat_atTSTAR, delta, time_star, time_i){
  
  S_hat <- pdffz(S_hat)
  
  G_hat_atT <- pdffz(G_hat_atT)
  
  G_hat_atTSTAR <- pdffz(G_hat_atTSTAR)
  
  Term1 <- ( S_hat )^2 / ( G_hat_atT )
  
  Term2 <- ( 1 - S_hat )^2 / ( G_hat_atTSTAR )
  
  Case1 <- as.numeric((time_i < time_star) & delta == 1)
  
  Case2 <- as.numeric( (time_i >= time_star) )
  
  BS_i <- Term1*Case1 + Term2*Case2
  
  ### NEW CHECK
  BS_i <- ifelse(BS_i > 1, 1, BS_i)
  
  BS <- mean(BS_i)
  
  return(BS)
  
}

### This function runs the simulations, the arguments are hopefully self-explanatory:
# seed                            <- seed for initialising random numbers

# p                               <- number of covariates in the data (10, 500, 1000)

# corr_toe                        <- correlation of the covariance matrix used to generate Toeplitz covariates  

# censoring                       <- censoring rate of the margins (ONLY mild)

# TryMSTOP_Margins                <- Initial, non-optimal number of fitting iterations for the margins
# TryMSTOP_Dependence_SemiComp    <- Initial, non-optimal number of fitting iterations for the SCR estimation loss
# TryMSTOP_Dependence             <- Initial, non-optimal number of fitting iterations for the dependence parameter

# n.train, n.mstop, n.test,       <- number of observations used for fitting, tuning mstop and evaluate the performance metrics

# boost.nu.steplength             <- step length of boosting algorithm.

# COPFAM                          <- Copula family according to VineCopula's terminology (COPFAM = 4 (GUMBEL) for semi-competing risks) 
#
#
sim_TwoStage_Estimation_TOEPLITZ <- function(seed, p, corr_toe = 0.5, 
                                             censoring = "mild", 
                                             TryMSTOP_Margins = 400, 
                                             TryMSTOP_Dependence_SemiComp = 1000, 
                                             TryMSTOP_Dependence = 2500, 
                                             n.train, n.mstop, n.test, 
                                             boost.nu.steplength = 0.1, 
                                             COPFAM = 3){
  
  
  
  data.gen.biv <- function(n, mu1, mu2, sigma1, sigma2, theta, nu1 = NULL, nu2 = NULL, FAM = NULL, timecut){
    
    ######### 
    U_sim <- runif(n, 0, 1)
    V_sim <- VineCopula::BiCopCondSim(N = n, cond.val = U_sim, cond.var = 1, family = FAM, par = theta)
    

    
    # transform into non-uniform random variables:
    y1 <- qweibull(p = U_sim, scale = mu1, shape = sigma1, lower.tail = FALSE)
    # range(y1)
    # plot(y1)
    
    y2 <- qfisk(p = V_sim, scale = mu2, shape1.a = sigma2, lower.tail = FALSE)
    # range(y2)
    # plot(y2)
    
    plot(y1, y2)
    abline(a = 0, b = 1, col = "red")
    
    # Adding random, non-informative censoring:
    censoring_time <- runif(length(mu1), min = 0, max = timecut)
    
    ### All three times: 
    all_three_times <- data.frame(y1 = y1, y2 = y2, cens = censoring_time)
    
    only_two_times <- data.frame(y2 = y2, cens = censoring_time)
    
    ####### Determine censoring status of non-terminal event: 
    minimum_value <- c()
    
    margin1 <- c()
    margin2 <- c()
    
    status_m1 <- c()
    status_m2 <- c()
    
    
    minimum_value <- c()
    
    margin1 <- c()
    margin2 <- c()
    
    status_m1 <- c()
    status_m2 <- c()
    
    minimum_of_all_three <- apply(all_three_times, 1, which.min)
    minimum_of_remaining_two <- apply(only_two_times, 1, which.min)
    
    
    for(i in 1:length(mu1)){
      
      # Margin 1 is an event:
      if(minimum_of_all_three[i] == 1){
        
        margin1[i] <- all_three_times$y1[i]
        
        status_m1[i] <- 1
        
        
        if(minimum_of_remaining_two[i] == 1){
          
          margin2[i] <- all_three_times$y2[i]
          
          status_m2[i] <- 1
          
        }else{
          
          margin2[i] <- all_three_times$cens[i]
          
          status_m2[i] <- 0
          
        }
        
      }
      
      
      # Margin 2 is an event:  
      if(minimum_of_all_three[i] == 2){
        
        
        margin1[i] <- all_three_times$y2[i]
        margin2[i] <- all_three_times$y2[i]
        
        # death censored the non-terminal event
        status_m1[i] <- 0
        status_m2[i] <- 1
        
      }
      
      
      if(minimum_of_all_three[i] == 3){
        
        margin1[i] <- all_three_times$cens[i]  
        margin2[i] <- all_three_times$cens[i]
        
        status_m1[i] <- 0
        status_m2[i] <- 0
      }  
      
    }
    
    
    
    
    # # censoring rate margin 1: 
    table(status_m1)/n
    # 
    # # censoring rate margin 2: 
    table(status_m2)/n
    
    # Replace TRUE event times with observed times: 
    y1 <- margin1
    y2 <- margin2
    
    plot(y1, y2)
    
    dat <- data.frame(y1, y2, status_m1, status_m2)
    
    return(dat)
  }
  
  
  loss <- function(mu1, mu2, sigma1, sigma2, rho, y, delta1, delta2, nu1 = NULL, nu2 = NULL, FAM = NULL){
    
    
    # Margin 1: Weibull
    # Survival function
    S1 <- pdffz( 1 - pweibull(q = y[,1], scale = mu1, shape = sigma1) )
    
    # density
    f1 <- dweibull(x = y[,1], scale = mu1, shape = sigma1)
    f1 <- ifelse(f1 < 1e-40, 1e-40, f1)
    
    # Margin 2: LogLogistic
    S2 <- pdffz( 1 - pfisk(q = y[,2], scale = mu2, shape1.a = sigma2) )
    
    # density
    f2 <- dfisk(x = y[,2], scale = mu2, shape1.a = sigma2)
    f2 <- ifelse(f2 < 1e-40, 1e-40, f2)
    
    ###################################################
    # Copula parameter 
    thet <- rho
    
    ### Compute the copula terms: 
    copdensity <- VineCopula::BiCopPDF(u1 = S1, u2 = S2, family = FAM, par = thet)
    
    hfunc_m1 <- VineCopula::BiCopHfunc1(u1 = S1, u2 = S2, family = FAM, par = thet)
    
    hfunc_m2 <- VineCopula::BiCopHfunc1(u2 = S1, u1 = S2, family = FAM, par = thet) # I switch them here because hfunc1 is allegedly faster... 
    
    copcdf <- VineCopula::BiCopCDF(u1 = S1, u2 = S2, family = FAM, par = thet)
    ################################
    
    # both uncensored
    uncens <- ( log( f1 ) + log( f2 ) + log(  copdensity  ) )
    
    # margin one is censored
    onecens <- ( log( f2 ) +  log( pdffz( hfunc_m2 ) )  )
    
    # margin two is censored
    twocens <- ( log( f1 ) +  log( pdffz( hfunc_m1 ) )  )
    
    # both are censored
    bothcens <- log( pdffz( copcdf ) )
    
    L_Total <-   delta1 * delta2             * ( uncens ) + 
      ( 1 - delta1 ) * delta2     * ( onecens ) + 
      delta1 * ( 1 - delta2 )     * ( twocens ) +
      ( 1 - delta1 ) * ( 1 - delta2 ) * ( bothcens ) 
    
    # Negative log-likelihood
    return(- sum( L_Total ) )
    
  }
  
  loss_M1 <- function(mu1, sigma1, y1, delta1){
    
    # Margin 1: Weibull
    # Survival function
    S1 <- pdffz( 1 - pweibull(q = y1, scale = mu1, shape = sigma1) )
    
    # density
    f1 <- dweibull(x = y1, scale = mu1, shape = sigma1)
    f1 <- ifelse(f1 < 1e-40, 1e-40, f1)
    
    ################################
    
    L_Total <- delta1 * log( f1 ) + ( 1 - delta1 ) * log( S1 )
    
    # Negative log-likelihood
    return(- sum( L_Total ) )
    
  }
  
  loss_M2 <- function(mu2, sigma2, y2, delta2){
    
    
    # Margin 2: LogLogistic
    S2 <- pdffz( 1 - pfisk(q = y2, scale = mu2, shape1.a = sigma2) )
    
    # density
    f2 <- dfisk(x = y2, scale = mu2, shape1.a = sigma2)
    f2 <- ifelse(f2 < 1e-40, 1e-40, f2)
    
    ################################
    
    
    L_Total <-  delta2 * log( f2 ) + ( 1 - delta2 ) * log( S2 )
    
    # Negative log-likelihood
    return(- sum( L_Total ) )
    
  }
  
  
  
  beta11     <- c( 2, 0, 0, 0, 0, 0, rep(0, p-6)) 
  
  beta12    <- c( 0, +1, 0, 1.5, 0, 0, rep(0, p-6)) 
  
  beta21    <- c( 1, 1.5, 0, 0, 0, 0, rep(0, p-6)) 
  
  beta22    <- c( 0, 0.75, 0, +0.75, 0, 0, rep(0, p-6)) 
  
  betarho   <- c( 0, -2, 0, -2, +0, 0, rep(0, p-6)) 
  
  
  
  n <- n.train + n.mstop
  
  weights.mstop <- c(rep(1, times = n.train), rep(0, times = n.mstop)) 
  
  set.seed(seed)
  
  #### sample design matrix for train:
  # - training data
  BigXMatrix_Train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr_toe^x)))
  BigXMatrix_Train <- apply(BigXMatrix_Train, 2, pnorm)
  
  x.train <- BigXMatrix_Train
  
  #### sample design matrix for test / evaluation:
  BigXMatrix_Test <- rmvnorm(n = n.test, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr_toe^x)))
  BigXMatrix_Test <- apply(BigXMatrix_Test, 2, pnorm) 
  
  x.test <- BigXMatrix_Test
  
  colnames(x.train) <- paste0("X", 1:p)
  colnames(x.test) <- paste0("X", 1:p)
  
  if(censoring == "mild"){
    
    ########################### 
    eta_11 <- x.train %*% beta11
    eta_12 <- x.train %*% beta12
    
    eta_21 <- x.train %*% beta21
    eta_21 <- 0 + eta_21 
    
    eta_22 <- x.train %*% beta22
    eta_22 <- +1 + eta_22
    
    ############################## test / evaluation
    eta_11_test <- x.test %*% beta11
    eta_12_test <- x.test %*% beta12
    
    eta_21_test <- x.test %*% beta21
    eta_21_test <- 0 + eta_21_test
    
    eta_22_test <- x.test %*% beta22
    eta_22_test <- +1 + eta_22_test
    
    
    
    ### for the cutting / right-censoring of uniform times:
    timecutM1 <- 7 #10 # the value 7 seems to give results similar to our data of the application 
    
  }
  
  
  ### There is no heavy censoring rate for this simulation study
  # if(censoring == "heavy"){
  #   
  #   ### Only the intercepts change for HEAVY CENSORING: 
  #   ############################# 70% censoring
  #   eta_11 <- x.train %*% beta11
  #   eta_11 <- + 1 + eta_11
  #   
  #   eta_12 <- x.train %*% beta12
  #   
  #   
  #   eta_21 <- x.train %*% beta21
  #   eta_21 <- +0.5 + eta_21
  #   
  #   eta_22 <- x.train %*% beta22
  #   eta_22 <- +2 + eta_22
  #   
  #   ############################## test / evaluation
  #   eta_11_test <- x.test %*% beta11
  #   eta_11_test <- + 1 + eta_11_test
  #   
  #   eta_12_test <- x.test %*% beta12
  #   
  #   
  #   eta_21_test <- x.test %*% beta21
  #   eta_21_test <- +0.5 + eta_21_test
  #   
  #   eta_22_test <- x.test %*% beta22
  #   eta_22_test <- +2 + eta_22_test
  #   
  #   
  #   
  #   ### for the cutting / right-censoring of uniform times:
  #   timecutM1 <- 8.5
  #   timecutM2 <- 8.5
  #   
  # }
  
  
  
  ################################ DEPENDENCE: 
  eta_rho <- x.train %*% betarho
  eta_rho <- 2.5 + eta_rho 
  
  
  eta_rho_test <- x.test %*% betarho
  eta_rho_test <- 2.5 + eta_rho_test
  
  
  
  ################################ distribution parameters: 
  mu1     <- exp(eta_11)
  sigma1  <- exp(eta_12)
  
  mu2     <- exp(eta_21)
  sigma2  <- exp(eta_22)
  
  rho     <- exp(eta_rho) + 1
  
  ###########################
  mu1_test      <- exp(eta_11_test)
  sigma1_test   <- exp(eta_12_test)
  
  mu2_test      <- exp(eta_21_test)
  sigma2_test   <- exp(eta_22_test)
  
  rho_test      <- exp(eta_rho_test) + 1
  
  
  plot(rho)
  plot(VineCopula::BiCopPar2Tau(family = 4, rho))
  
  TRUE_KENDALL_RANGE <- range(VineCopula::BiCopPar2Tau(family = COPFAM, rho))
  
  ################################################################# SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.biv(mu1 = mu1, 
                          mu2 = mu2, 
                          sigma1 = sigma1, 
                          sigma2 = sigma2, 
                          theta = rho, 
                          timecut = timecutM1,
                          n = length(mu1), 
                          FAM = as.numeric(COPFAM))
  
  
  y.test <- data.gen.biv(mu1 = mu1_test, 
                         mu2 = mu2_test, 
                         sigma1 = sigma1_test, 
                         sigma2 = sigma2_test, 
                         theta = rho_test, 
                         timecut = timecutM1,
                         n = length(mu1_test), 
                         FAM = as.numeric(COPFAM))
  
  colnames(y.train) <- c("y1", "y2", "delta1", "delta2")
  dat.train <- data.frame(y.train, x.train)
  
  colnames(y.test) <- c("y1", "y2", "delta1", "delta2")
  dat.test <- data.frame(y.test, x.test)
  
  ### censoring rates:
  CensoringRates <- list()
  CensoringRates$Margin1 <- table(y.train$delta1)/n
  CensoringRates$Margin2 <- table(y.train$delta2)/n
  CensoringRates
  
  CensoringRates$CensTable <- table(y.train$delta1, y.train$delta2)
  
  ################## Assemble the datasets
  data_margin1 <- as.data.frame(cbind(dat.train$y1, dat.train$delta1, x.train))
  
  data_margin2 <- as.data.frame(cbind(dat.train$y2, dat.train$delta2, x.train))
  
  preds <- colnames(x.train)
  
  data_margin1$Xint <- rep(1, nrow(data_margin1))
  
  x.train <- cbind(rep(1, length.out = nrow(x.train)), x.train)
  
  colnames(x.train)[1] <- "Xint"
  
  
  ############################################################################################### FITTING 1: MARGINS
  
  
  ## declare formulas: 
  formula_margin1 <- as.formula(paste("cbind(y1, delta1, delta2, SURV2, PDF2) ~", paste0(preds, collapse = "+")))
  
  formula_margin2 <- as.formula(paste("cbind(y2, delta2) ~", paste0(preds, collapse = "+")))
  
  
  Margin1_Formula <- formula_margin1
  
  Margin2_Formula <- list(mu = formula_margin2, 
                          sigma = formula_margin2)
  
  cat("Fitting terminal event (margin 2)...\n")
  
  ############ FIT SEPARATELY:
  
  # Bivariate copula model
  GLM2TRY <- myTryCatch( glmboostLSS(formula = Margin2_Formula, 
                                     data = dat.train, 
                                     families = Custom_LogLogistiFamily(stabilization = "L2"), 
                                     weights = weights.mstop, 
                                     method = "noncyclic",
                                     control = boost_control(mstop = TryMSTOP_Margins, 
                                                             nu = boost.nu.steplength, 
                                                             risk = "oobag", 
                                                             trace = TRUE)) )
  
  if(is.null(GLM2TRY$value)){
    
    output <- list(ERROR = GLM2TRY$the_error, 
                   WARNING = GLM2TRY$the_warning)
    
    return(output)
    
  }else{
    
    
    GLM_M2 <- GLM2TRY$value
    
    
    risk_duringtraining_M2 <- risk(GLM_M2, merge = TRUE)
    
    MSTOP_OPT_M2 <- which.min(risk(GLM_M2, merge = TRUE))
    
    MSTOP_OPT_M2
    
    rm(GLM_M2)
    
    GLM_M2 <- glmboostLSS(formula = Margin2_Formula, 
                          data = dat.train, 
                          families = Custom_LogLogisticFamily(stabilization = "L2"), 
                          method = "noncyclic",
                          control = boost_control(mstop = MSTOP_OPT_M2, 
                                                  nu = boost.nu.steplength, 
                                                  trace = TRUE))
    
    coef(GLM_M2)
    
    
    risk_final_M2 <- risk(GLM_M2, merge = TRUE)
    
    ############################### Compute survival and pdf
    mu2_hat <- exp( predict(GLM_M2$mu, type = "link") )
    sigma2_hat <- exp( predict(GLM_M2$sigma, type = "link") )
    
    SURV2 <- 1 - pLogLogistic(x = dat.train$y2, mu = mu2_hat, sigma = sigma2_hat)
    PDF2 <- dLogLogistic(x = dat.train$y2, mu = mu2_hat, sigma = sigma2_hat)
    
    
    dat.train$SURV2 <- SURV2
    dat.train$PDF2 <- PDF2
    
    
    ############################################################################################# Margin 2:
    
    cat("Fitting non-terminal event and copula (margin 1)...\n")
    
    GLM_M1 <- glmboostLSS(formula = Margin1_Formula, 
                          data = dat.train,
                          families = Gumbel_Copula_SemiCompeting(marg1 = "WEI", stabilization = "L2"
                          ), 
                          weights = weights.mstop, 
                          method = "noncyclic",
                          control = boost_control(mstop = TryMSTOP_Dependence_SemiComp, 
                                                  nu = boost.nu.steplength,
                                                  risk = "oobag", 
                                                  trace = TRUE))
    
    risk_duringtraining_M1 <- risk(GLM_M1, merge = TRUE)
    
    MSTOP_OPT_M1 <- which.min(risk(GLM_M1, merge = TRUE))
    
    MSTOP_OPT_M1
    
    rm(GLM_M1)
    
    GLM_M1 <- glmboostLSS(formula = Margin1_Formula, 
                          data = dat.train, 
                          families = Gumbel_Copula_SemiCompeting(marg1 = "WEI",
                                                                 stabilization = "L2"), 
                          method = "noncyclic",
                          control = boost_control(mstop = MSTOP_OPT_M1,
                                                  nu = boost.nu.steplength,
                                                  trace = TRUE))
    
    
    coef(GLM_M1)
    
    risk_final_M1 <- risk(GLM_M1, merge = TRUE)
    
    
    ################################################################################################################
    ################################################################################################################
    ##### Predict some quantities required for the copula parameter: 
    mu1_hat <- exp( predict(GLM_M1$mu1, type = "link") )
    sigma1_hat <- exp( predict(GLM_M1$sigma1, type = "link") )
    
    
    
    #######################################
    SURV1 <- pweibull(q = dat.train$y1, scale = mu1_hat, shape = sigma1_hat, lower.tail = FALSE) 
    PDF1 <- dweibull(x = dat.train$y1, scale = mu1_hat, shape = sigma1_hat)
    
    
    SURV2 <- 1 - pLogLogistic(x = dat.train$y2, mu = mu2_hat, sigma = sigma2_hat)
    PDF2 <- dLogLogistic(x = dat.train$y2, mu = mu2_hat, sigma = sigma2_hat)
    
    ################################################################################################################
    ################################################################################################################    FIT SEPARATELY: 
    
    
    cat("/n Fitting non-terminal event in two-steps /n")
    
    ### Margin 2 is ready. Fit margin 1: 
    Margin1_Formula_SOLO <- as.formula(paste("cbind(y1, delta1) ~", paste0(preds, collapse = "+")))
    
    GLM_M1_SOLO <- glmboostLSS(formula = Margin1_Formula_SOLO, 
                               data = dat.train, 
                               families = Custom_WeibullFamily(stabilization = "L2"), 
                               weights = weights.mstop, 
                               method = "noncyclic",
                               control = boost_control(mstop = TryMSTOP_Margins, 
                                                       nu = boost.nu.steplength,
                                                       risk = "oobag", 
                                                       trace = TRUE))
    
    risk_duringtraining_M1_SOLO <- risk(GLM_M1_SOLO, merge = TRUE)
    
    MSTOP_OPT_M1_SOLO <- which.min(risk(GLM_M1_SOLO, merge = TRUE))
    
    MSTOP_OPT_M1_SOLO
    
    rm(GLM_M1_SOLO)
    
    GLM_M1_SOLO <- glmboostLSS(formula = Margin1_Formula_SOLO, 
                               data = dat.train, 
                               families = Custom_WeibullFamily(stabilization = "L2"), 
                               method = "noncyclic",
                               control = boost_control(mstop = MSTOP_OPT_M1_SOLO,
                                                       nu = boost.nu.steplength,
                                                       trace = TRUE))
    
    
    coef(GLM_M1_SOLO)
    
    risk_final_M1_SOLO <- risk(GLM_M1_SOLO, merge = TRUE)
    
    
    #### Dependence model:
    ## Obtain survival function and pdf from non-terminal event: 
    mu1_hat_solo <- exp( predict(GLM_M1_SOLO$mu, type = "link") )
    sigma1_hat_solo <- exp( predict(GLM_M1_SOLO$sigma, type = "link") )
    
    SURV1_SOLO <- pweibull(q = dat.train$y1, scale = mu1_hat_solo, shape = sigma1_hat_solo, lower.tail = FALSE) 
    PDF1_SOLO <- dweibull(x = dat.train$y1, scale = mu1_hat_solo, shape = sigma1_hat_solo)
    
    
    
    ### Create new dataset for DEPENDENCE PARAMETER
    xmat_coppar <- x.train
    yvec_coppar <- cbind(SURV1_SOLO, PDF1_SOLO, dat.train$delta1,
                         SURV2, PDF2, dat.train$delta2)
    
    cat("Fitting dependence parameter (simple two-step)...\n")
    
    GLM_COPPAR <- glmboost(y = yvec_coppar,
                           x = xmat_coppar,
                           family = BivAFT_RC_GumbelCopula_RhoSoloFamily(),
                           weights = weights.mstop,
                           control = boost_control(mstop = (TryMSTOP_Dependence + 250),
                                                   nu = boost.nu.steplength,
                                                   risk = "oobag",
                                                   trace = TRUE),
                           center = FALSE)
    
    risk_duringtraining_DEPENDENCE <- risk(GLM_COPPAR)
    
    MSTOP_OPT_COPPAR <- which.min(risk(GLM_COPPAR))
    
    MSTOP_OPT_COPPAR
    
    rm(GLM_COPPAR)
    
    GLM_COPPAR <- glmboost(y = yvec_coppar,
                           x = xmat_coppar,
                           family = BivAFT_RC_GumbelCopula_RhoSoloFamily(),
                           control = boost_control(mstop = MSTOP_OPT_COPPAR,
                                                   nu = boost.nu.steplength,
                                                   trace = TRUE),
                           center = FALSE)
    
    coef(GLM_COPPAR)
    
    risk_final_DEPENDENCE <- risk(GLM_COPPAR)
    
    ############################################################################################# Margin 2 with COX:
    cat("Fitting terminal event (margin 2) using Cox model...\n")
    
    library(survival)
    
    if(p > 500){
      
      TryMSTOP_Margins_COX <- TryMSTOP_Margins + 1300
      
    }else{
      
      TryMSTOP_Margins_COX <- TryMSTOP_Margins + 1500
      
    }
    
    
    Margin2_Formula_COX <- as.formula(paste("Surv(y2, delta2) ~ ", paste0(preds, collapse = "+")))
    
    GLM_COX_M2 <- glmboost(y = Surv(dat.train$y2, dat.train$delta2), 
                           x = x.train[,-1], 
                           family = CoxPH(), 
                           weights = weights.mstop, 
                           control = boost_control(mstop = TryMSTOP_Margins_COX, 
                                                   nu = boost.nu.steplength,
                                                   risk = "oobag", 
                                                   trace = TRUE), 
                           center = FALSE)
    
    risk_duringtraining_M2_COX <- risk(GLM_COX_M2, merge = TRUE)
    
    MSTOP_OPT_M2_COX <- which.min(risk(GLM_COX_M2, merge = TRUE))
    
    MSTOP_OPT_M2_COX
    
    rm(GLM_COX_M2)
    
    GLM_COX_M2 <- glmboost(y = Surv(dat.train$y2, dat.train$delta2),
                           x = x.train[,-1],
                           family = CoxPH(), 
                           control = boost_control(mstop = MSTOP_OPT_M2_COX,
                                                   nu = boost.nu.steplength,
                                                   trace = TRUE),
                           center = FALSE)
    
    risk_final_M2_COX <- risk(GLM_COX_M2)
    
    
    SURVFIT_THING_MBOOST_TRAIN <- mboost::survFit(GLM_COX_M2)
    
    
    ################################################################################################################ Margin 2 with COX:
    cat("Fitting non-terminal event (margin 1) using Cox model...\n")
    
    TryMSTOP_Margins_COX <- TryMSTOP_Margins_COX 
    
    Margin1_Formula_COX <- as.formula(paste("Surv(y1, delta1) ~ ", paste0(preds, collapse = "+")))
    
    GLM_COX_M1 <- glmboost(y = Surv(dat.train$y1, dat.train$delta1), 
                           x = x.train[,-1], 
                           family = CoxPH(), 
                           weights = weights.mstop, 
                           control = boost_control(mstop = TryMSTOP_Margins_COX, 
                                                   nu = boost.nu.steplength,
                                                   risk = "oobag", 
                                                   trace = TRUE), 
                           center = FALSE)
    
    risk_duringtraining_M1_COX <- risk(GLM_COX_M1, merge = TRUE)
    
    MSTOP_OPT_M1_COX <- which.min(risk(GLM_COX_M1, merge = TRUE))
    
    MSTOP_OPT_M1_COX
    
    rm(GLM_COX_M1)
    
    GLM_COX_M1 <- glmboost(y = Surv(dat.train$y1, dat.train$delta1), 
                           x = x.train[,-1], 
                           family = CoxPH(), 
                           control = boost_control(mstop = MSTOP_OPT_M1_COX,
                                                   nu = boost.nu.steplength,
                                                   trace = TRUE),
                           center = FALSE)
    
    risk_final_M1_COX <- risk(GLM_COX_M1)
    
    SURVFIT_THING_MBOOST_TRAIN_M1 <- mboost::survFit(GLM_COX_M1)
    
    
    ################################################################################################################
    ################################################################################################################ Compute stuff
    ## Losses: 
    ListOfLOSS <- list()
    
    dat.test$Xint <- rep(1, nrow(dat.test))
    
    x.test_coppar <- cbind(rep(1, nrow(dat.test)), x.test)
    
    colnames(x.test_coppar)[1] <- "Xint"
    
    # Predict on test data:
    mu1_hat_OnTest      <- exp( predict(GLM_M1$mu1, dat.test, type = "link") )
    sigma1_hat_OnTest   <- exp( predict(GLM_M1$sigma1, dat.test, type = "link") )
    
    mu2_hat_OnTest      <- exp( predict(GLM_M2$mu, dat.test, type = "link") )
    sigma2_hat_OnTest   <- exp( predict(GLM_M2$sigma, dat.test, type = "link") )
    
    rho_hat_OnTest      <- exp( predict(GLM_M1$rho, dat.test, type = "link") ) + 1
    
    
    ######
    mu1_hat_OnTest_SOLO      <- exp( predict(GLM_M1_SOLO$mu, dat.test, type = "link") )
    sigma1_hat_OnTest_SOLO   <- exp( predict(GLM_M1_SOLO$sigma, dat.test, type = "link") )
    rho_hat_OnTest_SOLO    <- exp( predict(GLM_COPPAR, x.test_coppar, type = "link") ) + 1
    
    M1_LOSS   <- loss_M1(mu1 = mu1_hat_OnTest, sigma1 = sigma1_hat_OnTest, y1 = dat.test$y1, delta1 = dat.test$delta1)
    M2_LOSS   <- loss_M2(mu2 = mu2_hat_OnTest, sigma2 = sigma2_hat_OnTest, y2 = dat.test$y2, delta2 = dat.test$delta2)
    COP_LOSS  <- loss(mu1 = mu1_hat_OnTest, sigma1 = sigma1_hat_OnTest, 
                      mu2 = mu2_hat_OnTest, sigma2 = sigma2_hat_OnTest, 
                      delta1 = dat.test$delta1,
                      delta2 = dat.test$delta2,
                      rho = rho_hat_OnTest,
                      y = cbind(dat.test$y1, dat.test$y2), 
                      FAM = COPFAM
    ) 
    
    COP_LOSS_TSE <- loss(mu1 = mu1_hat_OnTest_SOLO, sigma1 = sigma1_hat_OnTest_SOLO, 
                         mu2 = mu2_hat_OnTest, sigma2 = sigma2_hat_OnTest, 
                         delta1 = dat.test$delta1,
                         delta2 = dat.test$delta2,
                         rho = rho_hat_OnTest_SOLO,
                         y = cbind(dat.test$y1, dat.test$y2), 
                         FAM = COPFAM
    ) 
    
    ListOfLOSS <- list(Univariate = (M1_LOSS + M2_LOSS), 
                       UnivariateCOX = (GLM_COX_M1$logLik() + GLM_COX_M2$logLik()),
                       Copula = COP_LOSS,
                       Copula_TWOSTAGE = COP_LOSS_TSE)
    
    
    ################## Survival functions and PDF on TEST
    SURV1_ONTEST_FORCONC <- pweibull(q = dat.test$y1, scale = mu1_hat_OnTest, shape = sigma1_hat_OnTest, lower.tail = FALSE) 
    PDF1_ONTEST_FORCONC <- dweibull(x = dat.test$y1, scale = mu1_hat_OnTest, shape = sigma1_hat_OnTest)
    
    
    SURV2_ONTEST_FORCONC <- 1 - pLogLogistic(x = dat.test$y2, mu = mu2_hat_OnTest, sigma = sigma2_hat_OnTest)
    PDF2_ONTEST_FORCONC <- dLogLogistic(x = dat.test$y2, mu = mu2_hat_OnTest, sigma = sigma2_hat_OnTest)
    
    
    
    
    #####################################################################################################################################
    ## Scores and stuff
    
    # Create grid of time values: 
    time_grid_M1 <- quantile(dat.test$y1, probs = seq(0.05, 0.95, by = 0.05))
    time_grid_M2 <- quantile(dat.test$y2, probs = seq(0.05, 0.95, by = 0.05)) 
    
    time_grid_matrix_M1 <- matrix(replicate(time_grid_M1, n = nrow(dat.test)), nrow = length(time_grid_M1), ncol = nrow(dat.test), byrow = FALSE)
    time_grid_matrix_M2 <- matrix(replicate(time_grid_M2, n = nrow(dat.test)), nrow = length(time_grid_M2), ncol = nrow(dat.test), byrow = FALSE)
    
    KM_Object_M1 <- survival::survfit(survival::Surv(y1, (1 - delta1) ) ~ 1, data = dat.train)
    KM_Object_M2 <- survival::survfit(survival::Surv(y2, (1 - delta2) ) ~ 1, data = dat.train)
    
    # Predict the survival probability at the time_grid
    ProbOfCens_M1 <- summary(KM_Object_M1, times = time_grid_M1, extend = TRUE)$surv
    ProbOfCens_M2 <- summary(KM_Object_M2, times = time_grid_M2, extend = TRUE)$surv
    
    
    ProbOfCens_Test_M1 <- summary(KM_Object_M1, times = dat.test$y1, extend = TRUE)$surv
    ProbOfCens_Test_M2 <- summary(KM_Object_M2, times = dat.test$y2, extend = TRUE)$surv
    
    ### We need to compute: S_hat( time_star ), G_hat( time_star ), G_hat_Time... (all using the test observations....)
    M1_S_hat_star <- t(sapply(1:length(time_grid_M1), function(i) 
      1 - pweibull(q = time_grid_matrix_M1[i,], scale = mu1_hat_OnTest, shape = sigma1_hat_OnTest), simplify = TRUE))
    
    M2_S_hat_star <- t(sapply(1:length(time_grid_M2), function(i) 
      1 - pfisk(q = time_grid_matrix_M2[i,], scale = mu2_hat_OnTest, shape1.a = sigma2_hat_OnTest), simplify = TRUE))
    
    
    #### Compute this too using the TWO-STAGE-MODEL (SOLO)
    M1_S_hat_star_SOLO <- t(sapply(1:length(time_grid_M1), function(i) 
      1 - pweibull(q = time_grid_matrix_M1[i,], scale = mu1_hat_OnTest_SOLO, shape = sigma1_hat_OnTest_SOLO), simplify = TRUE))
    
    
    ### Baseline fit using survfit to obtain Cox predictions: 
    M1_Baseline_Survival_COX <- survival::survfit(survival::Surv(y1, delta1) ~ 1, data = dat.train) 
    M2_Baseline_Survival_COX <- survival::survfit(survival::Surv(y2, delta2) ~ 1, data = dat.train) 
    
    M1_BaselineCUMUHAZ_TestTimes_COX <- summary(M1_Baseline_Survival_COX, times = time_grid_M1, extend = TRUE, cumhaz = TRUE)$cumhaz
    M2_BaselineCUMUHAZ_TestTimes_COX <- summary(M2_Baseline_Survival_COX, times = time_grid_M2, extend = TRUE, cumhaz = TRUE)$cumhaz
    
    M1_COX_ResponsePreds <- exp( predict(GLM_COX_M1, newdata = x.test, type = "link") )
    M2_COX_ResponsePreds <- exp( predict(GLM_COX_M2, newdata = x.test, type = "link") )
    
    M1_S_hat_star_COX <- matrix(0, ncol = nrow(dat.test), nrow = length(time_grid_M1))
    M2_S_hat_star_COX <- matrix(0, ncol = nrow(dat.test), nrow = length(time_grid_M2))
    
    ## For the Cox model, the predicted survival probability is: S( t ) = exp( - CumuHaz(t) * exp( x beta_hat ) )
    for(i in 1:length(time_grid_M1)){
      
      M1_S_hat_star_COX[i,] <- exp( - M1_Baseline_Survival_COX$cumhaz[i] * M1_COX_ResponsePreds )
      
      M2_S_hat_star_COX[i,] <- exp( - M2_Baseline_Survival_COX$cumhaz[i] * M2_COX_ResponsePreds )
      
    }
    
    
    
    M1_G_hat_star <- matrix(replicate(ProbOfCens_M1, n = nrow(dat.test)), nrow = length(time_grid_M1), ncol = nrow(dat.test), byrow = FALSE)
    M2_G_hat_star <- matrix(replicate(ProbOfCens_M2, n = nrow(dat.test)), nrow = length(time_grid_M2), ncol = nrow(dat.test), byrow = FALSE)
    
    M1_G_hat_test <- ProbOfCens_Test_M1
    M2_G_hat_test <- ProbOfCens_Test_M2
    
    BrierScores_M1 <- sapply(1:length(time_grid_M1), function(i) 
      ComputeBrierScore(S_hat = M1_S_hat_star[i,], 
                        G_hat_atT = M1_G_hat_test, 
                        G_hat_atTSTAR = M1_G_hat_star[i,], 
                        delta = dat.test$delta1, 
                        time_star = time_grid_matrix_M1[i,], 
                        time_i = dat.test$y1), 
      simplify = TRUE)
    
    BrierScores_M2 <- sapply(1:length(time_grid_M2), function(i) 
      ComputeBrierScore(S_hat = M2_S_hat_star[i,], 
                        G_hat_atT = M2_G_hat_test, 
                        G_hat_atTSTAR = M2_G_hat_star[i,], 
                        delta = dat.test$delta2, 
                        time_star = time_grid_matrix_M2[i,], 
                        time_i = dat.test$y2), 
      simplify = TRUE)
    
    
    ##### USING TWO-STEP:
    BrierScores_M1_SOLO <- sapply(1:length(time_grid_M1), function(i) 
      ComputeBrierScore(S_hat = M1_S_hat_star_SOLO[i,], 
                        G_hat_atT = M1_G_hat_test, 
                        G_hat_atTSTAR = M1_G_hat_star[i,], 
                        delta = dat.test$delta1, 
                        time_star = time_grid_matrix_M1[i,], 
                        time_i = dat.test$y1), 
      simplify = TRUE)
    
    
    BrierScores_M1_COX <- sapply(1:length(time_grid_M1), function(i) 
      ComputeBrierScore(S_hat = M1_S_hat_star_COX[i,], 
                        G_hat_atT = M1_G_hat_test, 
                        G_hat_atTSTAR = M1_G_hat_star[i,], 
                        delta = dat.test$delta1, 
                        time_star = time_grid_matrix_M1[i,], 
                        time_i = dat.test$y1), 
      simplify = TRUE)
    
    BrierScores_M2_COX <- sapply(1:length(time_grid_M2), function(i) 
      ComputeBrierScore(S_hat = M2_S_hat_star_COX[i,], 
                        G_hat_atT = M2_G_hat_test, 
                        G_hat_atTSTAR = M2_G_hat_star[i,], 
                        delta = dat.test$delta2, 
                        time_star = time_grid_matrix_M2[i,], 
                        time_i = dat.test$y2), 
      simplify = TRUE)
    
    
    ## Integrated brier: 
    ComputeBrierScore_ForIntegration <- function(x, time_i, 
                                                 SurvObject,
                                                 delta,  
                                                 muhat, sigmahat, margin = 1){
      
      x_star <- x
      
      
      
      # Obtain Survival probability at 
      G_hat_atT     <- summary(SurvObject, times = time_i, extend = TRUE)$surv
      G_hat_atTSTAR <- summary(SurvObject, times = x_star, extend = TRUE)$surv
      
      S_hat <- vector(mode = "numeric", length = length(muhat))
      
      for( i in 1:length(muhat)){
        
        if(margin == 1){
          
          S_hat[i] <- 1 - pweibull(q = x_star, scale = muhat[i], shape = sigmahat[i])
          
        }else{
          
          S_hat[i] <- 1 - pfisk(q = x_star, scale = muhat[i], shape1.a = sigmahat[i])
          
        }
      }
      
      
      S_hat <- pdffz(S_hat)
      G_hat_atT <- pdffz(G_hat_atT)
      G_hat_atTSTAR <- pdffz(G_hat_atTSTAR)
      
      Term1 <- ( S_hat )^2 / ( G_hat_atT )
      
      Term2 <- ( 1 - S_hat )^2 / ( G_hat_atTSTAR ) 
      
      Case1 <- as.numeric( ( time_i < x_star ) & delta == 1)
      
      Case2 <- as.numeric( (time_i >= x_star ) )
      
      IBS <- mean(Term1*Case1 + Term2*Case2)
      
      ## NEW CHECK
      IBS <- ifelse(IBS > 1, 1, IBS)
      
      return(IBS)
    }
    
    
    ComputeBrierScore_ForIntegration_fromCOX <- function(x, time_i, 
                                                         SurvObject,
                                                         delta, 
                                                         baselineCumuHazardObject,
                                                         CoxHazardsPreds){
      
      x_star <- x
      
      BaselineCUMUHAZ <- summary(baselineCumuHazardObject, times = x_star, extend = TRUE, cumhaz = TRUE)$cumhaz
      
      # Obtain Survival probability at 
      G_hat_atT     <- summary(SurvObject, times = time_i, extend = TRUE)$surv
      G_hat_atTSTAR <- summary(SurvObject, times = x_star, extend = TRUE)$surv
      
      S_hat <- vector(mode = "numeric", length = length(delta))
      
      for( i in 1:length(delta)){
        
        S_hat[i] <- exp( - BaselineCUMUHAZ * M1_COX_ResponsePreds[i] )
        
      }
      
      S_hat <- pdffz(S_hat)
      G_hat_atT <- pdffz(G_hat_atT)
      G_hat_atTSTAR <- pdffz(G_hat_atTSTAR)
      
      Term1 <- ( S_hat )^2 / ( G_hat_atT )
      
      Term2 <- ( 1 - S_hat )^2 / ( G_hat_atTSTAR ) 
      
      Case1 <- as.numeric( ( time_i < x_star ) & delta == 1)
      
      Case2 <- as.numeric( (time_i >= x_star ) )
      
      IBS <- mean(Term1*Case1 + Term2*Case2)
      
      ## NEW CHECK
      IBS <- ifelse(IBS > 1, 1, IBS)
      
      return(IBS)
      
    }
    
    IntegratedBrier_M1 <- integrate(f = Vectorize(ComputeBrierScore_ForIntegration, vectorize.args = "x"), 
                                    lower = 0, upper = max(time_grid_M1), 
                                    margin = 1, 
                                    muhat = mu1_hat_OnTest, 
                                    sigmahat = sigma1_hat_OnTest, 
                                    delta = dat.test$delta1,
                                    time_i = dat.test$y1,
                                    SurvObject = KM_Object_M1, 
                                    subdivisions = 1000)
    
    ## 1 / t_max, according to definition:
    IntegratedBrier_M1_Def <- 1 / max(time_grid_M1) * IntegratedBrier_M1$value
    
    
    
    
    IntegratedBrier_M2 <- integrate(f = Vectorize(ComputeBrierScore_ForIntegration, vectorize.args = "x"), 
                                    lower = 0, upper = max(time_grid_M2), 
                                    margin = 2, 
                                    muhat = mu2_hat_OnTest, 
                                    sigmahat = sigma2_hat_OnTest, 
                                    delta = dat.test$delta2,
                                    time_i = dat.test$y2,
                                    SurvObject = KM_Object_M2,
                                    subdivisions = 1250)
    
    ## 1 / t_max, according to definition:
    IntegratedBrier_M2_Def <- 1 / max(time_grid_M2) * IntegratedBrier_M2$value
    
    #################################################################################################################### USING SOLO:
    IntegratedBrier_M1_SOLO <- integrate(f = Vectorize(ComputeBrierScore_ForIntegration, vectorize.args = "x"), 
                                         lower = 0, upper = max(time_grid_M1), 
                                         margin = 1, 
                                         muhat = mu1_hat_OnTest_SOLO, 
                                         sigmahat = sigma1_hat_OnTest_SOLO, 
                                         delta = dat.test$delta1,
                                         time_i = dat.test$y1,
                                         SurvObject = KM_Object_M1, 
                                         subdivisions = 1000)
    
    ## 1 / t_max, according to definition:
    IntegratedBrier_M1_Def_SOLO <- 1 / max(time_grid_M1) * IntegratedBrier_M1_SOLO$value
    
    
    
    ############ USING COX MODELS: 
    IntegratedBrier_M1_COX <- integrate(f = Vectorize(ComputeBrierScore_ForIntegration_fromCOX, vectorize.args = "x"), 
                                        lower = 0, upper = max(time_grid_M1), 
                                        delta = dat.test$delta1,
                                        time_i = dat.test$y1,
                                        SurvObject = KM_Object_M1,
                                        baselineCumuHazardObject = M1_Baseline_Survival_COX,
                                        CoxHazardsPreds = M1_COX_ResponsePreds,
                                        subdivisions = 1000)
    
    ## 1 / t_max, according to definition:
    IntegratedBrier_M1_COX_Def <- 1 / max(time_grid_M1) * IntegratedBrier_M1_COX$value
    
    
    IntegratedBrier_M2_COX <- integrate(f = Vectorize(ComputeBrierScore_ForIntegration_fromCOX, vectorize.args = "x"), 
                                        lower = 0, upper = max(time_grid_M1), 
                                        delta = dat.test$delta2,
                                        time_i = dat.test$y2,
                                        SurvObject = KM_Object_M2,
                                        baselineCumuHazardObject = M2_Baseline_Survival_COX,
                                        CoxHazardsPreds = M2_COX_ResponsePreds,
                                        subdivisions = 1000)
    
    ## 1 / t_max, according to definition:
    IntegratedBrier_M2_COX_Def <- 1 / max(time_grid_M2) * IntegratedBrier_M2_COX$value
    
    
    #####################################################################################################################################
    # Concordance: 
    
    
    M1_Shat_OnTest <- 1 - pweibull(q = dat.test$y1, scale = mu1_hat_OnTest, shape = sigma1_hat_OnTest)
    M2_Shat_OnTest <- 1 - pfisk(q = dat.test$y2, scale = mu2_hat_OnTest, shape1.a = sigma2_hat_OnTest)
    
    ################ SOLO
    M1_Shat_OnTest_SOLO <- 1 - pweibull(q = dat.test$y1, scale = mu1_hat_OnTest_SOLO, shape = sigma1_hat_OnTest_SOLO)
    
    M1_ShatBaseline_OnTest_COX <- summary(M1_Baseline_Survival_COX, times = dat.test$y1, extend = TRUE, cumhaz = TRUE)
    M2_ShatBaseline_OnTest_COX <- summary(M2_Baseline_Survival_COX, times = dat.test$y2, extend = TRUE, cumhaz = TRUE)
    
    M1_ShatBaseline_OnTest_COX <- M1_ShatBaseline_OnTest_COX$cumhaz
    M2_ShatBaseline_OnTest_COX <- M2_ShatBaseline_OnTest_COX$cumhaz 
    
    M1_Shat_OnTest_COX <- exp( - M1_ShatBaseline_OnTest_COX * M1_COX_ResponsePreds )
    M2_Shat_OnTest_COX <- exp( - M2_ShatBaseline_OnTest_COX * M2_COX_ResponsePreds )
    
    
    CIND_M1 <- survcomp:::concordance.index(x = M1_Shat_OnTest, surv.time = dat.test$y1, surv.event = dat.test$delta1)
    CIND_M2 <- survcomp:::concordance.index(x = M2_Shat_OnTest, surv.time = dat.test$y2, surv.event = dat.test$delta2)
    
    ####### SOLO
    CIND_M1_SOLO <- survcomp:::concordance.index(x = M1_Shat_OnTest_SOLO, surv.time = dat.test$y1, surv.event = dat.test$delta1)
    
    CIND_M1_COX <- survcomp:::concordance.index(x = M1_Shat_OnTest_COX, surv.time = dat.test$y1, surv.event = dat.test$delta1)
    CIND_M2_COX <- survcomp:::concordance.index(x = M2_Shat_OnTest_COX, surv.time = dat.test$y2, surv.event = dat.test$delta2)
    
    
    
    CONCORDANCE_M1 <- list(CINDEX = CIND_M1$c.index, 
                           SE = CIND_M1$se,
                           INTERVAL = c(CIND_M1$lower, CIND_M1$upper),
                           PVALUE = CIND_M1$p.value,
                           COMPPAIRS = CIND_M1$comppairs
    )
    
    
    ##### SOLO
    CONCORDANCE_M1_SOLO <- list(CINDEX = CIND_M1_SOLO$c.index, 
                                SE = CIND_M1_SOLO$se,
                                INTERVAL = c(CIND_M1_SOLO$lower, CIND_M1$upper),
                                PVALUE = CIND_M1_SOLO$p.value,
                                COMPPAIRS = CIND_M1_SOLO$comppairs
    )
    
    CONCORDANCE_M2 <- list(CINDEX = CIND_M2$c.index, 
                           SE = CIND_M2$se,
                           INTERVAL = c(CIND_M2$lower, CIND_M2$upper),
                           PVALUE = CIND_M2$p.value,
                           COMPPAIRS = CIND_M2$comppairs
    )
    
    CONCORDANCE_M1_COX <- list(CINDEX = CIND_M1_COX$c.index, 
                               SE = CIND_M1_COX$se,
                               INTERVAL = c(CIND_M1_COX$lower, CIND_M1_COX$upper),
                               PVALUE = CIND_M1_COX$p.value,
                               COMPPAIRS = CIND_M1_COX$comppairs
    )
    
    CONCORDANCE_M2_COX <- list(CINDEX = CIND_M2_COX$c.index, 
                               SE = CIND_M2_COX$se,
                               INTERVAL = c(CIND_M2_COX$lower, CIND_M2_COX$upper),
                               PVALUE = CIND_M2_COX$p.value,
                               COMPPAIRS = CIND_M2_COX$comppairs
    )
    
    
    
    ### Normal concordance index:
    M1_SURV_CONCORDANCE <- list(CINDEX = survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ mu1_hat_OnTest)$concordance,
                                CINDEX_SE = sqrt(survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ mu1_hat_OnTest)$var),
                                #
                                #
                                CINDEX_UNO = survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ mu1_hat_OnTest, timewt = "n/G2")$concordance,
                                CINDEX_UNO_SE = sqrt(survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ mu1_hat_OnTest, timewt = "n/G2")$var),
                                #
                                #
                                # SOLO ESTIMATION:
                                CINDEX_SOLO = survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ mu1_hat_OnTest_SOLO)$concordance,
                                CINDEX_SE_SOLO = sqrt(survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ mu1_hat_OnTest_SOLO)$var),
                                #
                                # 
                                CINDEX_UNO_SOLO = survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ mu1_hat_OnTest_SOLO, timewt = "n/G2")$concordance,
                                CINDEX_UNO_SE_SOLO = sqrt(survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ mu1_hat_OnTest_SOLO, timewt = "n/G2")$var),
                                #
                                #
                                CINDEX_COX = survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ M1_COX_ResponsePreds)$concordance,
                                CINDEX_SE_COX = sqrt(survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ M1_COX_ResponsePreds)$var),
                                #
                                #
                                #
                                CINDEX_UNO_COX = survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ M1_COX_ResponsePreds, timewt = "n/G2")$concordance,
                                CINDEX_UNO_SE_COX = sqrt(survival::concordance(Surv(dat.test$y1, dat.test$delta1) ~ M1_COX_ResponsePreds, timewt = "n/G2")$var)
    )
    
    
    
    M2_SURV_CONCORDANCE <- list(CINDEX = survival::concordance(Surv(dat.test$y2, dat.test$delta2) ~ mu2_hat_OnTest)$concordance,
                                CINDEX_SE = sqrt(survival::concordance(Surv(dat.test$y2, dat.test$delta2) ~ mu2_hat_OnTest)$var),
                                #
                                #
                                CINDEX_UNO = survival::concordance(Surv(dat.test$y2, dat.test$delta2) ~ mu2_hat_OnTest, timewt = "n/G2")$concordance,
                                CINDEX_UNO_SE = sqrt(survival::concordance(Surv(dat.test$y2, dat.test$delta2) ~ mu2_hat_OnTest, timewt = "n/G2")$var),
                                #
                                #
                                CINDEX_COX = survival::concordance(Surv(dat.test$y2, dat.test$delta2) ~ M2_COX_ResponsePreds)$concordance,
                                CINDEX_SE_COX = sqrt(survival::concordance(Surv(dat.test$y2, dat.test$delta2) ~ M2_COX_ResponsePreds)$var),
                                #
                                #
                                #
                                CINDEX_UNO_COX = survival::concordance(Surv(dat.test$y2, dat.test$delta2) ~ M2_COX_ResponsePreds, timewt = "n/G2")$concordance,
                                CINDEX_UNO_SE_COX = sqrt(survival::concordance(Surv(dat.test$y2, dat.test$delta2) ~ M2_COX_ResponsePreds, timewt = "n/G2")$var)
    )
    
    
    
    
    
    
    
    
    
    
    
    #####################################################################################################################################
    # Integrated Absolute and Squared Error
    Diff_SurvFunctions <- function(x, muhat, sigmahat, muTRUE, sigmaTRUE, margin = 1, type = "abs"){
      
      if(margin == 1){
        
        STrue <- 1 - pweibull(q = x, scale = muTRUE, shape = sigmaTRUE)
        Shat  <- 1 - pweibull(q = x, scale = muhat, shape = sigmahat)
        
      }else{
        
        STrue <- 1 - pfisk(q = x, scale = muTRUE, shape1.a = sigmaTRUE)
        Shat  <- 1 - pfisk(q = x, scale = muhat, shape1.a = sigmahat)
      }
      
      
      STrue <- pdffz(STrue)
      Shat <- pdffz(Shat)
      
      if(type == "abs"){
        
        diff <- abs(STrue - Shat)
        
      }else{
        
        diff <- (STrue - Shat)^2
        
      }
      
      
      
      return(diff)
    }
    
    Diff_SurvFunctions_COX <- function(x, baselineObject, hazPred, muTRUE, sigmaTRUE, margin = 1, type = "abs"){
      
      baseCumu <-  summary(baselineObject, times = x, extend = TRUE, cumhaz = TRUE)
      
      baseCumu <- baseCumu$cumhaz
      
      if(margin == 1){
        
        STrue <- 1 - pweibull(q = x, scale = muTRUE, shape = sigmaTRUE)
        Shat  <- exp( - baseCumu * hazPred )
        
      }else{
        
        STrue <- 1 - pfisk(q = x, scale = muTRUE, shape1.a = sigmaTRUE)
        Shat  <- exp( - baseCumu * hazPred )
      }
      
      
      
      STrue <- pdffz(STrue)
      Shat <- pdffz(Shat)
      
      if(type == "abs"){
        
        diff <- abs(STrue - Shat)
        
      }else{
        
        diff <- (STrue - Shat)^2
        
      }
      
      
      
      return(diff)
    }
    
    M1_IAE <- t(sapply(1:length(mu1_test),function(i) integrate(f = Diff_SurvFunctions, 
                                                                lower = 0,
                                                                upper = max(time_grid_M1), 
                                                                muhat = mu1_hat_OnTest[i], sigmahat = sigma1_hat_OnTest[i],
                                                                muTRUE = mu1_test[i], sigmaTRUE = sigma1_test[i],
                                                                margin = 1), 
                       simplify = TRUE))
    
    M1_ISE <- t(sapply(1:length(mu1_test),function(i) integrate(f = Diff_SurvFunctions, 
                                                                lower = 0,
                                                                upper = max(time_grid_M1), 
                                                                muhat = mu1_hat_OnTest[i], sigmahat = sigma1_hat_OnTest[i],
                                                                muTRUE = mu1_test[i], sigmaTRUE = sigma1_test[i], 
                                                                type = "squared",
                                                                margin = 1), 
                       simplify = TRUE))
    
    ############################ SOLO
    M1_IAE_SOLO <- t(sapply(1:length(mu1_test),function(i) integrate(f = Diff_SurvFunctions, 
                                                                     lower = 0,
                                                                     upper = max(time_grid_M1), 
                                                                     muhat = mu1_hat_OnTest_SOLO[i], sigmahat = sigma1_hat_OnTest_SOLO[i],
                                                                     muTRUE = mu1_test[i], sigmaTRUE = sigma1_test[i],
                                                                     margin = 1), 
                            simplify = TRUE))
    
    M1_ISE_SOLO <- t(sapply(1:length(mu1_test),function(i) integrate(f = Diff_SurvFunctions, 
                                                                     lower = 0,
                                                                     upper = max(time_grid_M1), 
                                                                     muhat = mu1_hat_OnTest_SOLO[i], sigmahat = sigma1_hat_OnTest_SOLO[i],
                                                                     muTRUE = mu1_test[i], sigmaTRUE = sigma1_test[i], 
                                                                     type = "squared",
                                                                     margin = 1), 
                            simplify = TRUE))
    
    
    
    M2_IAE <- t(sapply(1:length(mu1_test),function(i) integrate(f = Diff_SurvFunctions, 
                                                                lower = 0,
                                                                upper = max(time_grid_M2), 
                                                                muhat = mu2_hat_OnTest[i], sigmahat = sigma2_hat_OnTest[i],
                                                                muTRUE = mu2_test[i], sigmaTRUE = sigma2_test[i],
                                                                margin = 2), 
                       simplify = TRUE))
    
    M2_ISE <- t(sapply(1:length(mu1_test),function(i) integrate(f = Diff_SurvFunctions, 
                                                                lower = 0,
                                                                upper = max(time_grid_M2), 
                                                                muhat = mu2_hat_OnTest[i], sigmahat = sigma2_hat_OnTest[i],
                                                                muTRUE = mu2_test[i], sigmaTRUE = sigma2_test[i], 
                                                                type = "squared",
                                                                margin = 2), 
                       simplify = TRUE))
    
    
    
    ############################################################################################### INTEGRATION BY HAND
    
    # time grid is chosen based on the individual maximum time point seen in the test data.
    # should not make a difference because we are integrating out time.... 
    a <- 0
    b <- max(dat.test$y1) 
    
    b2 <- max(dat.test$y2) 
    
    # number of subintervals
    n_int <- 5000
    
    # width of each subinterval
    h <- (b - a) / n_int
    
    x <- seq(a, b, length.out = n_int + 1)
    
    x2 <- seq(a, b2, length.out = n_int + 1)
    
    ## Obtain some predictions from the cox models:
    Cox_ResponsePreds <- exp( predict(GLM_COX_M1, newdata = x.test, type = "link") )
    baseline_surv <- survival::survfit(survival::Surv(y1, delta1) ~ 1, data = dat.train) 
    baseline_surv_times <- summary(baseline_surv, times = x, extend = TRUE, cumhaz = TRUE)
    
    Cox_ResponsePreds_M2 <- exp( predict(GLM_COX_M2, newdata = x.test, type = "link") )
    baseline_surv_M2 <- survival::survfit(survival::Surv(y2, delta2) ~ 1, data = dat.train) 
    baseline_surv_times_M2 <- summary(baseline_surv_M2, times = x2, extend = TRUE, cumhaz = TRUE)
    
    
    #### Initialise
    val_b_cop_squared <- vector()
    val_b_cox_squared <- vector()
    
    val_b_cop_abs <- vector()
    val_b_cox_abs <- vector()
    
    val_b_cop_squared_M2 <- vector()
    val_b_cox_squared_M2 <- vector()
    
    val_b_cop_abs_M2 <- vector()
    val_b_cox_abs_M2 <- vector()
    
    ########## SOLO 
    val_b_cop_squared_SOLO <- vector()
    val_b_cop_abs_SOLO <- vector()
    
    
    
    for(i in 1:length(mu1_test)){
      
      
      ### MARGIN 1
      y_true <- sapply(x, function(x) pweibull(q = x, scale = mu1_test[i], shape = sigma1_test[i], lower.tail = FALSE))
      
      y_est  <- sapply(x, function(x) pweibull(q = x, scale = mu1_hat_OnTest[i], shape = sigma1_hat_OnTest[i], lower.tail = FALSE))
      
      y_est_cox  <- sapply(1:length(x), function(x) exp( - baseline_surv_times$cumhaz[x] * Cox_ResponsePreds[i] ))
      
      
      ### MARGIN 1 SOLO:
      y_est_SOLO  <- sapply(x, function(x) pweibull(q = x, scale = mu1_hat_OnTest_SOLO[i], shape = sigma1_hat_OnTest_SOLO[i], lower.tail = FALSE))
      
      
      
      
      #### MARGIN 2
      y_true_M2 <- sapply(x2, function(x) pfisk(q = x, scale = mu2_test[i], shape1.a = sigma2_test[i], lower.tail = FALSE))
      
      y_est_M2  <- sapply(x2, function(x) pfisk(q = x, scale = mu2_hat_OnTest[i], shape1.a = sigma2_hat_OnTest[i], lower.tail = FALSE))
      
      y_est_cox_M2  <- sapply(1:length(x2), function(x) exp( - baseline_surv_times_M2$cumhaz[x] * Cox_ResponsePreds_M2[i] ))
      
      
      
      
      ############################################# SQUARED ERROR
      squared_y_diff_cop <- (y_true - y_est)^2
      squared_y_diff_cox <- (y_true - y_est_cox)^2
      
      squared_y_diff_cop_M2 <- (y_true_M2 - y_est_M2)^2
      squared_y_diff_cox_M2 <- (y_true_M2 - y_est_cox_M2)^2
      
      ##### SOLO
      squared_y_diff_cop_SOLO <- (y_true - y_est_SOLO)^2
      
      
      integral_cop <- sum( (squared_y_diff_cop[-1] + squared_y_diff_cop[-length(squared_y_diff_cop)] ) * h / 2)
      integral_cox <- sum( (squared_y_diff_cox[-1] + squared_y_diff_cox[-length(squared_y_diff_cox)] ) * h / 2)
      
      integral_cop_M2 <- sum( (squared_y_diff_cop_M2[-1] + squared_y_diff_cop_M2[-length(squared_y_diff_cop_M2)] ) * h / 2)
      integral_cox_M2 <- sum( (squared_y_diff_cox_M2[-1] + squared_y_diff_cox_M2[-length(squared_y_diff_cox_M2)] ) * h / 2)
      
      
      #### SOLO
      integral_cop_SOLO <- sum( (squared_y_diff_cop_SOLO[-1] + squared_y_diff_cop_SOLO[-length(squared_y_diff_cop_SOLO)] ) * h / 2)
      
      
      ############################################# ABSOLUTE ERROR
      abs_y_diff_cop <- abs(y_true - y_est)
      abs_y_diff_cox <- abs(y_true - y_est_cox)
      
      abs_y_diff_cop_M2 <- abs(y_true_M2 - y_est_M2)
      abs_y_diff_cox_M2 <- abs(y_true_M2 - y_est_cox_M2)
      
      ##### SOLO
      abs_y_diff_cop_SOLO <- abs(y_true - y_est_SOLO)
      
      
      integral_cop_abs <- sum( (abs_y_diff_cop[-1] + abs_y_diff_cop[-length(abs_y_diff_cop)] ) * h / 2)
      integral_cox_abs <- sum( (abs_y_diff_cox[-1] + abs_y_diff_cox[-length(abs_y_diff_cox)] ) * h / 2)
      
      
      integral_cop_abs_M2 <- sum( (abs_y_diff_cop_M2[-1] + abs_y_diff_cop_M2[-length(abs_y_diff_cop_M2)] ) * h / 2)
      integral_cox_abs_M2 <- sum( (abs_y_diff_cox_M2[-1] + abs_y_diff_cox_M2[-length(abs_y_diff_cox_M2)] ) * h / 2)
      
      
      ######## SOLO 
      integral_cop_abs_SOLO <- sum( (abs_y_diff_cop_SOLO[-1] + abs_y_diff_cop_SOLO[-length(abs_y_diff_cop_SOLO)] ) * h / 2)
      
      
      val_b_cop_squared <- c(val_b_cop_squared, integral_cop)
      val_b_cox_squared <- c(val_b_cox_squared, integral_cox)
      
      val_b_cop_abs <- c(val_b_cop_abs, integral_cop_abs)
      val_b_cox_abs <- c(val_b_cox_abs, integral_cox_abs)
      
      
      ######## SOLO
      val_b_cop_squared_SOLO <- c(val_b_cop_squared_SOLO, integral_cop_SOLO)
      val_b_cop_abs_SOLO <- c(val_b_cop_abs_SOLO, integral_cop_abs_SOLO)
      
      
      val_b_cop_squared_M2 <- c(val_b_cop_squared_M2, integral_cop_M2)
      val_b_cox_squared_M2 <- c(val_b_cox_squared_M2, integral_cox_M2)
      
      val_b_cop_abs_M2 <- c(val_b_cop_abs_M2, integral_cop_abs_M2)
      val_b_cox_abs_M2 <- c(val_b_cox_abs_M2, integral_cox_abs_M2)
      
      
    }
    
    #########################################################################################################
    PerformanceMetrics <- list()
    
    PerformanceMetrics$BRIERSCORE <- list(MARGIN1 = BrierScores_M1,
                                          #
                                          MARGIN1_SOLO = BrierScores_M1_SOLO,
                                          #
                                          MARGIN1_COX = BrierScores_M1_COX,
                                          MARGIN2 = BrierScores_M2,
                                          MARGIN2_COX = BrierScores_M2_COX
    )
    
    PerformanceMetrics$INTEGRATED_BRIERSCORE <- list(MARGIN1 = IntegratedBrier_M1, 
                                                     MARGIN1_DEF = IntegratedBrier_M1_Def,
                                                     #
                                                     MARGIN1_SOLO = IntegratedBrier_M1_SOLO, 
                                                     MARGIN1_DEF_SOLO = IntegratedBrier_M1_Def_SOLO,
                                                     #
                                                     MARGIN1_COX = IntegratedBrier_M1_COX,
                                                     MARGIN1_COX_DEF = IntegratedBrier_M1_COX_Def,
                                                     MARGIN2 = IntegratedBrier_M2,
                                                     MARGIN2_DEF = IntegratedBrier_M2_Def,
                                                     MARGIN2_COX = IntegratedBrier_M2_COX,
                                                     MARGIN2_COX_DEF = IntegratedBrier_M2_COX_Def
    )
    
    PerformanceMetrics$INTEGRATED_ABSOLUTE <- list(MARGIN1 = M1_IAE,
                                                   MARGIN1_SOLO = M1_IAE_SOLO,
                                                   MARGIN2 = M2_IAE)
    
    PerformanceMetrics$INTEGRATED_SQUARE <- list(MARGIN1 = M1_ISE,
                                                 MARGIN1_SOLO = M1_ISE,
                                                 MARGIN2 = M2_ISE)
    
    PerformanceMetrics$CONCORDANCE <- list(MARGIN1 = CONCORDANCE_M1, 
                                           MARGIN1_SOLO = CONCORDANCE_M1_SOLO, 
                                           MARGIN1_COX = CONCORDANCE_M1_COX,
                                           MARGIN2 = CONCORDANCE_M2,
                                           MARGIN2_COX = CONCORDANCE_M2_COX
    )
    
    PerformanceMetrics$SURV_CONCORDANCE <- list(MARGIN1 = M1_SURV_CONCORDANCE, 
                                                MARGIN2 = M2_SURV_CONCORDANCE
    )
    
    
    PerformanceMetrics$INTEGRATED_SQUARED_ERRORS_BYHAND <- list(MARGIN1 = mean(val_b_cop_squared),
                                                                MARGIN1 = mean(val_b_cop_squared_SOLO),
                                                                MARGIN1_COX = mean(val_b_cox_squared),
                                                                MARGIN2 =  mean(val_b_cop_squared_M2),
                                                                MARGIN2_COX = mean(val_b_cox_squared_M2)
    )
    
    PerformanceMetrics$INTEGRATED_ABSOLUTE_ERRORS_BYHAND <- list(MARGIN1 = mean(val_b_cop_abs),
                                                                 MARGIN1 = mean(val_b_cop_abs_SOLO),
                                                                 MARGIN1_COX = mean(val_b_cox_abs),
                                                                 MARGIN2 =  mean(val_b_cop_abs_M2),
                                                                 MARGIN2_COX = mean(val_b_cox_abs_M2)
    )
    
    
    ################################################################################################################
    ################################################################################################################ Extract stuff
    
    ## 
    ListOfCoefficients <- list()
    ListOfRisks <- list()
    ListOfMSTOP <- list()
    ListOfDIAGNOSTICS <- list()
    
    
    
    ListOfCoefficients$Margin1 <- list(MU1 = coef(GLM_M1$mu1, which = ""),
                                       SIGMA1 = coef(GLM_M1$sigma1, which = ""))
    
    ListOfCoefficients$Margin2 <- list(MU2 = coef(GLM_M2$mu, which = ""),
                                       SIGMA2 = coef(GLM_M2$sigma, which = ""))
    
    ListOfCoefficients$DependenceParam <- coef(GLM_M1$rho, which = "")
    
    ######## SOLO:
    ListOfCoefficients$Margin1_SOLO <- list(MU1 = coef(GLM_M1_SOLO$mu, which = ""),
                                            SIGMA1 = coef(GLM_M1_SOLO$sigma, which = ""))
    
    ListOfCoefficients$DependenceParam_SOLO <- coef(GLM_COPPAR, which = "")
    
    
    
    ListOfCoefficients$COX_1 <- coef(GLM_COX_M1, which = "")
    ListOfCoefficients$COX_2 <- coef(GLM_COX_M2, which = "")
    
    
    ListOfRisks$Margin1 <- list(DuringTraining = risk_duringtraining_M1,
                                AllData = risk_final_M1) 
    
    ListOfRisks$Margin2 <- list(DuringTraining = risk_duringtraining_M2,
                                AllData = risk_final_M2) 
    
    ListOfRisks$DependenceParam <- list(DuringTraining = risk_duringtraining_DEPENDENCE,
                                        AllData = risk_final_DEPENDENCE) 
    
    ######## SOLO
    ListOfRisks$Margin1_SOLO <- list(DuringTraining = risk_duringtraining_M1_SOLO,
                                     AllData = risk_final_M1_SOLO) 
    
    ListOfRisks$DependenceParam_SOLO <- list(DuringTraining = risk_duringtraining_DEPENDENCE,
                                             AllData = risk_final_DEPENDENCE) 
    
    
    ListOfRisks$COX_1 <- list(DuringTraining = risk_duringtraining_M1_COX,
                              AllData = risk_final_M1_COX) 
    
    ListOfRisks$COX_2 <- list(DuringTraining = risk_duringtraining_M2_COX,
                              AllData = risk_final_M2_COX) 
    
    ListOfMSTOP$Margin1 <- mstop(GLM_M1)
    ListOfMSTOP$Margin2 <- mstop(GLM_M2)
    ListOfMSTOP$DependenceParam <- mstop(GLM_M1$rho)
    ListOfMSTOP$COX_1 <- mstop(GLM_COX_M1)
    ListOfMSTOP$COX_2 <- mstop(GLM_COX_M2)
    
    
    ################## SOLO
    ListOfMSTOP$Margin1_SOLO <- mstop(GLM_M1_SOLO)
    
    
    
    ListOfDIAGNOSTICS$Margin1 <- list(SurvRange = range(SURV1),
                                      PDFRange = range(PDF1))
    
    ListOfDIAGNOSTICS$Margin2 <- list(SurvRange = range(SURV2),
                                      PDFRange = range(PDF2))
    
    ListOfDIAGNOSTICS$Dependence <- list(RangeOfDependence = range( exp( predict(GLM_M1$rho, type = "link") ) + 1 ),
                                         RangeOfTau = VineCopula::BiCopPar2Tau(family = COPFAM, 
                                                                               par = range( exp(predict(GLM_M1$rho, type = "link") ) + 1 ) ) )
    
    
    ######### SOLO
    ListOfDIAGNOSTICS$Margin1_SOLO <- list(SurvRange = range(SURV1_SOLO),
                                           PDFRange = range(PDF1_SOLO))
    
    ListOfDIAGNOSTICS$Dependence_SOLO <- list(RangeOfDependence = range(exp( predict(GLM_COPPAR, type = "link") ) + 1 ),
                                              RangeOfTau = VineCopula::BiCopPar2Tau(family = COPFAM, 
                                                                                    par = range(exp(predict(GLM_COPPAR, type = "link") ) + 1 ) ) )
    
    
    
    OtherConfigurations <- list(CensoringRates = CensoringRates, 
                                NumberOfCovariates = p, 
                                n.train = n.train,
                                n.mstop = n.mstop,
                                n.test = n.test,
                                Range_TrueKendall = TRUE_KENDALL_RANGE,
                                TrueCoefficients = list(MU1 = beta11, 
                                                        SIGMA1 = beta12,
                                                        MU2 = beta21,
                                                        SIGMA2 = beta22,
                                                        DEPENDENCE = betarho)
    )
    
    
    output <- list(COEFFICIENTS = ListOfCoefficients, 
                   PEFORMANCEMETRICS = PerformanceMetrics,
                   LOSS_VALUES = ListOfLOSS,
                   RISKS = ListOfRisks,
                   MSTOP = ListOfMSTOP, 
                   DIAGNOSTICS = ListOfDIAGNOSTICS, 
                   CONFIGURATIONS = OtherConfigurations)
    
    
    
  }
  
  return(output)
  
}