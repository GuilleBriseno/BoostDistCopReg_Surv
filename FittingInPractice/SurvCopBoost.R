###


SurvCopBoost <- function(formulas, margins = c("WEIBULL", "WEIBULL"), copula = c("CLAYTON"), 
                         response_1, response_2, data,
                         oobag_weights, mstops, s_step, stabilization, trace = TRUE){
  
  output <- list()
  
  output$ModelName <- paste("Bivariate copula model for right-censored data:", copula, "copula with", margins[1], "and", margins[2], "margins.")
  
  fitting_data <- data
  
  Margin1_Formula <- formulas[[1]]
  Margin2_Formula <- formulas[[2]]
  Dependence_Formula <- formulas[[3]]
  
  
  fitting_data_M1 <- as.data.frame(cbind(response_1, fitting_data))
  fitting_data_M2 <- as.data.frame(cbind(response_2, fitting_data))
  
  #### Assign families: 
  if(marings[1] == "WEIBULL"){
    
    Margin1Family <- Custom_WeibullFamily
    
  }
  if(marings[1] == "LOGNORMAL"){
    
    Margin1Family <- Custom_LogNormalFamily
    
  }
  if(marings[1] == "LOGLOGISTIC"){
    
    Margin1Family <- Custom_LogLogisticFamily
    
  }
   
  
  if(marings[2] == "WEIBULL"){
    
    Margin2Family <- Custom_WeibullFamily
    
  }
  if(marings[2] == "LOGNORMAL"){
    
    Margin2Family <- Custom_LogNormalFamily
    
  }
  if(marings[2] == "LOGLOGISTIC"){
    
    Margin2Family <- Custom_LogLogisticFamily
    
  }
  
  
  if(copula == "CLAYTON"){
    
    DependenceFamily <- BivAFT_RC_ClaytonCopula_RhoSoloFamily
    
  }
  if(copula == "CLAYTON90"){
    
    DependenceFamily <- BivAFT_RC_ClaytonCopula_90_RhoSoloFamily
    
  }
  if(copula == "CLAYTON180"){
    
    DependenceFamily <- BivAFT_RC_ClaytonCopula_180_RhoSoloFamily
    
  }
  if(copula == "CLAYTON270"){
    
    DependenceFamily <- BivAFT_RC_ClaytonCopula_270_RhoSoloFamily
    
  }
  if(copula == "GUMBEL"){
    
    DependenceFamily <- BivAFT_RC_GumbelCopula_RhoSoloFamily
    
  }
  if(copula == "GUMBEL90"){
    
    DependenceFamily <- BivAFT_RC_GumbelCopula_90_RhoSoloFamily
    
  }
  if(copula == "GUMBEL180"){
    
    DependenceFamily <- BivAFT_RC_GumbelCopula_180_RhoSoloFamily
    
  }
  if(copula == "GUMBEL270"){
    
    DependenceFamily <- BivAFT_RC_GumbelCopula_270_RhoSoloFamily
    
  }
  if(copula == "JOE"){
    
    DependenceFamily <- BivAFT_RC_JoeCopula_RhoSoloFamily
    
  }
  if(copula == "JOE90"){
    
    DependenceFamily <- BivAFT_RC_JoeCopula_90_RhoSoloFamily
    
  }
  if(copula == "JOE180"){
    
    DependenceFamily <- BivAFT_RC_JoeCopula_180_RhoSoloFamily
    
  }
  if(copula == "JOE270"){
    
    DependenceFamily <- BivAFT_RC_JoeCopula_270_RhoSoloFamily
    
  }
  if(copula == "NORMAL"){
    
    DependenceFamily <- BivAFT_RC_GaussCopula_RhoSoloFamily
    
  }
  if(copula == "FRANK"){
    
    DependenceFamily <- BivAFT_RC_FrankCopula_RhoSoloFamily
    
  }
  
  
  
  ####################################################
  cat("Boosting margins...")
  GLM1TRY <- myTryCatch( gamboostLSS(formula = Margin1_Formula, 
                                     data = fitting_data_M1,
                                     families = Margin1Family(stabilization = stabilization), 
                                     weights = oobag_weights, 
                                     method = "noncyclic",
                                     control = boost_control(mstop = mstops[1], 
                                                             nu = s_step, 
                                                             risk = "oobag", 
                                                             trace = trace)) )
  
  if(is.null(GLM1TRY$value)){
    
    
    
    output <- list(ERROR = GLM1TRY$the_error, 
                   WARNING = GLM1TRY$the_warning)
    
    
    warning("Could not fit model of margin 1. \nTry: changing the step length (s_step) or the number of fitting iterations (mstop).")
    return(output)
    
  }else{
    
    
    GLM_M1 <- GLM1TRY$value
    
    MSTOP_OPT_M1 <- which.min(risk(GLM_M1, merge = TRUE))
    
    Margin1Model <- gamboostLSS(formula = Margin1_Formula, 
                          data = fitting_data_M1, 
                          families = Margin1Family(stabilization = stabilization), 
                          method = "noncyclic",
                          control = boost_control(mstop = MSTOP_OPT_M1, 
                                                  nu = s_step, 
                                                  trace = trace))
  }
   
  
  
  Margin2Model <- myTryCatch( gamboostLSS(formula = Margin2_Formula, 
                                     data = fitting_data_M2,
                                     families = Margin2Family(stabilization = stabilization), 
                                     weights = oobag_weights, 
                                     method = "noncyclic",
                                     control = boost_control(mstop = mstops[2], 
                                                             nu = s_step, 
                                                             risk = "oobag", 
                                                             trace = trace)) )
  
  if(is.null(GLM2TRY$value)){
    
    
    
    output <- list(ERROR = GLM2TRY$the_error, 
                   WARNING = GLM2TRY$the_warning)
    
    
    warning("Could not fit model of margin 2. \nTry: changing the step length (s_step) or the number of fitting iterations (mstop).")
    return(output)
    
  }else{
    
    
    GLM_M2 <- GLM2TRY$value
    
    MSTOP_OPT_M2 <- which.min(risk(GLM_M2, merge = TRUE))
    
    GLM_M2 <- gamboostLSS(formula = Margin2_Formula, 
                          data = fitting_data_M2, 
                          families = Margin2Family(stabilization = stabilization), 
                          method = "noncyclic",
                          control = boost_control(mstop = MSTOP_OPT_M2, 
                                                  nu = s_step, 
                                                  trace = trace))
  }
  
  
  
  
  ###### COMPUTE THE FUNCTIONS: 
  
  ##### Predict some quantities required for the copula parameter: 
  mu1_hat <- predict(Margin1Model$mu, type = "response")
  sigma1_hat <- predict(Margin1Model$sigma, type = "response")
   
  mu2_hat <- predict(Margin2Model$mu, type = "response")
  sigma2_hat <- predict(Margin2Model$sigma, type = "response")
  
  
  if(marings[1] == "WEIBULL"){
    
    SURV1 <- pweibull(q = response1[,1], scale = mu1_hat, shape = sigma1_hat, lower.tail = FALSE) 
    PDF1 <- dweibull(x = response1[,1], scale = mu1_hat, shape = sigma1_hat)
    
  }
  if(marings[1] == "LOGNORMAL"){
    
    SURV1 <- plnorm(q = response_1[,1], meanlog = mu1_hat, sdlog = sigma1_hat)
    PDF1 <- dlnorm(x = response_1[,1], mealong = mu1_hat, sdlog = sigma1_hat)
    
  }
  if(marings[1] == "LOGLOGISTIC"){
    
    SURV1 <- 1 - pLogLogistic(x = response_1[,1], mu = mu1_hat, sigma = sigma1_hat)
    PDF1 <- dLogLogistic(x = response_1[,1], mu = mu1_hat, sigma = sigma1_hat)
  }
  
  
  if(marings[2] == "WEIBULL"){
    
    SURV2 <- pweibull(q = response2[,1], scale = mu2_hat, shape = sigma2_hat, lower.tail = FALSE) 
    PDF2 <- dweibull(x = response2[,1], scale = mu2_hat, shape = sigma2_hat)
    
  }
  if(marings[2] == "LOGNORMAL"){
    
    SURV2 <- plnorm(q = response_2[,1], meanlog = mu2_hat, sdlog = sigma2_hat)
    PDF2 <- dlnorm(x = response_2[,1], mealong = mu2_hat, sdlog = sigma2_hat)
    
  }
  if(marings[2] == "LOGLOGISTIC"){
    
    SURV2 <- 1 - pLogLogistic(x = response_2[,1], mu = mu2_hat, sigma = sigma2_hat)
    PDF2 <- dLogLogistic(x = response_2[,1], mu = mu2_hat, sigma = sigma2_hat)
  }
  
  
  

  dependence_data <- data.frame(SURV1 = SURV1, PDF1 = PDF1, delta1 = response_1[,2], 
                           SURV2 = SURV2, PDF2 = PDF2, delta2 = response_2[,2])
  
  dependence_data <- cbind(dependence_data, fitting_data)
  
  
  cat("\nBoosting copula...")
  GLM_COPPAR <- gamboost(Dependence_Formula,  
                         dependence_data,
                         family = DependenceFamily(), 
                         weights = oobag_weights,
                         control = boost_control(mstop = mstops[3],
                                                 nu = s_step,
                                                 risk = "oobag", 
                                                 trace = trace), 
                         center = FALSE)
  
  risk_duringtraining_DEPENDENCE <- risk(GLM_COPPAR)
  
  MSTOP_OPT_COPPAR <- which.min(risk(GLM_COPPAR))
  
  
  rm(GLM_COPPAR)
  
  CopulaModel <- gamboost(Dependence_Formula, 
                         dependence_data, 
                         family = DependenceFamily(),
                         control = boost_control(mstop = MSTOP_OPT_COPPAR, 
                                                 nu = s_step,
                                                 trace = trace), 
                         center = FALSE)
  
 
  
  
  
  ### Attach models: 
  output$Model_Margin1 <- Margin1Model
  output$Model_Margin2 <- Margin2Model
  output$Model_Copula <- CopulaModel 
  
  cat("\nDone")
  
  
  return(output)
  
  
}


