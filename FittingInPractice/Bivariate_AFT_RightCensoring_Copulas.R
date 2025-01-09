
############### bivariate copula functions for right-censored data:  

# BIG credit to GJRM for the derivatives of the copula terms (copula density, h-functions and copula CDF) wrt the
# dependence parameter

# These functions are the negative log-likelihood (loss function) for bivariate continuous time-to-event data 
# with independent right-censoring. These functions depend on ONE parameter only, the copula dependence parameter. 

############################################################################################# Gaussian & Frank copulas

BivAFT_RC_GaussCopula_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- tanh(f) 
    
    
    ###### Copula expressions: 
    copdensity <-  1/sqrt(1 - thet^2)*exp(  - (thet^2*( qnorm(S1)^2 +  qnorm(S2)^2 ) - 2*thet*qnorm(S1)*qnorm(S2) ) / (2*(1 - thet^2)) )
    
    hfunc_m1 <- pnorm( (qnorm(S2) - thet*qnorm(S1))/sqrt(1 - thet^2)   )
    
    hfunc_m2 <- pnorm( (qnorm(S1) - thet*qnorm(S2))/sqrt(1 - thet^2)   )
    
    copcdf <-  VGAM:::pbinorm( qnorm(S1), qnorm(S2), cov12 = thet)
    
    ################################## Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits:
    
    dercopdensity_dtheta <- exp(-(thet * (thet * (qnorm(S1)^2 + qnorm(S2)^2) - 2 * (qnorm(S1) * qnorm(S2)))/(2 * (1 - thet^2)))) * (qnorm(S1) * qnorm(S2) + thet * 
                                                                                                                                      (1 - (qnorm(S1)^2 + qnorm(S2)^2 + thet * (thet * (qnorm(S1)^2 + qnorm(S2)^2) - 2 * (qnorm(S1) * qnorm(S2)))/(1 - thet^2))))/((1 - thet^2) * sqrt(1 - thet^2))
    
    
    derhfuncm1_dtheta <- (-(dnorm((qnorm(S2) - tanh(f) * qnorm(S1))/sqrt(1 - tanh(f)^2)) *  (1/cosh(f)^2 * qnorm(S1)/sqrt(1 - tanh(f)^2) - (qnorm(S2) - 
                                                                                                                                              tanh(f) * qnorm(S1)) * (0.5 * (2 * (1/cosh(f)^2 * tanh(f)) * (1 - tanh(f)^2)^-0.5))/sqrt(1 - tanh(f)^2)^2)))/(1/cosh(f)^2)
    
    
    
    
    derhfuncm2_dtheta <-  (-(dnorm((qnorm(S1) - tanh(f) * qnorm(S2))/sqrt(1 - tanh(f)^2)) *  (1/cosh(f)^2 * qnorm(S2)/sqrt(1 - tanh(f)^2) - (qnorm(S1) - tanh(f) * qnorm(S2)) * (0.5 * (2 * (1/cosh(f)^2 * 
                                                                                                                                                                                               tanh(f)) * (1 - tanh(f)^2)^-0.5))/sqrt(1 - tanh(f)^2)^2)))/(1/cosh(f)^2)
    
    
    dercopcdf_dtheta <- VGAM:::dbinorm(qnorm(S1),qnorm(S2), cov12=thet)
    
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * ( 1/cosh(f)^2 ) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * ( 1/cosh(f)^2 ) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * ( 1/cosh(f)^2 ) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * ( 1/cosh(f)^2 ) )
    
    
    return(ngr)
    
  },
  
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- tanh(f) 
    
    
    ###### Copula expressions: 
    copdensity <-  1/sqrt(1 - thet^2)*exp(  - (thet^2*( qnorm(S1)^2 +  qnorm(S2)^2 ) - 2*thet*qnorm(S1)*qnorm(S2) ) / (2*(1 - thet^2)) )
    
    hfunc_m1 <- pnorm( (qnorm(S2) - thet*qnorm(S1))/sqrt(1 - thet^2)   )
    
    hfunc_m2 <- pnorm( (qnorm(S1) - thet*qnorm(S2))/sqrt(1 - thet^2)   )
    
    copcdf <-  VGAM:::pbinorm( qnorm(S1), qnorm(S2), cov12 = thet)
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ tanh(f) },
  
  name = "Gaussian copula for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_FrankCopula_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- f 
    
    
    ###### Copula expressions: 
    copdensity <- (exp((1 + S1 + S2)* thet)* (-1 + exp(thet))* thet)/(exp((S1 + S2)* thet) - exp(thet)* (-1 + exp(S1* thet) + exp(S2* thet)))^2
    
    hfunc_m1 <- (exp(thet)* (-1 + exp(S2* thet)))/(-exp((S1 + S2)* thet) + exp(thet)* (-1 + exp(S1* thet) + exp(S2* thet)))
    
    hfunc_m2 <- (exp(thet)* (-1 + exp(S1* thet)))/(-exp((S1 + S2)* thet) + exp(thet)* (-1 + exp(S2* thet) + exp(S1* thet)))
    
    bit <- -expm1(-thet) # 1 - exp(-par1)
    
    copcdf <-  -(1/thet)*log( (bit - (1 - exp(-thet*S1))*(1 - exp(-thet*S2)))/bit ) 
    
    ################################### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t2 <- exp(thet);
    t3 <- t2-1.0;
    t4 <- thet*S2;
    t5 <- thet*S1;
    t7 <- exp(t4+t5+thet);
    t10 <- exp(t4+t5);
    t12 <- exp(t4+thet);
    t14 <- exp(t5+thet);
    t15 <- t10-t12-t14+t2;
    t16 <- t15*t15;
    t17 <- 1/t16;
    t21 <- thet*t3;
    
    dercopdensity_dtheta <- t3*t7*t17+thet*t2*t7*t17+t21*(S2+S1+1.0)*t7*t17-2.0*t21*t7/t15/t16*((S2+S1)*t10-(S2+1.0)*t12-(S1+1.0)*t14+t2);
    
    derhfuncm1_dtheta <- (exp(thet + 
                                S1 *thet)* (exp(2* S2* thet)* (-1 + S1) + exp(thet)* S1 - 
                                              exp(S2* thet)* (-1 + S1 + exp(thet)* (S1 - S2) + S2)))/(exp((S1 + 
                                                                                                             S2)* thet) - exp(thet)* (-1 + exp(S1* thet) + exp(S2* thet)))^2
    
    derhfuncm2_dtheta <- (exp(thet + 
                                S2 *thet)* (exp(2* S1* thet)* (-1 + S2) + exp(thet)* S2 - 
                                              exp(S1* thet)* (-1 + S2 + exp(thet)* (S2 - S1) + S1)))/(exp((S2 + 
                                                                                                             S1)* thet) - exp(thet)* (-1 + exp(S2* thet) + exp(S1* thet)))^2
    
    dercopcdf_dtheta <- (exp(thet)* (1/(-1 + exp(thet)) + (-1 - exp(S2* thet)* (-1 + S1) + S1 - exp(S1* thet)* (-1 + S2) + S2)/(exp((S1 + S2)* thet) - exp(thet)* (-1 + exp(S1* thet) + exp(S2* thet))))* thet + log((exp(-(S1 + S2)* thet)* (-exp((S1 + S2)* thet) + exp(thet)* (-1 + exp(S1* thet) + exp(S2* thet))))/(-1 + exp(thet))))/thet^2
    
    
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * 1  ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * 1  ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * 1  )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * 1 )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- f 
    
    
    ###### Copula expressions: 
    copdensity <- (exp((1 + S1 + S2)* thet)* (-1 + exp(thet))* thet)/(exp((S1 + S2)* thet) - exp(thet)* (-1 + exp(S1* thet) + exp(S2* thet)))^2
    
    hfunc_m1 <- (exp(thet)* (-1 + exp(S2* thet)))/(-exp((S1 + S2)* thet) + exp(thet)* (-1 + exp(S1* thet) + exp(S2* thet)))
    
    hfunc_m2 <- (exp(thet)* (-1 + exp(S1* thet)))/(-exp((S1 + S2)* thet) + exp(thet)* (-1 + exp(S2* thet) + exp(S1* thet)))
    
    bit <- -expm1(-thet) # 1 - exp(-par1)
    
    copcdf <-  -(1/thet)*log( (bit - (1 - exp(-thet*S1))*(1 - exp(-thet*S2)))/bit ) 
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ f },
  
  name = "Frank copula for bivariate time-to-event data with right-censoring scheme")
}

############################################################################################# Clayton & rotations
BivAFT_RC_ClaytonCopula_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f)
    
    
    ###### Copula expressions: 
    copdensity <- S1^(-1 - thet)* S2^(-1 - thet)* (-1 + S1^-thet + S2^-thet)^(-2 - 1/thet) *(1 + thet)
    
    hfunc_m1 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S1^((-thet) - 1) * (-thet)))
    
    hfunc_m2 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S2^((-thet) - 1) * (-thet)))
    
    copcdf <- ( S1^-thet + S2^-thet - 1 )^( -1/thet ) 
    
    
    #################################### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits:
    
    t1 = S1*S2;
    t2 = -thet-1.0;
    t3 = t1^(1.0*t2);
    t4 = S1^(-1.0*thet);
    t5 = S2^(-1.0*thet);
    t6 = t4+t5-1.0;
    t7 = -2.0-1/thet;
    t8 = t6^(1.0*t7);
    t9 = -t2*t3;
    t10 = log(t1);
    t11 = thet^2;
    t12 = log(t6);
    t13 = log(S1);
    t14 = log(S2);
    
    dercopdensity_dtheta <-  t3*t8-t9*t10*t8+t9*t8*(1/t11*t12+t7*(-t4*t13-t5*t14)/t6);
    
    derhfuncm1_dtheta <-  ((1 + 1/thet) * (log(S1)/S1^thet + log(S2)/S2^thet)/(1/S1^thet + 1/S2^thet - 1)^(1/thet + 2) + log(1/S1^thet + 1/S2^thet - 1)/(thet^2 * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)))/S1^(1 + thet) + (1/S1^(1 + thet) - (1/S1^(1 + thet) + thet * log(S1)/S1^(1 + thet)))/(thet * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)) 
    
    
    
    
    derhfuncm2_dtheta <- ((1 + 1/thet) * (log(S1)/S1^thet + log(S2)/S2^thet)/(1/S1^thet + 1/S2^thet - 1)^(1/thet + 2) + log(1/S1^thet + 1/S2^thet - 1)/(thet^2 * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)))/S2^(1 + thet) + (1/S2^(1 + thet) - (1/S2^(1 + thet) + thet * log(S2)/S2^(1 + thet)))/(thet * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet))
    
    
    
    
    dercopcdf_dtheta <-  ((-1 + S1^-thet + S2^-thet)^(-1/thet) *((thet *(S2^thet *log(S1) + S1^thet* log(S2)))/(S2^thet - S1^thet* (-1 + S2^thet)) + log(-1 + S1^-thet + S2^-thet)))/thet^2
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f)
    
    ###### Copula expressions: 
    copdensity <- S1^(-1 - thet)* S2^(-1 - thet)* (-1 + S1^-thet + S2^-thet)^(-2 - 1/thet) *(1 + thet)
    
    hfunc_m1 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S1^((-thet) - 1) * (-thet)))
    
    hfunc_m2 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S2^((-thet) - 1) * (-thet)))
    
    copcdf <- ( S1^-thet + S2^-thet - 1 )^( -1/thet ) 
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) },
  
  name = "Clayton copula for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_ClaytonCopula_180_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz( y[,1] )
    
    S2    <- pdffz( y[,4] )
    
    ####### ROTATION TERMS (used outside of the copula function)
    S1_Tilde <- pdffz(S1)
    
    S2_Tilde <- pdffz(S2)
    
    S1 <- pdffz(1 - S1)
    S2 <- pdffz(1 - S2)
    
    thet <- exp(f)
    
    
    ###### Copula expressions: 
    copdensity <- S1^(-1 - thet)* S2^(-1 - thet)* (-1 + S1^-thet + S2^-thet)^(-2 - 1/thet) *(1 + thet)
    
    hfunc_m1 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S1^((-thet) - 1) * (-thet)))
    
    hfunc_m2 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S2^((-thet) - 1) * (-thet)))
    
    copcdf <- ( S1^-thet + S2^-thet - 1 )^( -1/thet ) 
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz(copdensity)
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S1_Tilde + S2_Tilde - 1 + copcdf)
    
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t1 = S1*S2;
    t2 = -thet-1.0;
    t3 = t1^(1.0*t2);
    t4 = S1^(-1.0*thet);
    t5 = S2^(-1.0*thet);
    t6 = t4+t5-1.0;
    t7 = -2.0-1/thet;
    t8 = t6^(1.0*t7);
    t9 = -t2*t3;
    t10 = log(t1);
    t11 = thet^2;
    t12 = log(t6);
    t13 = log(S1);
    t14 = log(S2);
    
    dercopdensity_dtheta <-  t3*t8-t9*t10*t8+t9*t8*(1/t11*t12+t7*(-t4*t13-t5*t14)/t6);
    
    derhfuncm1_dtheta <-  ((1 + 1/thet) * (log(S1)/S1^thet + log(S2)/S2^thet)/(1/S1^thet + 1/S2^thet - 1)^(1/thet + 2) + log(1/S1^thet + 1/S2^thet - 1)/(thet^2 * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)))/S1^(1 + thet) + (1/S1^(1 + thet) - (1/S1^(1 + thet) + thet * log(S1)/S1^(1 + thet)))/(thet * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)) 
    
    
    
    
    derhfuncm2_dtheta <- ((1 + 1/thet) * (log(S1)/S1^thet + log(S2)/S2^thet)/(1/S1^thet + 1/S2^thet - 1)^(1/thet + 2) + log(1/S1^thet + 1/S2^thet - 1)/(thet^2 * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)))/S2^(1 + thet) + (1/S2^(1 + thet) - (1/S2^(1 + thet) + thet * log(S2)/S2^(1 + thet)))/(thet * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet))
    
    
    
    
    dercopcdf_dtheta <-  ((-1 + S1^-thet + S2^-thet)^(-1/thet) *((thet *(S2^thet *log(S1) + S1^thet* log(S2)))/(S2^thet - S1^thet* (-1 + S2^thet)) + log(-1 + S1^-thet + S2^-thet)))/thet^2
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <- - derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <- - derhfuncm2_dtheta
    
    dercopcdf_dtheta        <- dercopcdf_dtheta
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    ####### ROTATION TERMS (used outside of the copula function)
    S1_Tilde <- pdffz(S1)
    S2_Tilde <- pdffz(S2)
    
    S1 <- pdffz(1 - S1)
    S2 <- pdffz(1 - S2)
    
    thet <- exp(f)
    
    ###### Copula expressions: 
    copdensity <- S1^(-1 - thet)* S2^(-1 - thet)* (-1 + S1^-thet + S2^-thet)^(-2 - 1/thet) *(1 + thet)
    
    hfunc_m1 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S1^((-thet) - 1) * (-thet)))
    
    hfunc_m2 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S2^((-thet) - 1) * (-thet)))
    
    copcdf <- ( S1^-thet + S2^-thet - 1 )^( -1/thet ) 
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz(copdensity)
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S1_Tilde + S2_Tilde - 1 + copcdf)
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) },
  
  name = "Clayton copula (180° rotation) for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_ClaytonCopula_90_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz( y[,1] )
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f)
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S1 <- pdffz(1 - S1)
    
    
    
    ###### Copula expressions: 
    copdensity <- S1^(-1 - thet)* S2^(-1 - thet)* (-1 + S1^-thet + S2^-thet)^(-2 - 1/thet) *(1 + thet)
    
    hfunc_m1 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S1^((-thet) - 1) * (-thet)))
    
    hfunc_m2 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S2^((-thet) - 1) * (-thet)))
    
    copcdf <- ( S1^-thet + S2^-thet - 1 )^( -1/thet ) 
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S2 - copcdf)
    
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t1 = S1*S2;
    t2 = -thet-1.0;
    t3 = t1^(1.0*t2);
    t4 = S1^(-1.0*thet);
    t5 = S2^(-1.0*thet);
    t6 = t4+t5-1.0;
    t7 = -2.0-1/thet;
    t8 = t6^(1.0*t7);
    t9 = -t2*t3;
    t10 = log(t1);
    t11 = thet^2;
    t12 = log(t6);
    t13 = log(S1);
    t14 = log(S2);
    
    dercopdensity_dtheta <-  t3*t8-t9*t10*t8+t9*t8*(1/t11*t12+t7*(-t4*t13-t5*t14)/t6);
    
    derhfuncm1_dtheta <-  ((1 + 1/thet) * (log(S1)/S1^thet + log(S2)/S2^thet)/(1/S1^thet + 1/S2^thet - 1)^(1/thet + 2) + log(1/S1^thet + 1/S2^thet - 1)/(thet^2 * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)))/S1^(1 + thet) + (1/S1^(1 + thet) - (1/S1^(1 + thet) + thet * log(S1)/S1^(1 + thet)))/(thet * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)) 
    
    
    
    
    derhfuncm2_dtheta <- ((1 + 1/thet) * (log(S1)/S1^thet + log(S2)/S2^thet)/(1/S1^thet + 1/S2^thet - 1)^(1/thet + 2) + log(1/S1^thet + 1/S2^thet - 1)/(thet^2 * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)))/S2^(1 + thet) + (1/S2^(1 + thet) - (1/S2^(1 + thet) + thet * log(S2)/S2^(1 + thet)))/(thet * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet))
    
    
    
    
    dercopcdf_dtheta <-  ((-1 + S1^-thet + S2^-thet)^(-1/thet) *((thet *(S2^thet *log(S1) + S1^thet* log(S2)))/(S2^thet - S1^thet* (-1 + S2^thet)) + log(-1 + S1^-thet + S2^-thet)))/thet^2
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <-   derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <- - derhfuncm2_dtheta
    
    dercopcdf_dtheta        <-  - dercopcdf_dtheta
    
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f)
    
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S1 <- pdffz(1 - S1)
    
    
    
    ###### Copula expressions: 
    copdensity <- S1^(-1 - thet)* S2^(-1 - thet)* (-1 + S1^-thet + S2^-thet)^(-2 - 1/thet) *(1 + thet)
    
    hfunc_m1 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S1^((-thet) - 1) * (-thet)))
    
    hfunc_m2 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S2^((-thet) - 1) * (-thet)))
    
    copcdf <- ( S1^-thet + S2^-thet - 1 )^( -1/thet ) 
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S2 - copcdf)
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) },
  
  name = "Clayton copula (90° rotation) for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_ClaytonCopula_270_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz( y[,1] )
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f)
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S2 <- pdffz(1 - S2)
    
    
    
    ###### Copula expressions: 
    copdensity <- S1^(-1 - thet)* S2^(-1 - thet)* (-1 + S1^-thet + S2^-thet)^(-2 - 1/thet) *(1 + thet)
    
    hfunc_m1 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S1^((-thet) - 1) * (-thet)))
    
    hfunc_m2 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S2^((-thet) - 1) * (-thet)))
    
    copcdf <- ( S1^-thet + S2^-thet - 1 )^( -1/thet ) 
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(hfunc_m2)
    
    copcdf      <- pdffz(S1 - copcdf)
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t1 = S1*S2;
    t2 = -thet-1.0;
    t3 = t1^(1.0*t2);
    t4 = S1^(-1.0*thet);
    t5 = S2^(-1.0*thet);
    t6 = t4+t5-1.0;
    t7 = -2.0-1/thet;
    t8 = t6^(1.0*t7);
    t9 = -t2*t3;
    t10 = log(t1);
    t11 = thet^2;
    t12 = log(t6);
    t13 = log(S1);
    t14 = log(S2);
    
    dercopdensity_dtheta <-  t3*t8-t9*t10*t8+t9*t8*(1/t11*t12+t7*(-t4*t13-t5*t14)/t6);
    
    derhfuncm1_dtheta <-  ((1 + 1/thet) * (log(S1)/S1^thet + log(S2)/S2^thet)/(1/S1^thet + 1/S2^thet - 1)^(1/thet + 2) + log(1/S1^thet + 1/S2^thet - 1)/(thet^2 * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)))/S1^(1 + thet) + (1/S1^(1 + thet) - (1/S1^(1 + thet) + thet * log(S1)/S1^(1 + thet)))/(thet * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)) 
    
    
    
    
    derhfuncm2_dtheta <- ((1 + 1/thet) * (log(S1)/S1^thet + log(S2)/S2^thet)/(1/S1^thet + 1/S2^thet - 1)^(1/thet + 2) + log(1/S1^thet + 1/S2^thet - 1)/(thet^2 * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet)))/S2^(1 + thet) + (1/S2^(1 + thet) - (1/S2^(1 + thet) + thet * log(S2)/S2^(1 + thet)))/(thet * (1/S1^thet + 1/S2^thet - 1)^(1 + 1/thet))
    
    
    
    
    dercopcdf_dtheta <-  ((-1 + S1^-thet + S2^-thet)^(-1/thet) *((thet *(S2^thet *log(S1) + S1^thet* log(S2)))/(S2^thet - S1^thet* (-1 + S2^thet)) + log(-1 + S1^-thet + S2^-thet)))/thet^2
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <-  - derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <-    derhfuncm2_dtheta
    
    dercopcdf_dtheta        <-   - dercopcdf_dtheta
    
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f)
    
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - v
    S2 <- pdffz(1 - S2)
    
    
    
    ###### Copula expressions: 
    copdensity <- S1^(-1 - thet)* S2^(-1 - thet)* (-1 + S1^-thet + S2^-thet)^(-2 - 1/thet) *(1 + thet)
    
    hfunc_m1 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S1^((-thet) - 1) * (-thet)))
    
    hfunc_m2 <- (S1^(-thet) + S2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (S2^((-thet) - 1) * (-thet)))
    
    copcdf <- ( S1^-thet + S2^-thet - 1 )^( -1/thet ) 
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(hfunc_m2)
    
    copcdf      <- pdffz(S1 - copcdf)
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) },
  
  name = "Clayton copula (270° rotation) for bivariate time-to-event data with right-censoring scheme")
}


############################################################################################# Gumbel & rotations
BivAFT_RC_GumbelCopula_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f) + 1
    
    
    ###### Copula expressions: 
    copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
    
    hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
    
    hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
    
    copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
    ################################
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t3 = log(S1);
    t4 = (-t3)^thet
    t5 = log(S2);
    t6 = (-t5)^thet
    t7 = t4+t6;
    t8 = 1/thet;
    t9 = t7^t8
    t10 = thet^2;
    t12 = log(t7);
    t13 = 1/t10*t12;
    t14 = log(-t3);
    t16 = log(-t5);
    t18 = t4*t14+t6*t16;
    t20 = 1/t7;
    t22 = -t13+t8*t18*t20;
    t24 = exp(-t9);
    t26 = t24/S1;
    t28 = 1/S2;
    t29 = -1.0+t8;
    t30 = t7^(2*t29)
    t32 = t3*t5;
    t33 = thet-1.0;
    t34 = t32^t33
    t35 = t7^(-1.0*t8)
    t36 = t33*t35;
    t17 = 1.0+t36;
    t15 = t34*t17;
    t11 = t26*t28;
    t2 = t30*t34;
    t1 = log(t32);
    
    
    dercopdensity_dtheta <- -t9*t22*t26*t28*t30*t15+t11*t30*(-2.0*t13+2.0*t29*t18*t20)*t15+t11*t2*t1*t17+t11*t2*(t35-t36*t22);
    
    
    derhfuncm1_dtheta <- ((-log(S1))^(-1 + thet) * (((-log(S1))^thet * log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 2) * (1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 1) * log((-log(S1))^thet + (-log(S2))^thet)/thet^2) + ((-log(S1))^(-1 + thet) * log(-log(S1)) - (-log(S1))^(-1 + thet) * (((-log(S1))^thet * 
                                                                                                                                                                                                                                                                                                                                                                       log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                             (-log(S2))^thet)^(1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet) * 
                                                                                                                                                                                                                                                                                                                                                                      log((-log(S1))^thet + (-log(S2))^thet)/thet)/thet) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                              (-log(S2))^thet)^(1/thet - 1)) * exp(-((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (-log(S2))^thet)^(1/thet))/S1
    
    
    
    
    derhfuncm2_dtheta <- ((-log(S2))^(-1 + thet) * (((-log(S1))^thet * log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 2) * (1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 1) * log((-log(S1))^thet + (-log(S2))^thet)/thet^2) + ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                   (-log(S2))^thet)^(1/thet - 1) * ((-log(S2))^(-1 + thet) * 
                                                                                                                                                                                                                                                                                                                                      log(-log(S2)) - (-log(S2))^(-1 + thet) * (((-log(S1))^thet * 
                                                                                                                                                                                                                                                                                                                                                                                   log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                         (-log(S2))^thet)^(1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet) * 
                                                                                                                                                                                                                                                                                                                                                                                  log((-log(S1))^thet + (-log(S2))^thet)/thet)/thet)) * exp(-((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                (-log(S2))^thet)^(1/thet))/S2
    
    
    
    
    dercopcdf_dtheta <-  (1/(thet^2))*exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet)* (-thet* (-log(S1))^thet* log(-log(S1)) + ((-log(S1))^thet + (-log(S2))^thet)* log((-log(S1))^thet + (-log(S2))^thet) - thet *(-log(S2))^thet* log(-log(S2)))
    
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f) + 1
    
    ###### Copula expressions: 
    copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
    
    hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
    
    hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
    
    copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) + 1 },
  
  name = "Gumbel copula for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_GumbelCopula_180_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    S1_Tilde <- pdffz(S1)
    
    S2_Tilde <- pdffz(S2)
    
    S1 <- pdffz(1 - S1)
    S2 <- pdffz(1 - S2)
    
    
    ###### Copula expressions: 
    copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
    
    hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
    
    hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
    
    copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz(copdensity)
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S1_Tilde + S2_Tilde - 1 + copcdf)
    
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t3 = log(S1);
    t4 = (-t3)^thet
    t5 = log(S2);
    t6 = (-t5)^thet
    t7 = t4+t6;
    t8 = 1/thet;
    t9 = t7^t8
    t10 = thet^2;
    t12 = log(t7);
    t13 = 1/t10*t12;
    t14 = log(-t3);
    t16 = log(-t5);
    t18 = t4*t14+t6*t16;
    t20 = 1/t7;
    t22 = -t13+t8*t18*t20;
    t24 = exp(-t9);
    t26 = t24/S1;
    t28 = 1/S2;
    t29 = -1.0+t8;
    t30 = t7^(2*t29)
    t32 = t3*t5;
    t33 = thet-1.0;
    t34 = t32^t33
    t35 = t7^(-1.0*t8)
    t36 = t33*t35;
    t17 = 1.0+t36;
    t15 = t34*t17;
    t11 = t26*t28;
    t2 = t30*t34;
    t1 = log(t32);
    
    
    dercopdensity_dtheta <- -t9*t22*t26*t28*t30*t15+t11*t30*(-2.0*t13+2.0*t29*t18*t20)*t15+t11*t2*t1*t17+t11*t2*(t35-t36*t22);
    
    
    derhfuncm1_dtheta <- ((-log(S1))^(-1 + thet) * (((-log(S1))^thet * log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 2) * (1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 1) * log((-log(S1))^thet + (-log(S2))^thet)/thet^2) + ((-log(S1))^(-1 + thet) * log(-log(S1)) - (-log(S1))^(-1 + thet) * (((-log(S1))^thet * 
                                                                                                                                                                                                                                                                                                                                                                       log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                             (-log(S2))^thet)^(1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet) * 
                                                                                                                                                                                                                                                                                                                                                                      log((-log(S1))^thet + (-log(S2))^thet)/thet)/thet) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                              (-log(S2))^thet)^(1/thet - 1)) * exp(-((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (-log(S2))^thet)^(1/thet))/S1
    
    
    
    
    derhfuncm2_dtheta <- ((-log(S2))^(-1 + thet) * (((-log(S1))^thet * log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 2) * (1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 1) * log((-log(S1))^thet + (-log(S2))^thet)/thet^2) + ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                   (-log(S2))^thet)^(1/thet - 1) * ((-log(S2))^(-1 + thet) * 
                                                                                                                                                                                                                                                                                                                                      log(-log(S2)) - (-log(S2))^(-1 + thet) * (((-log(S1))^thet * 
                                                                                                                                                                                                                                                                                                                                                                                   log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                         (-log(S2))^thet)^(1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet) * 
                                                                                                                                                                                                                                                                                                                                                                                  log((-log(S1))^thet + (-log(S2))^thet)/thet)/thet)) * exp(-((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                (-log(S2))^thet)^(1/thet))/S2
    
    
    
    
    dercopcdf_dtheta <-  (1/(thet^2))*exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet)* (-thet* (-log(S1))^thet* log(-log(S1)) + ((-log(S1))^thet + (-log(S2))^thet)* log((-log(S1))^thet + (-log(S2))^thet) - thet *(-log(S2))^thet* log(-log(S2)))
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <- - derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <- - derhfuncm2_dtheta
    
    dercopcdf_dtheta        <- dercopcdf_dtheta
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    S1_Tilde <- pdffz(S1)
    S2_Tilde <- pdffz(S2)
    
    S1 <- pdffz(1 - S1)
    S2 <- pdffz(1 - S2)
    
    
    ###### Copula expressions: 
    copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
    
    hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
    
    hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
    
    copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz(copdensity)
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S1_Tilde + S2_Tilde - 1 + copcdf)
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) + 1 },
  
  name = "Gumbel copula (180° rotation) for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_GumbelCopula_90_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S1 <- pdffz(1 - S1)
    
    
    ###### Copula expressions: 
    copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
    
    hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
    
    hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
    
    copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S2 - copcdf)
    
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t3 = log(S1);
    t4 = (-t3)^thet
    t5 = log(S2);
    t6 = (-t5)^thet
    t7 = t4+t6;
    t8 = 1/thet;
    t9 = t7^t8
    t10 = thet^2;
    t12 = log(t7);
    t13 = 1/t10*t12;
    t14 = log(-t3);
    t16 = log(-t5);
    t18 = t4*t14+t6*t16;
    t20 = 1/t7;
    t22 = -t13+t8*t18*t20;
    t24 = exp(-t9);
    t26 = t24/S1;
    t28 = 1/S2;
    t29 = -1.0+t8;
    t30 = t7^(2*t29)
    t32 = t3*t5;
    t33 = thet-1.0;
    t34 = t32^t33
    t35 = t7^(-1.0*t8)
    t36 = t33*t35;
    t17 = 1.0+t36;
    t15 = t34*t17;
    t11 = t26*t28;
    t2 = t30*t34;
    t1 = log(t32);
    
    
    dercopdensity_dtheta <- -t9*t22*t26*t28*t30*t15+t11*t30*(-2.0*t13+2.0*t29*t18*t20)*t15+t11*t2*t1*t17+t11*t2*(t35-t36*t22);
    
    
    derhfuncm1_dtheta <- ((-log(S1))^(-1 + thet) * (((-log(S1))^thet * log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 2) * (1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 1) * log((-log(S1))^thet + (-log(S2))^thet)/thet^2) + ((-log(S1))^(-1 + thet) * log(-log(S1)) - (-log(S1))^(-1 + thet) * (((-log(S1))^thet * 
                                                                                                                                                                                                                                                                                                                                                                       log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                             (-log(S2))^thet)^(1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet) * 
                                                                                                                                                                                                                                                                                                                                                                      log((-log(S1))^thet + (-log(S2))^thet)/thet)/thet) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                              (-log(S2))^thet)^(1/thet - 1)) * exp(-((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (-log(S2))^thet)^(1/thet))/S1
    
    
    
    
    derhfuncm2_dtheta <- ((-log(S2))^(-1 + thet) * (((-log(S1))^thet * log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 2) * (1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 1) * log((-log(S1))^thet + (-log(S2))^thet)/thet^2) + ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                   (-log(S2))^thet)^(1/thet - 1) * ((-log(S2))^(-1 + thet) * 
                                                                                                                                                                                                                                                                                                                                      log(-log(S2)) - (-log(S2))^(-1 + thet) * (((-log(S1))^thet * 
                                                                                                                                                                                                                                                                                                                                                                                   log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                         (-log(S2))^thet)^(1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet) * 
                                                                                                                                                                                                                                                                                                                                                                                  log((-log(S1))^thet + (-log(S2))^thet)/thet)/thet)) * exp(-((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                (-log(S2))^thet)^(1/thet))/S2
    
    
    
    
    dercopcdf_dtheta <-  (1/(thet^2))*exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet)* (-thet* (-log(S1))^thet* log(-log(S1)) + ((-log(S1))^thet + (-log(S2))^thet)* log((-log(S1))^thet + (-log(S2))^thet) - thet *(-log(S2))^thet* log(-log(S2)))
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <-   derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <- - derhfuncm2_dtheta
    
    dercopcdf_dtheta        <-  - dercopcdf_dtheta
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S1 <- pdffz(1 - S1)
    
    
    ###### Copula expressions: 
    copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
    
    hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
    
    hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
    
    copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S2 - copcdf)
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) + 1 },
  
  name = "Gumbel copula (90° rotation) for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_GumbelCopula_270_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S2 <- pdffz(1 - S2)
    
    
    ###### Copula expressions: 
    copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
    
    hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
    
    hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
    
    copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(hfunc_m2)
    
    copcdf      <- pdffz(S1 - copcdf)
    
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t3 = log(S1);
    t4 = (-t3)^thet
    t5 = log(S2);
    t6 = (-t5)^thet
    t7 = t4+t6;
    t8 = 1/thet;
    t9 = t7^t8
    t10 = thet^2;
    t12 = log(t7);
    t13 = 1/t10*t12;
    t14 = log(-t3);
    t16 = log(-t5);
    t18 = t4*t14+t6*t16;
    t20 = 1/t7;
    t22 = -t13+t8*t18*t20;
    t24 = exp(-t9);
    t26 = t24/S1;
    t28 = 1/S2;
    t29 = -1.0+t8;
    t30 = t7^(2*t29)
    t32 = t3*t5;
    t33 = thet-1.0;
    t34 = t32^t33
    t35 = t7^(-1.0*t8)
    t36 = t33*t35;
    t17 = 1.0+t36;
    t15 = t34*t17;
    t11 = t26*t28;
    t2 = t30*t34;
    t1 = log(t32);
    
    
    dercopdensity_dtheta <- -t9*t22*t26*t28*t30*t15+t11*t30*(-2.0*t13+2.0*t29*t18*t20)*t15+t11*t2*t1*t17+t11*t2*(t35-t36*t22);
    
    
    derhfuncm1_dtheta <- ((-log(S1))^(-1 + thet) * (((-log(S1))^thet * log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 2) * (1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 1) * log((-log(S1))^thet + (-log(S2))^thet)/thet^2) + ((-log(S1))^(-1 + thet) * log(-log(S1)) - (-log(S1))^(-1 + thet) * (((-log(S1))^thet * 
                                                                                                                                                                                                                                                                                                                                                                       log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                             (-log(S2))^thet)^(1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet) * 
                                                                                                                                                                                                                                                                                                                                                                      log((-log(S1))^thet + (-log(S2))^thet)/thet)/thet) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                              (-log(S2))^thet)^(1/thet - 1)) * exp(-((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (-log(S2))^thet)^(1/thet))/S1
    
    
    
    
    derhfuncm2_dtheta <- ((-log(S2))^(-1 + thet) * (((-log(S1))^thet * log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 2) * (1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet - 1) * log((-log(S1))^thet + (-log(S2))^thet)/thet^2) + ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                   (-log(S2))^thet)^(1/thet - 1) * ((-log(S2))^(-1 + thet) * 
                                                                                                                                                                                                                                                                                                                                      log(-log(S2)) - (-log(S2))^(-1 + thet) * (((-log(S1))^thet * 
                                                                                                                                                                                                                                                                                                                                                                                   log(-log(S1)) + (-log(S2))^thet * log(-log(S2))) * ((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                         (-log(S2))^thet)^(1/thet - 1) - ((-log(S1))^thet + (-log(S2))^thet)^(1/thet) * 
                                                                                                                                                                                                                                                                                                                                                                                  log((-log(S1))^thet + (-log(S2))^thet)/thet)/thet)) * exp(-((-log(S1))^thet + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                (-log(S2))^thet)^(1/thet))/S2
    
    
    
    
    dercopcdf_dtheta <-  (1/(thet^2))*exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet)* (-thet* (-log(S1))^thet* log(-log(S1)) + ((-log(S1))^thet + (-log(S2))^thet)* log((-log(S1))^thet + (-log(S2))^thet) - thet *(-log(S2))^thet* log(-log(S2)))
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <-  - derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <-    derhfuncm2_dtheta
    
    dercopcdf_dtheta        <-   - dercopcdf_dtheta
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - v
    S2 <- pdffz(1 - S2)
    
    
    ###### Copula expressions: 
    copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
    
    hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
    
    hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
    
    copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(hfunc_m2)
    
    copcdf      <- pdffz(S1 - copcdf)
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) + 1 },
  
  name = "Gumbel copula (270° rotation) for bivariate time-to-event data with right-censoring scheme")
}


############################################################################################# Joe & rotations
BivAFT_RC_JoeCopula_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f) + 1
    
    
    ###### Copula expressions: 
    copdensity <- (1 - S1)^(-1 + thet)* ((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)^(-2 + 1/thet) *(1 - S2)^(-1 + thet)* (-(-1 + (1 - S1)^thet)* (-1 + (1 - S2)^thet) + thet)
    
    hfunc_m1 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S1)^(thet - 1) * thet - (1 - S1)^(thet - 1) * thet * (1 - S2)^thet))
    
    hfunc_m2 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S2)^(thet - 1) * thet - (1 - S1)^thet * ((1 - S2)^(thet - 1) * thet)))
    
    bit1 <- (1 - S1)^thet
    
    bit2 <- (1 - S2)^thet
    
    copcdf <-  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)
    ################################
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t1 = 1.0-S1;
    t2 = t1^(1.0*thet);
    t3 = 1.0-S2;
    t4 = t3^(1.0*thet);
    t5 = t2*t4;
    t6 = t2+t4-t5;
    t8 = 1/thet-2.0;
    t9 = t6^(1.0*t8);
    t10 = thet^2;
    t11 = log(t6);
    t12 = log(t1);
    t13 = t2*t12;
    t14 = log(t3);
    t15 = t4*t14;
    t16 = t13*t4;
    t19 = t5*t14;
    t21 = thet-1.0;
    t27 = t1^(1.0*t21);
    t28 = t3^(1.0*t21);
    t30 = thet-1.0+t2+t4-t5;
    t33 = t9*t27;
    
    dercopdensity_dtheta <- t9*(-1/t10*t11+t8*(t13+t15-t16-t19)/t6)*t27*t28*t30+t33*t12*t28*t30+t33*t28*t14*t30+t33*t28*(1.0+t13+t15-t16-t19);
    
    
    derhfuncm1_dtheta <- ((((1 - S1)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - S1) + 
                             ((1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - 
                                                                                     S2)) * ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                               (1 - S2)^thet)^(1/thet - 2) * (1/thet - 1) - ((1 - S1)^thet + 
                                                                                                                                               (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                 1) * log((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                                                                                                                                            (1 - S2)^thet)/thet^2) * ((1 - S1)^(thet - 1) - (1 - S1)^(thet - 
                                                                                                                                                                                                                                                                        1) * (1 - S2)^thet) + ((1 - S1)^(thet - 1) + (1 - S1)^(thet - 
                                                                                                                                                                                                                                                                                                                                 1) * (1 - S2)^thet + thet * ((1 - S1)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                                                          S1) - (1 - S1)^(thet - 1) * (1 - S2)^thet * log(1 - S2)) - 
                                                                                                                                                                                                                                                                                                 (((1 - S1)^(thet - 1) + thet * (1 - S1)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                            S1)) * (1 - S2)^thet + (1 - S1)^(thet - 1))) * ((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                               S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            1)/thet
    
    
    
    
    derhfuncm2_dtheta <- ((((1 - S1)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - S1) + 
                             ((1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - 
                                                                                     S2)) * ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                               (1 - S2)^thet)^(1/thet - 2) * (1/thet - 1) - ((1 - S1)^thet + 
                                                                                                                                               (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                 1) * log((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                                                                                                                                            (1 - S2)^thet)/thet^2) * ((1 - S2)^(thet - 1) - (1 - S1)^thet * 
                                                                                                                                                                                                                                        (1 - S2)^(thet - 1)) + ((1 - S1)^thet * (1 - S2)^(thet - 
                                                                                                                                                                                                                                                                                            1) + (1 - S2)^(thet - 1) + thet * ((1 - S2)^(thet - 1) * 
                                                                                                                                                                                                                                                                                                                                 log(1 - S2) - (1 - S1)^thet * (1 - S2)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                                                           S1)) - (((1 - S2)^(thet - 1) + thet * (1 - S2)^(thet - 1) * 
                                                                                                                                                                                                                                                                                                                                                                                                      log(1 - S2)) * (1 - S1)^thet + (1 - S2)^(thet - 1))) * ((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                 S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              1)/thet
    
    
    
    
    dercopcdf_dtheta <-  (((1 - S1)^thet - (-1 + (1 - S1)^thet) *(1 - S2)^thet)^(1/thet)* (log((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet) + (
      thet* ((1 - S1)^thet* (-1 + (1 - S2)^thet)* log(1 - S1) + (-1 + (1 - S1)^thet) *(1 - S2)^thet* log(1 - S2)))/((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)))/thet^2
    
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f) + 1
    
    ###### Copula expressions: 
    copdensity <- (1 - S1)^(-1 + thet)* ((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)^(-2 + 1/thet) *(1 - S2)^(-1 + thet)* (-(-1 + (1 - S1)^thet)* (-1 + (1 - S2)^thet) + thet)
    
    hfunc_m1 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S1)^(thet - 1) * thet - (1 - S1)^(thet - 1) * thet * (1 - S2)^thet))
    
    hfunc_m2 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S2)^(thet - 1) * thet - (1 - S1)^thet * ((1 - S2)^(thet - 1) * thet)))
    
    bit1 <- (1 - S1)^thet
    
    bit2 <- (1 - S2)^thet
    
    copcdf <-  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) + 1 },
  
  name = "Joe copula for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_JoeCopula_180_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    S1_Tilde <- pdffz(S1)
    
    S2_Tilde <- pdffz(S2)
    
    S1 <- pdffz(1 - S1)
    S2 <- pdffz(1 - S2)
    
    ###### Copula expressions: 
    copdensity <- (1 - S1)^(-1 + thet)* ((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)^(-2 + 1/thet) *(1 - S2)^(-1 + thet)* (-(-1 + (1 - S1)^thet)* (-1 + (1 - S2)^thet) + thet)
    
    hfunc_m1 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S1)^(thet - 1) * thet - (1 - S1)^(thet - 1) * thet * (1 - S2)^thet))
    
    hfunc_m2 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S2)^(thet - 1) * thet - (1 - S1)^thet * ((1 - S2)^(thet - 1) * thet)))
    
    bit1 <- (1 - S1)^thet
    
    bit2 <- (1 - S2)^thet
    
    copcdf <-  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz(copdensity)
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S1_Tilde + S2_Tilde - 1 + copcdf)
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t1 = 1.0-S1;
    t2 = t1^(1.0*thet);
    t3 = 1.0-S2;
    t4 = t3^(1.0*thet);
    t5 = t2*t4;
    t6 = t2+t4-t5;
    t8 = 1/thet-2.0;
    t9 = t6^(1.0*t8);
    t10 = thet^2;
    t11 = log(t6);
    t12 = log(t1);
    t13 = t2*t12;
    t14 = log(t3);
    t15 = t4*t14;
    t16 = t13*t4;
    t19 = t5*t14;
    t21 = thet-1.0;
    t27 = t1^(1.0*t21);
    t28 = t3^(1.0*t21);
    t30 = thet-1.0+t2+t4-t5;
    t33 = t9*t27;
    
    dercopdensity_dtheta <- t9*(-1/t10*t11+t8*(t13+t15-t16-t19)/t6)*t27*t28*t30+t33*t12*t28*t30+t33*t28*t14*t30+t33*t28*(1.0+t13+t15-t16-t19);
    
    
    derhfuncm1_dtheta <- ((((1 - S1)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - S1) + 
                             ((1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - 
                                                                                     S2)) * ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                               (1 - S2)^thet)^(1/thet - 2) * (1/thet - 1) - ((1 - S1)^thet + 
                                                                                                                                               (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                 1) * log((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                                                                                                                                            (1 - S2)^thet)/thet^2) * ((1 - S1)^(thet - 1) - (1 - S1)^(thet - 
                                                                                                                                                                                                                                                                        1) * (1 - S2)^thet) + ((1 - S1)^(thet - 1) + (1 - S1)^(thet - 
                                                                                                                                                                                                                                                                                                                                 1) * (1 - S2)^thet + thet * ((1 - S1)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                                                          S1) - (1 - S1)^(thet - 1) * (1 - S2)^thet * log(1 - S2)) - 
                                                                                                                                                                                                                                                                                                 (((1 - S1)^(thet - 1) + thet * (1 - S1)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                            S1)) * (1 - S2)^thet + (1 - S1)^(thet - 1))) * ((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                               S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            1)/thet
    
    
    
    
    derhfuncm2_dtheta <- ((((1 - S1)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - S1) + 
                             ((1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - 
                                                                                     S2)) * ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                               (1 - S2)^thet)^(1/thet - 2) * (1/thet - 1) - ((1 - S1)^thet + 
                                                                                                                                               (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                 1) * log((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                                                                                                                                            (1 - S2)^thet)/thet^2) * ((1 - S2)^(thet - 1) - (1 - S1)^thet * 
                                                                                                                                                                                                                                        (1 - S2)^(thet - 1)) + ((1 - S1)^thet * (1 - S2)^(thet - 
                                                                                                                                                                                                                                                                                            1) + (1 - S2)^(thet - 1) + thet * ((1 - S2)^(thet - 1) * 
                                                                                                                                                                                                                                                                                                                                 log(1 - S2) - (1 - S1)^thet * (1 - S2)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                                                           S1)) - (((1 - S2)^(thet - 1) + thet * (1 - S2)^(thet - 1) * 
                                                                                                                                                                                                                                                                                                                                                                                                      log(1 - S2)) * (1 - S1)^thet + (1 - S2)^(thet - 1))) * ((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                 S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              1)/thet
    
    
    
    
    dercopcdf_dtheta <-  (((1 - S1)^thet - (-1 + (1 - S1)^thet) *(1 - S2)^thet)^(1/thet)* (log((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet) + (
      thet* ((1 - S1)^thet* (-1 + (1 - S2)^thet)* log(1 - S1) + (-1 + (1 - S1)^thet) *(1 - S2)^thet* log(1 - S2)))/((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)))/thet^2
    
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <- - derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <- - derhfuncm2_dtheta
    
    dercopcdf_dtheta        <- dercopcdf_dtheta
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    S1_Tilde <- pdffz(S1)
    S2_Tilde <- pdffz(S2)
    
    S1 <- pdffz(1 - S1)
    S2 <- pdffz(1 - S2)
    
    
    ###### Copula expressions: 
    copdensity <- (1 - S1)^(-1 + thet)* ((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)^(-2 + 1/thet) *(1 - S2)^(-1 + thet)* (-(-1 + (1 - S1)^thet)* (-1 + (1 - S2)^thet) + thet)
    
    hfunc_m1 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S1)^(thet - 1) * thet - (1 - S1)^(thet - 1) * thet * (1 - S2)^thet))
    
    hfunc_m2 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S2)^(thet - 1) * thet - (1 - S1)^thet * ((1 - S2)^(thet - 1) * thet)))
    
    bit1 <- (1 - S1)^thet
    
    bit2 <- (1 - S2)^thet
    
    copcdf <-  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz(copdensity)
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S1_Tilde + S2_Tilde - 1 + copcdf)
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) + 1 },
  
  name = "Joe copula (180° rotation) for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_JoeCopula_90_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S1 <- pdffz(1 - S1)
    
    ###### Copula expressions: 
    copdensity <- (1 - S1)^(-1 + thet)* ((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)^(-2 + 1/thet) *(1 - S2)^(-1 + thet)* (-(-1 + (1 - S1)^thet)* (-1 + (1 - S2)^thet) + thet)
    
    hfunc_m1 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S1)^(thet - 1) * thet - (1 - S1)^(thet - 1) * thet * (1 - S2)^thet))
    
    hfunc_m2 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S2)^(thet - 1) * thet - (1 - S1)^thet * ((1 - S2)^(thet - 1) * thet)))
    
    bit1 <- (1 - S1)^thet
    
    bit2 <- (1 - S2)^thet
    
    copcdf <-  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S2 - copcdf)
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t1 = 1.0-S1;
    t2 = t1^(1.0*thet);
    t3 = 1.0-S2;
    t4 = t3^(1.0*thet);
    t5 = t2*t4;
    t6 = t2+t4-t5;
    t8 = 1/thet-2.0;
    t9 = t6^(1.0*t8);
    t10 = thet^2;
    t11 = log(t6);
    t12 = log(t1);
    t13 = t2*t12;
    t14 = log(t3);
    t15 = t4*t14;
    t16 = t13*t4;
    t19 = t5*t14;
    t21 = thet-1.0;
    t27 = t1^(1.0*t21);
    t28 = t3^(1.0*t21);
    t30 = thet-1.0+t2+t4-t5;
    t33 = t9*t27;
    
    dercopdensity_dtheta <- t9*(-1/t10*t11+t8*(t13+t15-t16-t19)/t6)*t27*t28*t30+t33*t12*t28*t30+t33*t28*t14*t30+t33*t28*(1.0+t13+t15-t16-t19);
    
    
    derhfuncm1_dtheta <- ((((1 - S1)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - S1) + 
                             ((1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - 
                                                                                     S2)) * ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                               (1 - S2)^thet)^(1/thet - 2) * (1/thet - 1) - ((1 - S1)^thet + 
                                                                                                                                               (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                 1) * log((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                                                                                                                                            (1 - S2)^thet)/thet^2) * ((1 - S1)^(thet - 1) - (1 - S1)^(thet - 
                                                                                                                                                                                                                                                                        1) * (1 - S2)^thet) + ((1 - S1)^(thet - 1) + (1 - S1)^(thet - 
                                                                                                                                                                                                                                                                                                                                 1) * (1 - S2)^thet + thet * ((1 - S1)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                                                          S1) - (1 - S1)^(thet - 1) * (1 - S2)^thet * log(1 - S2)) - 
                                                                                                                                                                                                                                                                                                 (((1 - S1)^(thet - 1) + thet * (1 - S1)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                            S1)) * (1 - S2)^thet + (1 - S1)^(thet - 1))) * ((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                               S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            1)/thet
    
    
    
    
    derhfuncm2_dtheta <- ((((1 - S1)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - S1) + 
                             ((1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - 
                                                                                     S2)) * ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                               (1 - S2)^thet)^(1/thet - 2) * (1/thet - 1) - ((1 - S1)^thet + 
                                                                                                                                               (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                 1) * log((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                                                                                                                                            (1 - S2)^thet)/thet^2) * ((1 - S2)^(thet - 1) - (1 - S1)^thet * 
                                                                                                                                                                                                                                        (1 - S2)^(thet - 1)) + ((1 - S1)^thet * (1 - S2)^(thet - 
                                                                                                                                                                                                                                                                                            1) + (1 - S2)^(thet - 1) + thet * ((1 - S2)^(thet - 1) * 
                                                                                                                                                                                                                                                                                                                                 log(1 - S2) - (1 - S1)^thet * (1 - S2)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                                                           S1)) - (((1 - S2)^(thet - 1) + thet * (1 - S2)^(thet - 1) * 
                                                                                                                                                                                                                                                                                                                                                                                                      log(1 - S2)) * (1 - S1)^thet + (1 - S2)^(thet - 1))) * ((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                 S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              1)/thet
    
    
    
    
    dercopcdf_dtheta <-  (((1 - S1)^thet - (-1 + (1 - S1)^thet) *(1 - S2)^thet)^(1/thet)* (log((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet) + (
      thet* ((1 - S1)^thet* (-1 + (1 - S2)^thet)* log(1 - S1) + (-1 + (1 - S1)^thet) *(1 - S2)^thet* log(1 - S2)))/((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)))/thet^2
    
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <-   derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <- - derhfuncm2_dtheta
    
    dercopcdf_dtheta        <-  - dercopcdf_dtheta
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S1 <- pdffz(1 - S1)
    
    
    ###### Copula expressions: 
    copdensity <- (1 - S1)^(-1 + thet)* ((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)^(-2 + 1/thet) *(1 - S2)^(-1 + thet)* (-(-1 + (1 - S1)^thet)* (-1 + (1 - S2)^thet) + thet)
    
    hfunc_m1 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S1)^(thet - 1) * thet - (1 - S1)^(thet - 1) * thet * (1 - S2)^thet))
    
    hfunc_m2 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S2)^(thet - 1) * thet - (1 - S1)^thet * ((1 - S2)^(thet - 1) * thet)))
    
    bit1 <- (1 - S1)^thet
    
    bit2 <- (1 - S2)^thet
    
    copcdf <-  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(hfunc_m1)
    
    hfunc_m2    <- pdffz(1 - hfunc_m2)
    
    copcdf      <- pdffz(S2 - copcdf)
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) + 1 },
  
  name = "Joe copula (90° rotation) for bivariate time-to-event data with right-censoring scheme")
}

BivAFT_RC_JoeCopula_270_RhoSoloFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    S1    <- pdffz(y[,1])
    
    S2    <- pdffz( y[,4] )
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S2 <- pdffz(1 - S2)
    
    ###### Copula expressions: 
    copdensity <- (1 - S1)^(-1 + thet)* ((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)^(-2 + 1/thet) *(1 - S2)^(-1 + thet)* (-(-1 + (1 - S1)^thet)* (-1 + (1 - S2)^thet) + thet)
    
    hfunc_m1 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S1)^(thet - 1) * thet - (1 - S1)^(thet - 1) * thet * (1 - S2)^thet))
    
    hfunc_m2 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S2)^(thet - 1) * thet - (1 - S1)^thet * ((1 - S2)^(thet - 1) * thet)))
    
    bit1 <- (1 - S1)^thet
    
    bit2 <- (1 - S2)^thet
    
    copcdf <-  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)
    ################################
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(hfunc_m2)
    
    copcdf      <- pdffz(S1 - copcdf)
    #### Derivatives of copula expressions w.r.t. copula parameter
    #### some other bits: BIG credit to GJRM!
    
    t1 = 1.0-S1;
    t2 = t1^(1.0*thet);
    t3 = 1.0-S2;
    t4 = t3^(1.0*thet);
    t5 = t2*t4;
    t6 = t2+t4-t5;
    t8 = 1/thet-2.0;
    t9 = t6^(1.0*t8);
    t10 = thet^2;
    t11 = log(t6);
    t12 = log(t1);
    t13 = t2*t12;
    t14 = log(t3);
    t15 = t4*t14;
    t16 = t13*t4;
    t19 = t5*t14;
    t21 = thet-1.0;
    t27 = t1^(1.0*t21);
    t28 = t3^(1.0*t21);
    t30 = thet-1.0+t2+t4-t5;
    t33 = t9*t27;
    
    dercopdensity_dtheta <- t9*(-1/t10*t11+t8*(t13+t15-t16-t19)/t6)*t27*t28*t30+t33*t12*t28*t30+t33*t28*t14*t30+t33*t28*(1.0+t13+t15-t16-t19);
    
    
    derhfuncm1_dtheta <- ((((1 - S1)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - S1) + 
                             ((1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - 
                                                                                     S2)) * ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                               (1 - S2)^thet)^(1/thet - 2) * (1/thet - 1) - ((1 - S1)^thet + 
                                                                                                                                               (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                 1) * log((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                                                                                                                                            (1 - S2)^thet)/thet^2) * ((1 - S1)^(thet - 1) - (1 - S1)^(thet - 
                                                                                                                                                                                                                                                                        1) * (1 - S2)^thet) + ((1 - S1)^(thet - 1) + (1 - S1)^(thet - 
                                                                                                                                                                                                                                                                                                                                 1) * (1 - S2)^thet + thet * ((1 - S1)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                                                          S1) - (1 - S1)^(thet - 1) * (1 - S2)^thet * log(1 - S2)) - 
                                                                                                                                                                                                                                                                                                 (((1 - S1)^(thet - 1) + thet * (1 - S1)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                            S1)) * (1 - S2)^thet + (1 - S1)^(thet - 1))) * ((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                               S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            1)/thet
    
    
    
    
    derhfuncm2_dtheta <- ((((1 - S1)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - S1) + 
                             ((1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet) * log(1 - 
                                                                                     S2)) * ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                               (1 - S2)^thet)^(1/thet - 2) * (1/thet - 1) - ((1 - S1)^thet + 
                                                                                                                                               (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                 1) * log((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * 
                                                                                                                                                                                                            (1 - S2)^thet)/thet^2) * ((1 - S2)^(thet - 1) - (1 - S1)^thet * 
                                                                                                                                                                                                                                        (1 - S2)^(thet - 1)) + ((1 - S1)^thet * (1 - S2)^(thet - 
                                                                                                                                                                                                                                                                                            1) + (1 - S2)^(thet - 1) + thet * ((1 - S2)^(thet - 1) * 
                                                                                                                                                                                                                                                                                                                                 log(1 - S2) - (1 - S1)^thet * (1 - S2)^(thet - 1) * log(1 - 
                                                                                                                                                                                                                                                                                                                                                                                           S1)) - (((1 - S2)^(thet - 1) + thet * (1 - S2)^(thet - 1) * 
                                                                                                                                                                                                                                                                                                                                                                                                      log(1 - S2)) * (1 - S1)^thet + (1 - S2)^(thet - 1))) * ((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                 S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^(1/thet - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              1)/thet
    
    
    
    
    dercopcdf_dtheta <-  (((1 - S1)^thet - (-1 + (1 - S1)^thet) *(1 - S2)^thet)^(1/thet)* (log((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet) + (
      thet* ((1 - S1)^thet* (-1 + (1 - S2)^thet)* log(1 - S1) + (-1 + (1 - S1)^thet) *(1 - S2)^thet* log(1 - S2)))/((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)))/thet^2
    
    
    ################################                      ROTATIONS TERMS OF DERIVATIVES
    dercopdensity_dtheta    <-  dercopdensity_dtheta
    
    derhfuncm1_dtheta       <-  - derhfuncm1_dtheta
    
    derhfuncm2_dtheta       <-    derhfuncm2_dtheta
    
    dercopcdf_dtheta        <-   - dercopcdf_dtheta
    
    #### some checks again: 
    dercopdensity_dtheta  <- dvffz( dercopdensity_dtheta )
    
    derhfuncm1_dtheta     <- dvffz( derhfuncm1_dtheta )
    
    derhfuncm2_dtheta     <- dvffz( derhfuncm2_dtheta )
    
    dercopcdf_dtheta      <- dvffz( dercopcdf_dtheta )
    
    ################################
    uncens <- ( dvffz( copdensity ) )
    
    twocens <- ( pdffz( hfunc_m1 )  )
    
    onecens <- ( pdffz( hfunc_m2 )  )
    
    bothcens <- pdffz( copcdf )
    
    
    ngr <- y[,3] * y[,6] *( uncens^(-1) * ( dercopdensity_dtheta ) * exp(f) ) +
      #
      y[,3] * ( 1 - y[,6] ) * ( twocens^(-1) * ( derhfuncm1_dtheta ) * exp(f) ) +
      #
      ( 1 - y[,3] ) * y[,6] * ( onecens^(-1) * ( derhfuncm2_dtheta ) * exp(f) )  +
      #
      ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens^(-1) * ( dercopcdf_dtheta ) * exp(f) )
    
    
    return(ngr)
    
  },
  
  loss = function(y, f){ 
    
    S1    <- pdffz( y[,1] )
    PDF1  <- dvffz( y[,2] )
    PDF1  <- ifelse(PDF1 < 1e-40, 1e-40, PDF1)
    
    S2    <- pdffz( y[,4] )
    PDF2  <- dvffz( y[,5] )
    PDF2  <- ifelse(PDF2 < 1e-40, 1e-40, PDF2)
    
    thet <- exp(f) + 1
    
    ####### ROTATION TERMS (used outside of the copula function)
    ## Into the copula goes: 1 - u
    S2 <- pdffz(1 - S2)
    
    
    ###### Copula expressions: 
    copdensity <- (1 - S1)^(-1 + thet)* ((1 - S1)^thet - (-1 + (1 - S1)^thet)* (1 - S2)^thet)^(-2 + 1/thet) *(1 - S2)^(-1 + thet)* (-(-1 + (1 - S1)^thet)* (-1 + (1 - S2)^thet) + thet)
    
    hfunc_m1 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S1)^(thet - 1) * thet - (1 - S1)^(thet - 1) * thet * (1 - S2)^thet))
    
    hfunc_m2 <- ((1 - S1)^thet + (1 - S2)^thet - (1 - S1)^thet * (1 - S2)^thet)^((1/thet) - 1) * ((1/thet) * ((1 - S2)^(thet - 1) * thet - (1 - S1)^thet * ((1 - S2)^(thet - 1) * thet)))
    
    bit1 <- (1 - S1)^thet
    
    bit2 <- (1 - S2)^thet
    
    copcdf <-  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)
    
    ################################                      ROTATIONS TERMS
    copdensity  <- dvffz( copdensity )
    
    hfunc_m1    <- pdffz(1 - hfunc_m1)
    
    hfunc_m2    <- pdffz(hfunc_m2)
    
    copcdf      <- pdffz(S1 - copcdf)
    
    ### SMALL CHECKS:
    copdensity <- dvffz( copdensity )
    
    hfunc_m1  <- pdffz( hfunc_m1 )
    
    hfunc_m2  <- pdffz( hfunc_m2 )
    
    copcdf    <- pdffz( copcdf )
    
    
    ################################
    
    
    uncens    <- ( log( PDF1 ) + log( PDF2 ) + log(  copdensity  ) )
    
    twocens   <- ( log( PDF1 ) +  log(  hfunc_m1  )  )
    
    onecens   <- ( log( PDF2 ) +  log(  hfunc_m2  )  )
    
    bothcens  <- log(  copcdf  )
    
    
    
    loglik <-  (
      # Uncensored
      y[,3] * y[,6] * ( uncens ) + 
        # y2 censored
        y[,3] * ( 1 - y[,6] ) * ( twocens ) +
        # y1 censored
        ( 1 - y[,3] ) * y[,6] * ( onecens ) + 
        # both censored
        ( 1 - y[,3] ) * ( 1 - y[,6] ) * ( bothcens )  )
    
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  response = function(f){ exp(f) + 1 },
  
  name = "Joe copula (270° rotation) for bivariate time-to-event data with right-censoring scheme")
}


