################################################################################
### DESCRIPTION
### Gauss_Cop() defines the function to create a Families object from gamboostLSS
### for the Gaussian copula with arbitrary marginal distributions.
### Marginal functions (like the ones implemented in the Marginals folder) 
### are handed over to the arguments marg1 and marg2. 
### Over the course of the function run, first generic bodies and arguments of the functions 
### for the loss, risk and negative gradient are defined for each parameter submodel, 
### i.e. the marginal distribution parameters and 
### the copula dependence parameter (also via the rlang package). 
### Thereafter, the parameter specific functions and Family objects (see mboost) 
### are created via the family_gen() function 
### (see Marginals folder) for each parameter sub-model. 
### Finally, a gamboostLSS Families object is created that suits the gamboostLSS 
### boosting algorithm implemented in the gamboostLSS package.


### load libraries 
library(mboost)
library(gamboostLSS)
library(rlang)


# load Marginal functions

source("Simulations/ContinuousMargins/family_gen.R")
source("Simulations/ContinuousMargins/Marginal_LogNormal.R")
source("Simulations/ContinuousMargins/Marginal_LogLogistic.R")
source("Simulations/ContinuousMargins/Marginal_Weibull.R")

# y[, 1] <- y1
# y[, 2] <- DELTA1
# y[, 3] <- DELTA2
# y[, 4] <- S2
# y[, 5] <- PDF2

Gumbel_Copula_SemiCompeting <- function(marg1 = NULL, 
                                     mu1 = NULL, sigma1 = NULL, nu1 = NULL, tau1 = NULL,
                                     rho = NULL,
                                     stabilization = c("none", "MAD", "L2")){
  
  ################################################################################
  ########################## marginal checks #####################################
  
  if(is.null(marg1)) stop("First marginal distribution for non-terminal event is not defined.")
  
  if(!(marg1 %in% c("NO", "LOGNO", "GA", "LOGLOG", "EXP", "WEI", "NORMAL"))) stop("First Marginal distribution not available.")
  
  
  ################################################################################
  #################### calling the marginal distributions ########################
  
  if(marg1 == "LOGNO") {
    marg1 <- LogNormal_Mar(loc = 1, offset_sigma = sigma1)
    
  } else if(marg1 == "LOGLOG") {
    marg1 <- LogLogistic_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1)
    
  } else if(marg1 == "WEI") {
    marg1 <- Weibull_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1)
  }
  
  
  
  
  
  ################################################################################
  ##################### check offsets of copula parameter ########################
  
  
  if ((!is.null(rho) && rho <= 1))
    stop(sQuote("rho"), " must be equal to or greater than 1.")
  
  
  ################################################################################
  ############################## Helpers #########################################
  ################### To be deleted when integration into package ################
  
  check_stabilization <- function(stabilization = c("none", "MAD", "L2")) {
    stabilization <- match.arg(stabilization)
    ## check if old stabilization interface is used and issue a warning
    if (getOption("gamboostLSS_stab_ngrad")) {
      warning("Usage of ", sQuote("options(gamboostLSS_stab_ngrad = TRUE)"),
              " is deprecated.\n", "Use argument ", sQuote("stabilization"),
              " in the fitting family. See ?Families for details.")
      if (stabilization == "none")
        warning(sQuote("stabilization"), " is set to ", dQuote("MAD"))
    }
    stabilization
  }
  
  
  ################################################################################
  ########################## Check stabilization #################################
  
  stabilization <- check_stabilization(stabilization)
  
  
  
  
  ################################################################################
  ############### Explanation of the automatic function creation: ################
  
  ### Aim: Automatically create loss, risk, gradient, offset and check_y function,
  ###      which are necessary for the family object construction (see mboost).
  
  ### Remark: Each parameter (marginals and copulas) has its individual functions.
  ###         The loss, risk and gradients consist partly of copula and partly of marginal elements,
  ###         which need to be combined appropriately.
  
  ### Procedure: All functions are created automatically. 
  ###            To do so, it is necessary to define appropriate arguments and appropriate bodies.
  ###            In the following, firstly the functions arguments are created, 
  ###            secondly the functions bodies (at least the generic functions) are created 
  ###            and finally for each parameter the appropriate functions and family object are constructed (via the family_gen function).
  
  ### for code implementation and the idea see body() and formals() function.
  
  
  
  ################################################################################
  ####################### Arguments for function creation ########################
  
  
  # Note: I decided to put the arguments creation in the copula function (and not in the family_gen function)
  #       because of potential problems with parameter specific loss creation. 
  #       Also conceptually it makes sense to have the whole function creation procedure in one place.
  
  # code explanation: Create the arguments for the loss, risk, gradient, offset and check_y function.
  # The risk, gradient, offset and check_y arguments are equal for all parameters.
  # For the loss function parameters specific arguments need to be created via the args_loss_creator. 
  # for further information see also https://stackoverflow.com/questions/17751862/create-a-variable-length-alist
  
  
  # arguments for gradient and risk functions
  args_grad_risk <- c("y", "f", "w")
  l_args_grad_risk <- rep(list(expr()), length(args_grad_risk)) 
  names(l_args_grad_risk) <- args_grad_risk
  l_args_grad_risk[["w"]] <- 1
  
  # arguments for offset functions
  args_offset <- c("y", "w")
  l_args_offset <- rep(list(expr()), length(args_offset)) 
  names(l_args_offset) <- args_offset
  
  # arguments for check_y functions
  args_ckeck_y <- c("y")
  l_args_check_y <- rep(list(expr()), length(args_ckeck_y)) 
  names(l_args_check_y) <- args_ckeck_y
  
  # generic arguments for loss functions
  args_loss_gen <- c("y", 
                     paste(marg1$parameter_names, "1", sep = ""),
                     "rho")
  
  # counter required to create parameter specific loss arguments via args_loss_creator()
  ind_arg <- 2
  
  # creates the parameter specific loss arguments
  args_loss_creator <- function(arg_names, ind_arg){    
    arg_names[ind_arg] <- "f"
    l_args <- rep(list(expr()), length(arg_names))
    names(l_args) <- arg_names
    ind_arg <<- ind_arg + 1     
    return(list(l_args,
                arg_names))
  }
  
  
  
  
  ################################################################################
  ############################# generic functions ################################
  
  # rho parameter definition: 
  # 1. rho_simp means simple rho for marginal parameters  
  # 2. rho_gen means rho generic for the f-expression 
  rho_simp <- expr(rho)
  rho_gen <- expr({   exp( f ) + 1 })
  
  # generic functions for the parameter specific creation of loss, risk and gradients 
  
  # y[, 1] <- y1
  # y[, 2] <- DELTA1
  # y[, 3] <- DELTA2
  # y[, 4] <- S2
  # y[, 5] <- PDF2
  gradient_marg_gen <- function(cdf1, dercdf1, derlpdf1){
    
    expr({ 
      
      S1 <- 1 - pdffz(!!cdf1)
      
      S2 <- pdffz(y[,4])
      
      dF1 <- dvffz(!!dercdf1)
      
      dlogpdf1 <- dvffz(!!derlpdf1)
      
      thet <- rho
      
      ###### Copula expressions: 
      copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
      
      twocens <- (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
      
      onecens <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
      
      bothcens <- exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
      ################################
      # Derivatives of copula expressions w.r.t. margin 1:
      
      t3 = log(S1);
      t4 = (-t3)^thet;
      t5 = log(S2);
      t6 = (-t5)^thet
      t7 = t4+t6;
      t8 = 1/thet;
      t9 = t7^t8
      t11 = S1^2;
      t12 = 1/t11;
      t13 = 1/t3;
      t15 = 1/t7;
      t18 = exp(-t9);
      t19 = 1/S2;
      t21 = -1 + t8;
      t22 = (t7)^(2*t21)
      t24 = thet - 1
      t25 = (t3*t5)^t24
      t27 = t7^-t8
      t28 = t24*t27;
      t29 = 1 + t28;
      t30 = t22*t25*t29;
      t33 = t18*t12;
      t36 = t19*t22;
      
      # derivative copula density w.r.t. margin 1
      dcopdensity_dS1 <-    dvffz(  -t9*t4*t12*t13*t15*t18*t19*t30-t33*t19*t30+2.0*t33*t36*t21*t4*thet*t13*t15*t25*t29+t33*t36*t25*t24*t13*t29-t33*t36*t25*t28*t4*t13*t15 );
      
      
      # derivative hfunction w.r.t. margin 1 w.r.t. margin 1
      dhfuncm1_dS1 <- dvffz(  (1/(S1^2))*exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-2 + thet) *((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* ((-log(S1))^thet* (log(S1) + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet)) + (1 - thet + log(S1))* (-log(S2))^thet) )
      
      # copula density:
      dhfuncm2_dS1 <- dvffz( (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2) )         
      
      # derivative copula CDF w.r.t. margin 1: dcopcdf_dS1
      dcopcdf_dS1 <- pdffz(  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1  )
      
      ##############
      
      
      ngr <- y[,2] * y[,3] * ( dlogpdf1 + dvffz(copdensity)^(-1) * dcopdensity_dS1 * (- dF1 ) )   +
        #     
        y[,2] * ( 1 - y[,3] ) * ( dlogpdf1 + pdffz(twocens)^(-1) * dhfuncm1_dS1 * (- dF1 )  ) + 
        
        ( 1 - y[,2] ) * y[,3] * (              pdffz(onecens)^(-1) * dhfuncm2_dS1 * (- dF1 )  ) + 
        
        ( 1 - y[,2] ) * ( 1 - y[,3] ) * ( pdffz( bothcens )^(-1) * dcopcdf_dS1 * (- dF1 )  )
      
      return(ngr)
    })
    
  }
  
  
  gradient_cop_gen <- function(cdf1, rho){
    
    expr({
      
      S1 <- 1 - pdffz(!!cdf1)
      
      S2 <- pdffz(y[,4])
      
      thet <- !!rho
      
      ###### Copula expressions: 
      copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
      
      hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
      
      hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
      
      copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
      ################################
      #### Derivatives of copula expressions w.r.t. copula parameter
      #### some other bits:
      
      
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
      
      
      ################################
      uncens <- ( dvffz( copdensity ) )
      
      twocens <- ( pdffz( hfunc_m1 )  )
      
      onecens <- ( pdffz( hfunc_m2 )  )
      
      bothcens <- pdffz( copcdf )
      
      
      ngr <- y[,2] * y[,3] *( uncens^(-1) * dvffz( dercopdensity_dtheta ) * exp(f) ) +
        #
        y[,2] * ( 1 - y[,3] ) * ( twocens^(-1) * dvffz( derhfuncm1_dtheta ) * exp(f) ) +
        #
        ( 1 - y[,2] ) * y[,3] * ( onecens^(-1) * dvffz( derhfuncm2_dtheta ) * exp(f) )  +
        #
        ( 1 - y[,2] ) * ( 1 - y[,3] ) * ( bothcens^(-1) * dvffz( dercopcdf_dtheta ) * exp(f) )
      
      
      return(ngr)
    })
  }
  
  
  loss_gen <- function(pdf1, cdf1, rho){
    
    expr({
      
      S1 <- 1 - pdffz(!!cdf1)
      
      f1 <- !!pdf1
      
      f1 <- ifelse(f1 < 1e-40, 1e-40, f1)
      
      f1 <- dvffz(f1)
      
      
      S2 <- pdffz(y[,4])
      
      f2 <- y[,5]
      
      f2 <- ifelse(f2 < 1e-40, 1e-40, f2)
      
      f2 <- dvffz(f2)
      
      thet <- !!rho
      
      ###### Copula expressions: 
      copdensity <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet) *(-1 + thet + ((-log(S1))^thet + (-log(S2))^thet)^(1/thet))* ((-log(S1))^thet + (-log(S2))^thet)^(-2 + 1/thet)* (-log(S2))^(-1 + thet))/(S1 *S2)    
      
      hfunc_m1 <-  (exp(-((-log(S1))^thet + (-log(S2))^thet)^((1/thet)))* (-log(S1))^(-1 + thet)* ((-log(S1))^thet + (-log(S2))^thet)^(-1 + 1/thet))/S1
      
      hfunc_m2 <- (exp(-((-log(S2))^thet + (-log(S1))^thet)^((1/thet)))* (-log(S2))^(-1 + thet)* ((-log(S2))^thet + (-log(S1))^thet)^(-1 + 1/thet))/S2
      
      copcdf <-  exp( - ( (-log(S1))^thet + (-log(S2))^thet )^(1/thet) )
      ################################
      
      
      uncens <- ( log( f1 ) + log( f2 ) + log( dvffz( copdensity ) ) )
      
      twocens <- ( log( f1 ) +  log( pdffz( hfunc_m1 ) )  )
      
      onecens <- ( log( f2 ) +  log( pdffz( hfunc_m2 ) )  )
      
      bothcens <- log( pdffz( copcdf ) )
      
      
      return( 
        -(
          # Uncensored
          y[,2] * y[,3] * ( uncens ) + 
            # y2 censored
            y[,2] * ( 1 - y[,3] ) * ( twocens ) +
            # y1 censored
            ( 1 - y[,2] ) * y[,3] * ( onecens ) + 
            # both censored
            ( 1 - y[,2] ) * ( 1 - y[,3] ) * ( bothcens ) 
        ))
      
    })
  }
  
  
  risk_gen <- function(param){
    
    a <- param
    b <- paste("=", param)
    param <- paste(a, b, collapse = ", ")
    
    loss <- parse_expr(paste0("loss(", param, ")"))  
    
    expr(sum(w*(!!loss)))
  }
  
  
  
  
  
  ################################################################################
  ############### initializing the gamboostLSS Families object ###################
  
  
  name_fam <- paste("Gumbel copula for bivariate right-censored time-to-event data (semi-competing risks):", marg1$marg_name)
  
  ### not yet required, rather for prediction, I think.
  GaussCopFam <- Families(name = name_fam)
  
  
  
  ################################################################################
  ############# creating the Family-objects for the marginals ####################
  
  ### first marginal
  
  for(i in marg1$parameter_names){
    
    #i <- marg1$parameter_names[1]
    
    # body and arguments for loss
    
    loss_body <- loss_gen(pdf1 = marg1[[i]]$pdf,
                          cdf1 = marg1[[i]]$cdf,
                          rho = rho_simp)
    
    loss_args <- args_loss_creator(arg_names = args_loss_gen,
                                   ind_arg = ind_arg)
    
    
    # body and arguments for risk
    
    risk_body <- risk_gen(loss_args[[2]]) 
    risk_args <- l_args_grad_risk
    
    
    # body and arguments for gradient
    
    grad_body <- gradient_marg_gen(derlpdf1 = marg1[[i]]$derlpdf, 
                                   cdf1 = marg1[[i]]$cdf, 
                                   dercdf1 = marg1[[i]]$dercdf
    ) 
    
    grad_body <- exprs(!!grad_body,
                       ngr <- stabilize_ngradient(ngr, w = w, stabilization),
                       return(ngr))
    
    grad_args <- l_args_grad_risk
    
    
    
    # creating the Family object
    fam_obj <- family_gen(mu1 = mu1, sigma1 = sigma1, nu1 = nu1, tau1 = tau1,
                          rho = rho,
                          stabilization = stabilization,
                          loss_body = loss_body, loss_args = loss_args[[1]], 
                          risk_body = risk_body, risk_args = risk_args,
                          grad_body = grad_body, grad_args = grad_args,
                          offset_body = marg1$offset[[i]], offset_args = l_args_offset,
                          response = marg1$response[[i]],
                          name = marg1$name[[i]],
                          check_y_body = marg1$check_y[[i]], check_y_args = l_args_check_y)
    
    
    # saving the parameter specific family object in overall families object
    GaussCopFam[[paste(i, "1", sep = "")]] <- fam_obj
    
    # removing family object
    rm(fam_obj)
    
  }
  
  
  ################################################################################
  ############# creating the Family-objects for the Copula Para ##################
  
  # body and arguments for loss
  loss_body <- loss_gen(pdf1 = marg1[["generic"]]$pdf,
                        cdf1 = marg1[["generic"]]$cdf,
                        rho = rho_gen
  )
  
  loss_args <- args_loss_creator(arg_names = args_loss_gen,
                                 ind_arg = ind_arg)
  
  
  # body and arguments for risk
  risk_body <- risk_gen(loss_args[[2]])
  
  risk_args <- l_args_grad_risk
  
  
  # body and arguments for gradient
  grad_body <- gradient_cop_gen(cdf1 = marg1[["generic"]]$cdf, 
                                rho = rho_gen 
  )
  
  grad_body <- exprs(!!grad_body,
                     ngr <- stabilize_ngradient(ngr, w = w, stabilization),
                     return(ngr))
  
  grad_args <- l_args_grad_risk
  
  
  # definition of offset, response, check_y and name for copula parameter
  offset_cop <- expr({
    if (!is.null(rho)){
      RET <- log(rho)
    } else {
      RET <- 0.01 # rho = 0/(sqrt(1 + 0)) = 0
      # pear_cor <- wdm(x = y[,1], y = y[,2], method = "pearson", weights = w, remove_missing = F) 
      # RET <- pear_cor/sqrt(1-pear_cor^2)
    }
    return(RET)
  })
  
  # offset_cop <- expr({
  #   if (!is.null(rho)){
  #     RET <- rho/sqrt(1 - rho^2)
  #   } else {
  #     if (is.null(rho))
  #     RET <- 0
  #   }
  #   return(RET)
  # })
  
  response_cop <- function(f) exp(f) + 1
  
  check_y_cop <- expr({
    y
  })
  
  name_cop <- "Gumbel copula for bivariate time-to-event semi-competing risks responses: rho(log(.-1) link)"                                       
  
  
  # creating the family object
  fam_obj <- family_gen(mu1 = mu1, sigma1 = sigma1, nu1 = nu1, tau1 = tau1,
                        rho = rho,
                        stabilization = stabilization,
                        loss_body = loss_body, loss_args = loss_args[[1]], 
                        risk_body = risk_body, risk_args = risk_args,
                        grad_body = grad_body, grad_args = grad_args,
                        offset_body = offset_cop, offset_args = l_args_offset,
                        response = response_cop,
                        name = name_cop,
                        check_y_body = check_y_cop, check_y_args = l_args_check_y)
  
  
  # saving the parameter specific family object in overall families object
  GaussCopFam[["rho"]] <- fam_obj
  
  # removing family object
  rm(fam_obj)
  
  
  ### return final Families object to boost with
  return(GaussCopFam)
  
}

# Gauss_Cop <- Gauss_Cop(marg1 = "LOGLOG", marg2 = "LOGLOG")
# Gauss_Cop$rho@offset

