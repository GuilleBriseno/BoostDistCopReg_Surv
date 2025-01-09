#### Fit independent COX models: 

library("VineCopula")
library("VGAM")
library("gamboostLSS")
library("mvtnorm")
library("survival")
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



############################################################################# DEFINE HYPERPARAMETERS: 
boost_nu <- 0.005


the_stabilization <- "L2"


###### These are for initialisation: 
TryMstop_TERMINAL_COX <- 10000


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
###############################++++++++++++++++++++++++++++++ FIT MODEL WITH CLINICAL COVARIATES AND SNPs
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

### re-scale the genetic covariates:
JUST_SNPS <- scale(as.matrix( JUST_SNPS ), center = TRUE, scale = TRUE)

dim(JUST_SNPS)



SCR_Terminal_Data           <- data.frame(terminal_event = SCR_DATA_FITTING$terminal_event, 
                                          terminal_status = SCR_DATA_FITTING$terminal_status,
                                          residual_tumor_size = SCR_DATA_FITTING$residual_tumor_size, 
                                          TumorStage_2 = SCR_DATA_FITTING$TumorStage_2,
                                          TumorStage_3 = SCR_DATA_FITTING$TumorStage_3,
                                          TumorStage_4 = SCR_DATA_FITTING$TumorStage_4,
                                          #
                                          JUST_SNPS[observations_FITTING, ] ) 


SCR_Terminal_Data_External  <- data.frame(terminal_event = SCR_DATA_EXTERNAL$terminal_event, 
                                          terminal_status = SCR_DATA_EXTERNAL$terminal_status,
                                          residual_tumor_size = SCR_DATA_EXTERNAL$residual_tumor_size, 
                                          TumorStage_2 = SCR_DATA_EXTERNAL$TumorStage_2,
                                          TumorStage_3 = SCR_DATA_EXTERNAL$TumorStage_3,
                                          TumorStage_4 = SCR_DATA_EXTERNAL$TumorStage_4,
                                          #
                                          JUST_SNPS[observations_EXTERNAL, ] ) 


head(colnames(SCR_Terminal_Data), 10)
tail(colnames(SCR_Terminal_Data), 10)


SCR_Terminal_Data$residual_tumor_size <- factor( SCR_Terminal_Data$residual_tumor_size )

SCR_Terminal_Data$TumorStage_2        <- factor( SCR_Terminal_Data$TumorStage_2 )
SCR_Terminal_Data$TumorStage_3        <- factor( SCR_Terminal_Data$TumorStage_3 )
SCR_Terminal_Data$TumorStage_4        <- factor( SCR_Terminal_Data$TumorStage_4 )


SCR_Terminal_Data_External$residual_tumor_size <- factor( SCR_Terminal_Data_External$residual_tumor_size )

SCR_Terminal_Data_External$TumorStage_2 <- factor( SCR_Terminal_Data_External$TumorStage_2 )
SCR_Terminal_Data_External$TumorStage_3 <- factor( SCR_Terminal_Data_External$TumorStage_3 )
SCR_Terminal_Data_External$TumorStage_4 <- factor( SCR_Terminal_Data_External$TumorStage_4 )

######################################################################################################### FOR MODEL OF GENETIC COVARIATES
terminal_formula_COX <- formula(Surv(terminal_event, terminal_status) ~ . )


#########################################################################################################
#########################################################################################################
#########################################################################################################     COX MODEL
#########################################################################################################
#########################################################################################################
#########################################################################################################


TERMINAL_EVENT_COX <- glmboost(terminal_formula_COX, 
                               SCR_Terminal_Data, 
                               family = CoxPH(), 
                               weights = weights_mstop, 
                               control = boost_control(mstop = TryMstop_TERMINAL_COX, 
                                                       nu = boost_nu,
                                                       risk = "oobag", 
                                                       trace = TRUE), 
                               center = FALSE)



plot(risk(TERMINAL_EVENT_COX), type = "l")
which.min(risk(TERMINAL_EVENT_COX))

COX_TERMINAL_MSTOP_OPT_OOBAG <- which.min(risk(TERMINAL_EVENT_COX))

TERMINAL_EVENT_COX <- glmboost(terminal_formula_COX, 
                               SCR_Terminal_Data, 
                               family = CoxPH(), 
                               control = boost_control(mstop = COX_TERMINAL_MSTOP_OPT_OOBAG, 
                                                       nu = boost_nu,
                                                       trace = TRUE), 
                               center = FALSE)


length(coef(TERMINAL_EVENT_COX))

COX_TERMINAL_MSTOP_OPT_OOBAG

names(coef(TERMINAL_EVENT_COX))


##################################################################. MARGIN 1:
#########################################################################################################
#########################################################################################################
#########################################################################################################
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


###### FORMULA
nonterminal_formula_COX <- formula(Surv(nonterminal_event, nonterminal_status) ~ . )



NON_TERMINAL_EVENT_COX <- glmboost(nonterminal_formula_COX, 
                               SCR_Nonterminal_Data, 
                               family = CoxPH(), 
                               weights = weights_mstop, 
                               control = boost_control(mstop = TryMstop_TERMINAL_COX, 
                                                       nu = boost_nu,
                                                       risk = "oobag", 
                                                       trace = TRUE), 
                               center = FALSE)



plot(risk(NON_TERMINAL_EVENT_COX), type = "l")
which.min(risk(NON_TERMINAL_EVENT_COX))

COX_NON_TERMINAL_MSTOP_OPT_OOBAG <- which.min(risk(NON_TERMINAL_EVENT_COX))

NON_TERMINAL_EVENT_COX <- glmboost(nonterminal_formula_COX, 
                                   SCR_Nonterminal_Data, 
                               family = CoxPH(), 
                               control = boost_control(mstop = COX_NON_TERMINAL_MSTOP_OPT_OOBAG, 
                                                       nu = boost_nu,
                                                       trace = TRUE), 
                               center = FALSE)


length(coef(NON_TERMINAL_EVENT_COX))

COX_NON_TERMINAL_MSTOP_OPT_OOBAG

names(coef(NON_TERMINAL_EVENT_COX))

mstops_COX_MODELS <- list(TERMINAL = COX_TERMINAL_MSTOP_OPT_OOBAG, 
                          NON_TERMINAL = COX_NON_TERMINAL_MSTOP_OPT_OOBAG)



COX_MODELS_OUTPUT <- list(TERMINAL_mstop = COX_TERMINAL_MSTOP_OPT_OOBAG,
                          NON_TERMINAL_mstop = COX_NON_TERMINAL_MSTOP_OPT_OOBAG,
                          TERMINAL_coef_names = names(coef(TERMINAL_EVENT_COX)), 
                          NON_TERMINAL_coef_names = names(coef(NON_TERMINAL_EVENT_COX))
                          )

#save(COX_MODELS_OUTPUT, file = "Application/SCR_COXMODELS_OUTPUT.RData")


#### Survival functions: 
SURVFIT_THING_MBOOST_M1 <- mboost::survFit(TERMINAL_EVENT_COX)

SURVFIT_THING_MBOOST_M2 <- mboost::survFit(NON_TERMINAL_EVENT_COX)

basehaz(SURVFIT_THING_MBOOST_M1$surv)


km_fit_TERMINAL <- survival::survfit(Surv(terminal_event, terminal_status) ~ 1, data = (SCR_Terminal_Data))
km_fit_NONTERMINAL <- survival::survfit(Surv(nonterminal_event, nonterminal_status) ~ 1, data = (SCR_Nonterminal_Data))

library(bshazard)


TERMINAL_hazardRate <- bshazard(Surv(terminal_event, terminal_status) ~ 1, data = (SCR_Terminal_Data), degree = 1)
NON_TERMINAL_hazardRate <- bshazard(Surv(nonterminal_event, nonterminal_status) ~ 1, data = (SCR_Nonterminal_Data), degree = 1)
plot(TERMINAL_hazardRate)
plot(NON_TERMINAL_hazardRate)

library(muhaz)

plot(kphaz.plot(kphaz.fit(SCR_Terminal_Data$terminal_event, SCR_Terminal_Data$terminal_status)))
plot(kphaz.plot(kphaz.fit(SCR_Nonterminal_Data$nonterminal_event, SCR_Nonterminal_Data$nonterminal_status)))

kphaz.fit(time,status,strata,q=1,method="nelson")

plot(muhaz(times = SCR_Terminal_Data$terminal_event, delta = SCR_Terminal_Data$terminal_status))

###
par(mfrow = c(1, 2))
plot(km_fit_TERMINAL)
plot(SURVFIT_THING_MBOOST_M1)

plot(km_fit_NONTERMINAL)
plot(SURVFIT_THING_MBOOST_M2)


library(tidymodels)
head(tidy(km_fit_TERMINAL, cumhaz = TRUE))


##################
head(km_fit_TERMINAL)


SCR_Terminal_ONLYTRAIN <- SCR_Terminal_Data[which(weights_mstop == 1),]
SCR_NonTerminal_ONLYTRAIN <- SCR_Nonterminal_Data[which(weights_mstop == 1),]


km_fit_TERMINAL <- survival::survfit(Surv(terminal_event, terminal_status) ~ 1, data = SCR_Terminal_Data)
km_fit_NONTERMINAL <- survival::survfit(Surv(nonterminal_event, nonterminal_status) ~ 1, data = SCR_Nonterminal_Data)


km_fit_TERMINAL
km_fit_NONTERMINAL

plot(km_fit_TERMINAL)
plot(km_fit_NONTERMINAL)


plot(tidy(km_fit_TERMINAL)$estimate, type = "l")
head(tidy(km_fit_TERMINAL))

plot(km_fit_NONTERMINAL$surv, type = "l")
plot(km_fit_TERMINAL$surv, type = "l")


plot(km_fit_NONTERMINAL$time, km_fit_NONTERMINAL$surv, type = "step")
plot(km_fit_TERMINAL$time, km_fit_TERMINAL$surv, type = "step")



MUHAZ_TERMINAL <- (muhaz(times = SCR_Terminal_Data$terminal_event, delta = SCR_Terminal_Data$terminal_status))
MUHAZ_NONTERMINAL <- (muhaz(times = SCR_Nonterminal_Data$nonterminal_event, delta = SCR_Nonterminal_Data$nonterminal_status))

length(MUHAZ_TERMINAL$haz.est)
length(MUHAZ_NONTERMINAL$haz.est)

#### Arrange the data from the Kaplan meier estimate: 
KaplanMeierEstimates <- data.frame(Time = c(km_fit_NONTERMINAL$time, km_fit_TERMINAL$time),
                                   Margin = c(rep("Tumour progression", length(km_fit_NONTERMINAL$time)), 
                                              rep("Death", length(km_fit_TERMINAL$time))), 
                                   Hazard = c(rep(0, length(km_fit_NONTERMINAL$time)), 
                                              rep(0, length(km_fit_TERMINAL$time))),
                                   Surv = c(km_fit_NONTERMINAL$surv, km_fit_TERMINAL$surv), 
                                   SurvMin = c(km_fit_NONTERMINAL$lower, km_fit_TERMINAL$lower), 
                                   SurvMax = c(km_fit_NONTERMINAL$upper, km_fit_TERMINAL$upper)
                                   )


MuHazEstimates <- data.frame(Time = c(MUHAZ_NONTERMINAL$est.grid, MUHAZ_TERMINAL$est.grid), 
                             Margin = c(rep("Tumour progression", length(MUHAZ_NONTERMINAL$haz.est)), 
                                        rep("Death", length(MUHAZ_TERMINAL$haz.est))),
                             Hazard = c(MUHAZ_NONTERMINAL$haz.est, MUHAZ_TERMINAL$haz.est)
                             )


DataRugs <- data.frame(EVENTS = SCR_Terminal_Data$terminal_event[which(SCR_Terminal_Data$terminal_status == 1)], 
                       Margin = rep("Death", length(which(SCR_Terminal_Data$terminal_status == 1))))

DataRugs_NONTERMINAL <- data.frame(EVENTS = SCR_Nonterminal_Data$nonterminal_event[which(SCR_Nonterminal_Data$nonterminal_status == 1)], 
                                   Margin = rep("Tumour progression", length(which(SCR_Nonterminal_Data$nonterminal_status == 1))))


ALL_RUGS <- rbind(DataRugs_NONTERMINAL, DataRugs)


MEDIANRUGS <- data.frame(EVENTS = c(570, 1353), 
                         Margin = c("Median", "Median"))


SURVIVALPLOT <- ggplot(MARGINAL_HAZARDS_DATA, aes(Time, Surv, group = Margin, col = Margin)) +
  geom_line(linewidth = 0.75, show.legend = T) +
  geom_step(data = KaplanMeierEstimates, aes(Time, Surv, group = Margin, col = Margin), linetype = "dashed", linewidth = 1, show.legend = T) +
  geom_rug(data = DataRugs, sides = "t", aes(x = EVENTS), inherit.aes = F, alpha = 0.75, col = "#E69F00") + 
  geom_rug(data = DataRugs_NONTERMINAL, sides = "b", aes(x = EVENTS), inherit.aes = F, alpha = 0.75, col = "#4F53B7") + 
  geom_segment(aes(x = 1353, y = 0, xend = 1353, yend = 0.09), col = "red", linewidth = 1) + 
  geom_segment(aes(x = 570, y = 0, xend = 570, yend = 0.09), col = "red", linewidth = 1) + 
  lims(x = c(0, 3000)) + 
  scale_color_manual(values = c("#E69F00","#4F53B7")) + 
  labs(x = "Time (days)", y = "Estimated survival function", col = "", linetype = "Model:") +
  scale_linetype_manual(name = "Model:",
                        values = c(1,2), 
                        labels = c("SurvCopBoost", "Non-parametric (Cox)")) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17)) 



HAZARDPLOT <- ggplot(MARGINAL_HAZARDS_DATA, aes(Time, Hazard, group = Margin, col = Margin)) +
  geom_line(linewidth = 0.75, show.legend = T) +
  geom_line(data = MuHazEstimates, aes(Time, Hazard, group = Margin, col = Margin), linetype = "dashed", linewidth = 1, show.legend = T) +
  geom_rug(data = DataRugs, sides = "t", aes(x = EVENTS), inherit.aes = F, alpha = 0.75, col = "#E69F00") + 
  geom_rug(data = DataRugs_NONTERMINAL, sides = "b", aes(x = EVENTS), inherit.aes = F, alpha = 0.75, col = "#4F53B7") + 
  geom_segment(aes(x = 1353, y = 0, xend = 1353, yend = 0.00015), col = "red", linewidth = 1) + 
  geom_segment(aes(x = 570, y = 0, xend = 570, yend = 0.00015), col = "red", linewidth = 1) + 
  lims(x = c(0, 3000)) + 
  scale_color_manual(values = c("#E69F00","#4F53B7")) + 
  labs(x = "Time (days)", y = "Estimated hazard rate", col = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17)) 



HAZARDPLOTWITHTITLE <- HAZARDPLOT + ggtitle("(a)") + theme(title = element_text(size = 20), 
                                                           plot.title = element_text(hjust = 0.5))

SURVIVALPLOTWITHTITLE <- SURVIVALPLOT + ggtitle("(b)") + theme(title = element_text(size = 20), 
                                                               plot.title = element_text(hjust = 0.5)) + 
  scale_linetype_manual(name = "Model:",
                        values = c(1,2), 
                        labels = c("SurvCopBoost", "Non-parametric (Cox)"))

ggpubr::ggarrange(HAZARDPLOTWITHTITLE, SURVIVALPLOTWITHTITLE, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")





### 
MARGINAL_HAZARDS_DATA_ALLLINES <- rbind(MARGINAL_HAZARDS_DATA, KaplanMeierEstimates[,1:4])

MARGINAL_HAZARDS_DATA_ALLLINES$Estimator <- factor(c(rep("SurvCopBoost", nrow(MARGINAL_HAZARDS_DATA)), 
                                                     rep("Semi-parametric (Cox)", nrow(KaplanMeierEstimates))), 
                                                   levels = c("SurvCopBoost", "Semi-parametric (Cox)"), 
                                                              order = TRUE)

SURVIVALPLOT <- ggplot(MARGINAL_HAZARDS_DATA_ALLLINES, aes(Time, Surv, group = interaction(Margin, Estimator), 
                                                           col = Margin, linetype = Estimator)) +
  geom_line(linewidth = 0.75) +
  geom_rug(data = DataRugs, sides = "t", aes(x = EVENTS), inherit.aes = F, alpha = 0.75, col = "#E69F00") + 
  geom_rug(data = DataRugs_NONTERMINAL, sides = "b", aes(x = EVENTS), inherit.aes = F, alpha = 0.75, col = "#4F53B7") + 
  geom_segment(aes(x = 1353, y = 0, xend = 1353, yend = 0.09), col = "red", linewidth = 1) + 
  geom_segment(aes(x = 570, y = 0, xend = 570, yend = 0.09), col = "red", linewidth = 1) + 
  lims(x = c(0, 3000)) + 
  scale_color_manual(values = c("#E69F00","#4F53B7")) + 
  labs(x = "Time (days)", y = "Estimated survival function", col = "", linetype = "Estimator:") +
  scale_linetype_manual(name = "Model:",
                        values = c(1,2), 
                        labels = c("SurvCopBoost", "Non-parametric (Cox)")) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17),
        legend.key.size = unit(1.5,"cm")) +
  guides(linetype = guide_legend(override.aes = list(color = "black")))


SURVIVALPLOT



MARGINAL_HAZARDS_DATA_ALLLINES_ONLYHAZARD <- rbind(MARGINAL_HAZARDS_DATA[,1:3], MuHazEstimates)


MARGINAL_HAZARDS_DATA_ALLLINES_ONLYHAZARD$Estimator <- factor(c(rep("SurvCopBoost", nrow(MARGINAL_HAZARDS_DATA)), 
                                                     rep("Semi-parametric (Cox)", nrow(MuHazEstimates))), 
                                                   levels = c("SurvCopBoost", "Semi-parametric (Cox)"), 
                                                   order = TRUE)

HAZARDPLOT <- ggplot(MARGINAL_HAZARDS_DATA_ALLLINES_ONLYHAZARD, aes(Time, Hazard, group = interaction(Margin, Estimator), col = Margin, linetype = Estimator)) +
  geom_line(linewidth = 0.75) +
  geom_rug(data = DataRugs, sides = "t", aes(x = EVENTS), inherit.aes = F, alpha = 0.75, col = "#E69F00") + 
  geom_rug(data = DataRugs_NONTERMINAL, sides = "b", aes(x = EVENTS), inherit.aes = F, alpha = 0.75, col = "#4F53B7") + 
  geom_segment(aes(x = 1353, y = 0, xend = 1353, yend = 0.00015), col = "red", linewidth = 1) + 
  geom_segment(aes(x = 570, y = 0, xend = 570, yend = 0.00015), col = "red", linewidth = 1) + 
  lims(x = c(0, 3000)) + 
  scale_color_manual(values = c("#E69F00","#4F53B7")) + 
  labs(x = "Time (days)", y = "Estimated hazard rate", col = "", linetype = "Estimator:") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17),
        legend.key.size = unit(1.5,"cm")) +
  guides(linetype = guide_legend(override.aes = list(color = "black")))


HAZARDPLOT



HAZARDPLOTWITHTITLE <- HAZARDPLOT + ggtitle("(a)") + theme(title = element_text(size = 20), 
                                                           plot.title = element_text(hjust = 0.5))

SURVIVALPLOTWITHTITLE <- SURVIVALPLOT + ggtitle("(b)") + theme(title = element_text(size = 20), 
                                                               plot.title = element_text(hjust = 0.5))


ggpubr::ggarrange(HAZARDPLOTWITHTITLE, SURVIVALPLOTWITHTITLE, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")






