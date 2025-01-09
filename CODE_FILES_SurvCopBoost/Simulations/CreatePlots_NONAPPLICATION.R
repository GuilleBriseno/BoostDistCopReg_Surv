

p <- 10

beta11     <- c( 2, 0, 0, 0, 0, 0, rep(0, p-6)) 

beta12    <- c( 0, +1, 0, 1.5, 0, 0, rep(0, p-6)) 

beta21    <- c( 1, 1.5, 0, 0, 0, 0, rep(0, p-6)) 

beta22    <- c( 0, 0.75, 0, +0.75, 0, 0, rep(0, p-6)) 

betarho   <- c( 0, -2, 0, -2, +0, 0, rep(0, p-6)) 



n <- 700

seed <- 1
set.seed(seed)

#### sample design matrix for train:
BigXMatrix_Train <- matrix(runif(n * p, 0, 1), nrow = n, ncol = p, byrow = T)
x.train <- BigXMatrix_Train

#### sample design matrix for test / evaluation:

colnames(x.train) <- paste0("X", 1:p)


########################### 30% censoring
eta_11 <- x.train %*% beta11
eta_12 <- x.train %*% beta12

eta_21 <- x.train %*% beta21
eta_21 <- -0.6 + eta_21

eta_22 <- x.train %*% beta22
eta_22 <- -1 + eta_22


### for the cutting / right-censoring of uniform times:
timecutM1 <- 9.5
timecutM2 <- 8.5


################################ DEPENDENCE: 
eta_rho <- x.train %*% betarho
eta_rho <- 3 + eta_rho


################################ distribution parameters: 
mu1     <- exp(eta_11)
sigma1  <- exp(eta_12)

mu2     <- (eta_21) # identity link
sigma2  <- exp(eta_22)

rho     <- exp(eta_rho)






u1u2 <- VineCopula::BiCopSim(n, family = 3, par = rho)

# transform into non-uniform random variables:
y1 <- qweibull(p = u1u2[,1], scale = mu1, shape = sigma1, lower.tail = FALSE)
# range(y1)
# plot(y1)

y2 <- qlnorm(p = u1u2[,2], meanlog = mu2, sdlog = sigma2, lower.tail = FALSE)
# range(y2)
# plot(y2)

plot(y1, y2)

# Adding random, non-informative censoring:
censoring_time1 <- runif(length(mu1), min = 0, max = timecutM1)
censoring_time2 <- runif(length(mu1), min = 0, max = timecutM2)

status_m1 <- 1 * (y1 <= censoring_time1)
status_m2 <- 1 * (y2 <= censoring_time2)

# # censoring rate margin 1: (around 20% sounds good for now)
table(status_m1)/n
# 
# # censoring rate margin 2: (around 20% sounds good for now)
table(status_m2)/n

# Replace TRUE event times with observed times: 
y1 <- pmin(y1, censoring_time1)
y2 <- pmin(y2, censoring_time2)

plot(y1, y2)

dat_BTE <- data.frame(y1, y2, status_m1, status_m2)


JOINT_CENS_BTE <- c()


for(i in 1:nrow(dat_BTE)){
  
  if(dat_BTE$status_m1[i] == 1 & dat_BTE$status_m2[i] == 1){
    
    JOINT_CENS_BTE[i] <- 1
    
  }
  
  if(dat_BTE$status_m1[i] == 1 & dat_BTE$status_m2[i] == 0){
    
    JOINT_CENS_BTE[i] <- 2
    
  }
  
  if(dat_BTE$status_m1[i] == 0 & dat_BTE$status_m2[i] == 1){
    
    JOINT_CENS_BTE[i] <- 3
    
  }
  
  if(dat_BTE$status_m1[i] == 0 & dat_BTE$status_m2[i] == 0){
    
    JOINT_CENS_BTE[i] <- 4
    
  }
  
}


dat_BTE$JointCens <- factor(JOINT_CENS_BTE, levels=c(1,2,3,4), labels = c("Both event", 
                                                                          "Y1 censored",
                                                                          "Y2 censored", 
                                                                          "Both censored"))

table(JOINT_CENS_BTE)

not_all_of_them <- floor(length(which(JOINT_CENS_BTE==4)) * 0.8)

BTE_CENSTIME <- runif(length(JOINT_CENS_BTE), min = 0, max = 8)

not_all_of_them

change_these <- sample(which(JOINT_CENS_BTE==4), not_all_of_them, replace=FALSE)

dat_BTE$y1[change_these] <- BTE_CENSTIME[change_these]
dat_BTE$y2[change_these] <- BTE_CENSTIME[change_these]





############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

beta11     <- c( 2, 0, 0, 0, 0, 0, rep(0, p-6)) 

beta12    <- c( 0, +1, 0, 1.5, 0, 0, rep(0, p-6)) 

beta21    <- c( 1, 1.5, 0, 0, 0, 0, rep(0, p-6)) 

beta22    <- c( 0, 0.75, 0, +0.75, 0, 0, rep(0, p-6)) 

betarho   <- c( 0, -2, 0, -2, +0, 0, rep(0, p-6)) 



n <- 700

seed <- 10
set.seed(seed)

#### sample design matrix for train:
BigXMatrix_Train <- matrix(runif(n * p, 0, 1), nrow = n, ncol = p, byrow = T)
x.train <- BigXMatrix_Train

#### sample design matrix for test / evaluation:

colnames(x.train) <- paste0("X", 1:p)


########################### 30% censoring
eta_11 <- 0 + x.train %*% beta11
eta_12 <- 0 + x.train %*% beta12

eta_21 <- x.train %*% beta21
eta_21 <- 0 + eta_21

eta_22 <- x.train %*% beta22
eta_22 <- 3 + eta_22


### for the cutting / right-censoring of uniform times:
timecutM1 <- 9.5
timecutM2 <- 8.5


################################ DEPENDENCE: 
eta_rho <- x.train %*% betarho
eta_rho <- 3 + eta_rho


################################ distribution parameters: 
mu1     <- exp(eta_11)
sigma1  <- exp(eta_12)

mu2     <- exp(eta_21) # identity link
sigma2  <- exp(eta_22)

rho     <- exp(eta_rho)



U_sim <- runif(n, 0, 1)
V_sim <- VineCopula::BiCopCondSim(N = n, cond.val = U_sim, cond.var = 1, family = 3, par = rho)

#u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)

# transform into non-uniform random variables:
y1 <- qweibull(p = U_sim, scale = mu1, shape = sigma1, lower.tail = FALSE)
# range(y1)
# plot(y1)

y2 <- qfisk(p = V_sim, scale = mu2, shape1.a = sigma2, lower.tail = FALSE)
#y2 <- qweibull(p = V_sim, scale = mu2, shape = sigma2, lower.tail = FALSE)
# range(y2)
# plot(y2)

plot(y1, y2)
abline(a = 0, b = 1, col = "red")
# Adding random, non-informative censoring:
censoring_time <- runif(length(mu1), min = 0, max = timecutM1)

### All three times: 
all_three_times <- data.frame(y1 = y1, y2 = y2, cens = censoring_time)

only_two_times <- data.frame(y2 = y2, cens = censoring_time)

####### Determine censoring status of non-terminal event: 
minimum_value <- c()

margin1 <- c()
margin2 <- c()

status_m1 <- c()
status_m2 <- c()

# for(i in 1:length(mu1)){
# 
#   minimum_value[i] <- which.min(all_three_times[i,])
# 
#   ####### Determine NON-TERMINAL EVENT: min(T1, T2, C)
#   margin1[i] <- min(all_three_times[i,])
# 
#   status_m1[i] <- ifelse(minimum_value[i] > 1, 0, 1)
# 
#   ####### Determine TERMINAL EVENT
#   margin2[i] <- min(only_two_times[i,])
# 
#   status_m2[i] <- ifelse(minimum_value[i] > 2, 0, 1)
# 
# }

minimum_of_all_three <- apply(all_three_times, 1, which.min)
minimum_of_remaining_two <- apply(only_two_times, 1, which.min)


#status_m1 <- ifelse(minimum_of_all_three == 1, 1, 0)
#status_m2 <- ifelse(minimum_of_all_three == 2, 1, 0)

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



#table(minimum_value)
#sum(table(minimum_value))

# # censoring rate margin 1: (around 20% sounds good for now)
table(status_m1)/n
# 
# # censoring rate margin 2: (around 20% sounds good for now)
table(status_m2)/n

# Replace TRUE event times with observed times: 
y1 <- margin1
y2 <- margin2

plot(y1, y2)

dat_SCR <- data.frame(y1, y2, status_m1, status_m2)

JOINT_CENS_SCR <- c()
 

for(i in 1:nrow(dat_SCR)){
  
  if(dat_SCR$status_m1[i] == 1 & dat_SCR$status_m2[i] == 1){
    
    JOINT_CENS_SCR[i] <- 1
    
  }
  
  if(dat_SCR$status_m1[i] == 1 & dat_SCR$status_m2[i] == 0){
    
    JOINT_CENS_SCR[i] <- 2
    
  }
  
  if(dat_SCR$status_m1[i] == 0 & dat_SCR$status_m2[i] == 1){
    
    JOINT_CENS_SCR[i] <- 3
    
  }
  
  if(dat_SCR$status_m1[i] == 0 & dat_SCR$status_m2[i] == 0){
    
    JOINT_CENS_SCR[i] <- 4
    
  }
  
}


table(JOINT_CENS_SCR)

#JOINT_CENS_SCR[which(JOINT_CENS_SCR==1)] <- sample(c(1,2), length(which(JOINT_CENS_SCR==1)), prob =c(0.7,0.3), replace=TRUE)

#JOINT_CENS_SCR[which(JOINT_CENS_SCR==1)] <- sample(c(1,4), length(which(JOINT_CENS_SCR==1)), prob =c(0.9,0.1), replace=TRUE)

#JOINT_CENS_SCR[which(JOINT_CENS_SCR==1)] <- sample(c(1,3), length(which(JOINT_CENS_SCR==1)), prob =c(0.8,0.2), replace=TRUE)


dat_SCR$JointCens <- factor(JOINT_CENS_SCR, levels=c(1,2,3,4), labels = c("Both event", 
                                                                          "Y2 censored",
                                                                          "Y1 censored", 
                                                                          "Both censored"))



table(JOINT_CENS_SCR)
table(JOINT_CENS_BTE)


#### 
library(ggplot2)


BTE_PLOT <- ggplot(dat_BTE, aes(y1, y2, col = JointCens)) +
  geom_point(size = 1) +
  ggtitle("(a)") + 
  labs(x = "Time to event 1", y = "Time to event 2", col = "Censoring status") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0,8)) + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0,8)) + 
  scale_colour_manual(values = c("#4F53B7","#E69F00", "#11C638", "#F85A3E")) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        #legend.key.size = unit(, "cm"), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size=18),
        plot.title = element_text(hjust = 0.5, size = 20)
        ) +
  guides(color = guide_legend(override.aes = list(size = 4))) 


SCR_PLOT <- ggplot(dat_SCR, aes(y1, y2, col = JointCens)) +
  geom_point(size = 1) +
  ggtitle("(b)") + 
  labs(x = "Time to non-terminal event (Y1)", y = "Time to terminal event (Y2)", col = "Censoring status") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0,8)) + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0,8)) + 
  scale_colour_manual(values = c("#4F53B7","#E69F00", "#11C638", "#F85A3E")) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        #legend.key.size = unit(, "cm"), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size=18),
        plot.title = element_text(hjust = 0.5, size = 20)) +
  guides(color = guide_legend(override.aes = list(size = 4))) 


ggpubr::ggarrange(BTE_PLOT, SCR_PLOT, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")



######### Check out Angela's transformation: 
AC_TransformMatrix <- matrix(c(1, -1, 0, 1), ncol = 2, nrow = 2, byrow = T)

times_vector <- cbind(dat_SCR$y2, dat_SCR$y1)

new_times <- (AC_TransformMatrix) %*% t(times_vector)
plot(t(new_times))



load("Applications/SCR_APPLICATION/SCR_ApplicationDataFiles_NOAGE_Standard.RData")
ALL_DATA <- rbind(SCR_DATA_FITTING, SCR_DATA_EXTERNAL)

head(colnames(ALL_DATA))

plot(ALL_DATA$nonterminal_event, ALL_DATA$terminal_event)

times_vector <- cbind(ALL_DATA$terminal_event, ALL_DATA$nonterminal_event)

new_times <- (AC_TransformMatrix) %*% t(times_vector)

plot(ALL_DATA$nonterminal_event, ALL_DATA$terminal_event)
plot(t(new_times))

par(mfrow = c(2,1))
plot(ALL_DATA$terminal_event, ALL_DATA$nonterminal_event, xlab = "Terminal event", ylab = "Non-terminal event")
plot(t(new_times),  xlab = "Terminal event", ylab = "Non-terminal event")


range(new_times[,-594])

