library(curatedOvarianData)


source(system.file("extdata", "patientselection.config",package="curatedOvarianData"))
ls()

sapply(ls(), function(x) if(!x %in% c("remove.samples", "duplicates")) print(get(x)))

source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))

names(esets) # From here we need to extract FOUR datasets: 

### indices
index_GSE17260 <- 4
index_GSE30161 <- 9
index_GSE9891  <- 13
index_TCGA     <- 15


#  numbers 4, 8, 11 and 14 corresponds to the study IDs from Ganzfried et al. (2013).
# "group=4" corresponds to 110 Japanese patients from the study of Yoshihara et al. (2010) (GEO accession number: GSE17260). 
# Other groupds are the studies of GSE30161 (58 patients), 
# GSE9891 (278 patients), 
#and TCGA (557 patients).


######## Here we do per dataset a cox model using the genetic marker....
SAMPLE1   <- as(phenoData(esets[[index_GSE17260]]), "data.frame")
SAMPLE2   <- as(phenoData(esets[[index_GSE30161]]), "data.frame")
SAMPLE3   <- as(phenoData(esets[[index_GSE9891]]), "data.frame")
SAMPLE4   <- as(phenoData(esets[[index_TCGA]]), "data.frame")




dim(SAMPLE1)
dim(SAMPLE2)
dim(SAMPLE3)
dim(SAMPLE4)

dim(rbind(SAMPLE1, SAMPLE2, SAMPLE3, SAMPLE4)) # n = 912



####################################################################################################### lets try to extract  these variables: 
probeset_GSE_1 <- featureNames(esets[index_GSE17260]) # these are ALL genetic variables' names
probeset_GSE_2 <- featureNames(esets[index_GSE30161])
probeset_GSE_3 <- featureNames(esets[index_GSE9891])
probeset_TCGA <- featureNames(esets[index_TCGA])

length(probeset_GSE_1)
length(probeset_GSE_2)
length(probeset_GSE_3)
length(probeset_TCGA) 

### Checking the dimension of the individual data frames corresponding to the studies: 
dim( t( exprs(esets[[index_GSE17260]])[featureNames(esets[index_GSE17260]), ]) )

dim( t( exprs(esets[[index_GSE30161]])[featureNames(esets[index_GSE30161]), ]) )

dim( t( exprs(esets[[index_GSE9891]])[featureNames(esets[index_GSE9891]), ]) )

dim( t( exprs(esets[[index_TCGA]])[featureNames(esets[index_TCGA]), ]) )



#### Obtain the variables all of these have in common:
INTERSECTION_OF_VARIABLES <- intersect(featureNames( esets[[index_TCGA]] ), featureNames( esets[[index_GSE17260]] ))

INTERSECTION_OF_VARIABLES_2 <- intersect(INTERSECTION_OF_VARIABLES, featureNames( esets[[index_GSE30161]] ))

INTERSECTION_OF_VARIABLES_3 <- intersect(INTERSECTION_OF_VARIABLES_2, featureNames( esets[[index_GSE9891]] ))


length(INTERSECTION_OF_VARIABLES)
length(INTERSECTION_OF_VARIABLES_2)
length(INTERSECTION_OF_VARIABLES_3)

##### Check that the datasets have the same number of columns:
dim( exprs(esets[[index_GSE17260]])[INTERSECTION_OF_VARIABLES, ] )
dim( exprs(esets[[index_GSE30161]])[INTERSECTION_OF_VARIABLES, ] )
dim( exprs(esets[[index_GSE9891]])[INTERSECTION_OF_VARIABLES, ] )
dim( exprs(esets[[index_TCGA]])[INTERSECTION_OF_VARIABLES, ] )

#### Declare new objects
SAMPLE1_GeneticVariables <- t( exprs(esets[[index_GSE17260]])[ INTERSECTION_OF_VARIABLES , ] )

SAMPLE2_GeneticVariables <- t( exprs(esets[[index_GSE30161]])[ INTERSECTION_OF_VARIABLES , ] )

SAMPLE3_GeneticVariables <- t( exprs(esets[[index_GSE9891]])[ INTERSECTION_OF_VARIABLES , ] )

SAMPLE4_GeneticVariables <- t( exprs(esets[[index_TCGA]])[ INTERSECTION_OF_VARIABLES , ] )

##### Some checks
dim(SAMPLE1_GeneticVariables)
dim(SAMPLE2_GeneticVariables)
dim(SAMPLE3_GeneticVariables)
dim(SAMPLE4_GeneticVariables)


identical(colnames(SAMPLE1_GeneticVariables), colnames(SAMPLE2_GeneticVariables))
identical(colnames(SAMPLE1_GeneticVariables), colnames(SAMPLE3_GeneticVariables))
identical(colnames(SAMPLE1_GeneticVariables), colnames(SAMPLE4_GeneticVariables))
identical(colnames(SAMPLE2_GeneticVariables), colnames(SAMPLE3_GeneticVariables))
identical(colnames(SAMPLE2_GeneticVariables), colnames(SAMPLE4_GeneticVariables))
identical(colnames(SAMPLE3_GeneticVariables), colnames(SAMPLE4_GeneticVariables))

##### remove batch effect: 
combine2 <- function(X1, X2){

  fids <- intersect(featureNames(X1), featureNames(X2))

  X1 <- X1[fids, ]
  X2 <- X2[fids, ]

  ExpressionSet(cbind(exprs(X1), exprs(X2)),
                AnnotatedDataFrame(rbind(as(phenoData(X1), "data.frame"),
                                         as(phenoData(X2), "data.frame"))))
}

data(GSE17260_eset)
data(GSE30161_eset)
data(GSE9891_eset)
data(TCGA_eset)

X <- combine2(GSE17260_eset, GSE30161_eset)
X2 <- combine2(X, GSE9891_eset)
X3 <- combine2(X2, TCGA_eset)

dim(X2)
dim(X3)

boxplot(exprs(X3)[,c(100:150, 169:200, 454:500)])



batch <- factor(c(rep("1", ncol(GSE17260_eset)),
                  rep("2", ncol(GSE30161_eset)),
                  rep("3", ncol(GSE9891_eset)),
                  rep("4", ncol(TCGA_eset))))


combat_data <- sva::ComBat(dat = as.matrix(X3), batch=batch, mod=NULL)

ncol(combat_data)

dim(combat_data)
dim(X3)

boxplot(combat_data[,c(100:150, 169:200, 454:500)])



sum(complete.cases( t( combat_data ) ))










#### STANDARDIZE ALL GENETIC COVARIATES: 
manual_standardise <- function(variable_x){
  
  mu_x <- mean(variable_x)
  
  sd_x <- sd(variable_x)
  
  standard_x <- ( variable_x - mu_x ) / sd_x
  
}


SAMPLE1_GeneticVariables <- apply(SAMPLE1_GeneticVariables, MARGIN = 2, FUN = manual_standardise)
SAMPLE2_GeneticVariables <- apply(SAMPLE2_GeneticVariables, MARGIN = 2, FUN = manual_standardise)
SAMPLE3_GeneticVariables <- apply(SAMPLE3_GeneticVariables, MARGIN = 2, FUN = manual_standardise)
SAMPLE4_GeneticVariables <- apply(SAMPLE4_GeneticVariables, MARGIN = 2, FUN = manual_standardise)


sum(colMeans(SAMPLE1_GeneticVariables))
sum(colMeans(SAMPLE2_GeneticVariables))
sum(colMeans(SAMPLE3_GeneticVariables))
sum(colMeans(SAMPLE4_GeneticVariables))

sum(sapply(1:ncol(SAMPLE1_GeneticVariables), function(i) var(SAMPLE1_GeneticVariables[,i]) ) )
sum(sapply(1:ncol(SAMPLE1_GeneticVariables), function(i) var(SAMPLE2_GeneticVariables[,i]) ) )
sum(sapply(1:ncol(SAMPLE1_GeneticVariables), function(i) var(SAMPLE3_GeneticVariables[,i]) ) )
sum(sapply(1:ncol(SAMPLE1_GeneticVariables), function(i) var(SAMPLE4_GeneticVariables[,i]) ) )


########################## Now we create new data.frames that contain terminal, nonterminal events and status, as well as some clinical covariates and the genetic covariates
# nonterminal_event, 
# terminal_event, 
# nonterminal_status,
# terminal_status, 
# tumorstage,
DATA1 <- data.frame(nonterminal_event = SAMPLE1$days_to_tumor_recurrence, 
                    terminal_event = SAMPLE1$days_to_death, 
                    nonterminal_status = SAMPLE1$recurrence_status, 
                    terminal_status = SAMPLE1$vital_status, 
                    tumorstage = SAMPLE1$tumorstage, 
                    #
                    StudyID = factor(rep("1", nrow(SAMPLE1)))
                    ,
                    #
                    #
                    residual_tumor_size = SAMPLE1$debulking
                    #
                    #
                    , 
                    ##
                    ##
                    SAMPLE1_GeneticVariables
)



#### ALL NAs here can be addressed.... regarding the variable nonterminal status
DATA2 <- data.frame(nonterminal_event = SAMPLE2$days_to_tumor_recurrence, 
                    terminal_event = SAMPLE2$days_to_death, 
                    nonterminal_status = SAMPLE2$recurrence_status, 
                    terminal_status = SAMPLE2$vital_status, 
                    tumorstage = SAMPLE2$tumorstage, 
                    #
                    StudyID = factor(rep("2", nrow(SAMPLE2)))
                    ,
                    #
                    #
                    residual_tumor_size = SAMPLE2$debulking
                    #
                    #
                    ,
                    ##
                    ##
                    SAMPLE2_GeneticVariables
)


## for data2: 
# observation 15 the status of nonterminal event is "norecurrence" (terminal censored and non-terminal censored)
# observation 18 the status of nonterminal event is "norecurrence" (terminal observed and non-terminal censored)
# observation 22 the status of nonterminal event is "norecurrence" (terminal observed and non-terminal censored)
# observation 57 the status of nonterminal event is "norecurrence" (terminal censored and non-terminal censored)

"deceased"
"living"
"norecurrence"
"recurrence"

DATA2$nonterminal_status[15] <- "norecurrence"
DATA2$nonterminal_status[18] <- "norecurrence"
DATA2$nonterminal_status[22] <- "norecurrence"
DATA2$nonterminal_status[57] <- "norecurrence"

sum(complete.cases(DATA2[,-c(5,6)]))
sum(complete.cases(DATA2))
nrow(DATA2)


### Here we have 3 NAs in the nonterminal event time... 
DATA3 <- data.frame(nonterminal_event = SAMPLE3$days_to_tumor_recurrence, 
                    terminal_event = SAMPLE3$days_to_death, 
                    nonterminal_status = SAMPLE3$recurrence_status, 
                    terminal_status = SAMPLE3$vital_status, 
                    tumorstage = SAMPLE3$tumorstage, 
                    #
                    StudyID = factor(rep("3", nrow(SAMPLE3)))
                    ,
                    #
                    #
                    residual_tumor_size = SAMPLE3$debulking
                    #
                    #
                    ,
                    ##
                    ##
                    SAMPLE3_GeneticVariables
)


## for DATA3: 
# observation 24 change the status of nonterminal event is "norecurrence" and set time to 810. (terminal censored and non-terminal censored)
# observation 60 the status of nonterminal event is "norecurrence"  and set time to 1260 (terminal censored and non-terminal censored)
# observation 69 the status of nonterminal event is "norecurrence" and set time to 540 (terminal censored and non-terminal censored)


DATA3$nonterminal_status[24] <- "norecurrence"
DATA3$nonterminal_event[24]  <-  DATA3$terminal_event[24]

DATA3$nonterminal_status[60] <- "norecurrence"
DATA3$nonterminal_event[60]  <-  DATA3$terminal_event[60]

DATA3$nonterminal_status[69] <- "norecurrence"
DATA3$nonterminal_event[69]  <-  DATA3$terminal_event[69]


sum(complete.cases(DATA3[,-c(5,6)]))
sum(complete.cases(DATA3)) ### 3 observations missing their ages or tumorstage or something...
nrow(DATA3)




DATA4 <- data.frame(nonterminal_event = SAMPLE4$days_to_tumor_recurrence, 
                    terminal_event = SAMPLE4$days_to_death, 
                    nonterminal_status = SAMPLE4$recurrence_status, 
                    terminal_status = SAMPLE4$vital_status, 
                    tumorstage = SAMPLE4$tumorstage, 
                    #
                    StudyID = factor(rep("4", nrow(SAMPLE4)))
                    ,
                    #
                    #
                    residual_tumor_size = SAMPLE4$debulking
                    #
                    #
                    ,
                    ##
                    ##
                    SAMPLE4_GeneticVariables
)


### There are 29 NAs here, which can be corrected:

## observation 16 is terminal censored and non-terminal censored
## observation 29 is terminal observed and non-terminal censored
## observation 58 is terminal observed and non-terminal censored
## observation 68 is terminal observed and non-terminal censored
## observation 93 is terminal observed and non-terminal censored
## observation 99 is terminal observed and non-terminal censored
## observation 160 is terminal observed and non-terminal censored
## observation 166 is terminal censored and non-terminal censored
## observation 176 is terminal censored and non-terminal censored
## observation 180 is terminal observed and non-terminal censored
## observation 218 is terminal censored and non-terminal censored
## observation 251 is terminal observed and non-terminal censored
## observation 291 is terminal observed and non-terminal censored
## observation 307 is terminal observed and non-terminal censored
## observation 336 is terminal censored and non-terminal censored
## observation 352 is terminal observed and non-terminal censored
## observation 392 is terminal observed and non-terminal censored
## observation 397 is terminal censored and non-terminal censored
## observation 418 is terminal observed and non-terminal censored
## observation 422 is terminal censored and non-terminal censored
## observation 428 is terminal censored and non-terminal censored
## observation 443 is terminal observed and non-terminal censored
## observation 447 is terminal censored and non-terminal censored
## observation 450 is terminal observed and non-terminal censored
## observation 472 is terminal observed and non-terminal censored
## observation 483 is terminal observed and non-terminal censored
## observation 484 is terminal censored and non-terminal censored
## observation 495 is terminal observed and non-terminal censored
## observation 498 is terminal censored and non-terminal censored
DATA4[16, 1:5]

DATA4$nonterminal_status[16] <- "norecurrence"
DATA4$nonterminal_event[16]  <-  DATA4$terminal_event[16]

DATA4$nonterminal_status[29] <- "norecurrence"
DATA4$nonterminal_event[29]  <-  DATA4$terminal_event[29]

DATA4$nonterminal_status[58] <- "norecurrence"
DATA4$nonterminal_event[58]  <-  DATA4$terminal_event[58]

DATA4$nonterminal_status[68] <- "norecurrence"
DATA4$nonterminal_event[68]  <-  DATA4$terminal_event[68]

DATA4$nonterminal_status[93] <- "norecurrence"
DATA4$nonterminal_event[93]  <-  DATA4$terminal_event[93]

DATA4$nonterminal_status[99] <- "norecurrence"
DATA4$nonterminal_event[99]  <-  DATA4$terminal_event[99]

DATA4$nonterminal_status[160] <- "norecurrence"
DATA4$nonterminal_event[160]  <-  DATA4$terminal_event[160]

DATA4$nonterminal_status[166] <- "norecurrence"
DATA4$nonterminal_event[166]  <-  DATA4$terminal_event[166]

DATA4$nonterminal_status[176] <- "norecurrence"
DATA4$nonterminal_event[176]  <-  DATA4$terminal_event[176]

DATA4$nonterminal_status[180] <- "norecurrence"
DATA4$nonterminal_event[180]  <-  DATA4$terminal_event[180]

DATA4$nonterminal_status[218] <- "norecurrence"
DATA4$nonterminal_event[218]  <-  DATA4$terminal_event[218]

DATA4$nonterminal_status[251] <- "norecurrence"
DATA4$nonterminal_event[251]  <-  DATA4$terminal_event[251]

DATA4$nonterminal_status[291] <- "norecurrence"
DATA4$nonterminal_event[291]  <-  DATA4$terminal_event[291]

DATA4$nonterminal_status[307] <- "norecurrence"
DATA4$nonterminal_event[307]  <-  DATA4$terminal_event[307]

DATA4$nonterminal_status[336] <- "norecurrence"
DATA4$nonterminal_event[336]  <-  DATA4$terminal_event[336]

DATA4$nonterminal_status[352] <- "norecurrence"
DATA4$nonterminal_event[352]  <-  DATA4$terminal_event[352]

DATA4$nonterminal_status[392] <- "norecurrence"
DATA4$nonterminal_event[392]  <-  DATA4$terminal_event[392]

DATA4$nonterminal_status[397] <- "norecurrence"
DATA4$nonterminal_event[397]  <-  DATA4$terminal_event[397]

DATA4$nonterminal_status[418] <- "norecurrence"
DATA4$nonterminal_event[418]  <-  DATA4$terminal_event[418]

DATA4$nonterminal_status[422] <- "norecurrence"
DATA4$nonterminal_event[422]  <-  DATA4$terminal_event[422]

DATA4$nonterminal_status[428] <- "norecurrence"
DATA4$nonterminal_event[428]  <-  DATA4$terminal_event[428]

DATA4$nonterminal_status[443] <- "norecurrence"
DATA4$nonterminal_event[443]  <-  DATA4$terminal_event[443]

DATA4$nonterminal_status[447] <- "norecurrence"
DATA4$nonterminal_event[447]  <-  DATA4$terminal_event[447]

DATA4$nonterminal_status[450] <- "norecurrence"
DATA4$nonterminal_event[450]  <-  DATA4$terminal_event[450]

DATA4$nonterminal_status[472] <- "norecurrence"
DATA4$nonterminal_event[472]  <-  DATA4$terminal_event[472]

DATA4$nonterminal_status[483] <- "norecurrence"
DATA4$nonterminal_event[483]  <-  DATA4$terminal_event[483]

DATA4$nonterminal_status[484] <- "norecurrence"
DATA4$nonterminal_event[484]  <-  DATA4$terminal_event[484]

DATA4$nonterminal_status[495] <- "norecurrence"
DATA4$nonterminal_event[495]  <-  DATA4$terminal_event[495]

DATA4$nonterminal_status[498] <- "norecurrence"
DATA4$nonterminal_event[498]  <-  DATA4$terminal_event[498]

sum(complete.cases(DATA4[,-c(5,6)]))
sum(complete.cases(DATA4))
nrow(DATA4)



dim(DATA1[,-c(6)])
dim(DATA1[complete.cases(DATA1[,-c(6)]), -c(6)])


dim(DATA2[,-c(6)])
dim(DATA2[complete.cases(DATA2[,-c(6)]), -c(6)])


dim(DATA3[,-c(6)])
dim(DATA3[complete.cases(DATA3[,-c(6)]), -c(6)])


dim(DATA4[,-c(6)])
dim(DATA4[complete.cases(DATA4[,-c(6)]), -c(6)])


#### 
sum(complete.cases(DATA1[,-c(5,6)])) + 
  sum(complete.cases(DATA2[,-c(5,6)])) + 
  sum(complete.cases(DATA3[,-c(5,6)])) + 
  sum(complete.cases(DATA4[,-c(5,6)])) 

sum(complete.cases(DATA1[,-c(5)])) + 
  sum(complete.cases(DATA2[,-c(5)])) + 
  sum(complete.cases(DATA3[,-c(5)])) + 
  sum(complete.cases(DATA4[,-c(5)])) 

#
sum(complete.cases(DATA1[,-c(6)])) + 
  sum(complete.cases(DATA2[,-c(6)])) + 
  sum(complete.cases(DATA3[,-c(6)])) + 
  sum(complete.cases(DATA4[,-c(6)])) 



SCR_CLINICAL <- rbind(DATA1,
                      DATA2,
                      DATA3,
                      DATA4)

SCR_CLINICAL <- SCR_CLINICAL[complete.cases(SCR_CLINICAL),]


SCR_CLINICAL$nonterminal_status <- ifelse(SCR_CLINICAL$nonterminal_status == "recurrence", 1, 0)
SCR_CLINICAL$terminal_status    <- ifelse(SCR_CLINICAL$terminal_status == "deceased", 1, 0)


##### Create joint censoring variable for visualisation: 
Joint_Censoring <- rep(NA, nrow(SCR_CLINICAL))

for(i in 1:nrow(SCR_CLINICAL)){
  
  if(SCR_CLINICAL$nonterminal_status[i] == 1 & SCR_CLINICAL$terminal_status[i] == 1){
    
    Joint_Censoring[i] <- 1
    
  }
  if(SCR_CLINICAL$nonterminal_status[i] == 1 & SCR_CLINICAL$terminal_status[i] == 0){
    
    Joint_Censoring[i] <- 2
    
  }
  if(SCR_CLINICAL$nonterminal_status[i] == 0 & SCR_CLINICAL$terminal_status[i] == 1){
    
    Joint_Censoring[i] <- 3
    
  }
  if(SCR_CLINICAL$nonterminal_status[i] == 0 & SCR_CLINICAL$terminal_status[i] == 0){
    
    Joint_Censoring[i] <- 4
    
  }
}

Joint_Censoring <- factor(Joint_Censoring, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Both uncensored", 
                                     "Y2 censored", 
                                     "Y1 censored", 
                                     "Both censored"))

SCR_CLINICAL_VIS <- data.frame(nonterminal_event = SCR_CLINICAL$nonterminal_event, 
                               terminal_event = SCR_CLINICAL$terminal_event, 
                               StudyID = SCR_CLINICAL$StudyID, 
                               Censoring = Joint_Censoring)

SCR_CLINICAL_VIS$StudyID <- factor(SCR_CLINICAL_VIS$StudyID, 
                                   levels = c(1, 2, 3, 4), 
                                   labels = c("GSE17260",
                                              "GSE30161",
                                              "GSE9891",
                                              "TCGA"))




#### Add median survival times to both margins: 

######### some visualisation
library(ggplot2)
library(viridis)
ggplot(SCR_CLINICAL_VIS, aes(x = nonterminal_event, y = terminal_event, col = StudyID)) +
  geom_point(aes(shape = Censoring), size = 2, alpha = 0.6) +
  scale_shape_manual(values = c(20, 17, 18, 15)) + 
  labs(x = "Time to tumor progression (days)", y = "Time to death (days)", col = "Study:", shape = "Censoring: ") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 19),
        axis.text = element_text(size = 11), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "cm")
  )

ggplot(SCR_CLINICAL_VIS, aes(x = nonterminal_event, y = terminal_event, col = Censoring)) +
  geom_point(size = 1.75, alpha = 0.6) +
  labs(x = "Time to tumor progression (days)", y = "Time to death (days)", col= "Censoring") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 19),
        axis.text = element_text(size = 11)
  )







SAMPLE1 <- DATA1[complete.cases(DATA1), ]
SAMPLE2 <- DATA2[complete.cases(DATA2), ]
SAMPLE3 <- DATA3[complete.cases(DATA3), ]
SAMPLE4 <- DATA4[complete.cases(DATA4), ]

dim(SAMPLE1)
dim(SAMPLE2)
dim(SAMPLE3)
dim(SAMPLE4)

nrow(SAMPLE1) + nrow(SAMPLE2) + nrow(SAMPLE3) + nrow(SAMPLE4)


################## Sample from the respective studies: # this guarantees that all studies are represented during fitting:
set.seed(1)
weights_stage1_study1 <- sample(c(0,1), size = nrow(SAMPLE1), replace = TRUE, prob = c(0.7, 0.3))
set.seed(1)
weights_stage1_study2 <- sample(c(0,1), size = nrow(SAMPLE2), replace = TRUE, prob = c(0.7, 0.3))
set.seed(1)
weights_stage1_study3 <- sample(c(0,1), size = nrow(SAMPLE3), replace = TRUE, prob = c(0.7, 0.3))
set.seed(1)
weights_stage1_study4 <- sample(c(0,1), size = nrow(SAMPLE4), replace = TRUE, prob = c(0.7, 0.3))


### Those that are = 0 will be used for fitting:
SAMPLE1_EXTERNAL_INDX <- which(weights_stage1_study1 != 0)
SAMPLE2_EXTERNAL_INDX <- which(weights_stage1_study2 != 0)
SAMPLE3_EXTERNAL_INDX <- which(weights_stage1_study3 != 0)
SAMPLE4_EXTERNAL_INDX <- which(weights_stage1_study4 != 0)



sum(c(weights_stage1_study1, 
  weights_stage1_study2, weights_stage1_study3, weights_stage1_study4)) # 


length(SAMPLE1_EXTERNAL_INDX) # 
length(SAMPLE2_EXTERNAL_INDX) # 
length(SAMPLE3_EXTERNAL_INDX) # 
length(SAMPLE4_EXTERNAL_INDX) # 

set.seed(1)
sample1_oobag_indx <- sample(c(0, 1), length(SAMPLE1_EXTERNAL_INDX), replace = TRUE, prob = c(0.5, 0.5))
set.seed(1)
sample2_oobag_indx <- sample(c(0, 1), length(SAMPLE2_EXTERNAL_INDX), replace = TRUE, prob = c(0.5, 0.5))
set.seed(1)
sample3_oobag_indx <- sample(c(0, 1), length(SAMPLE3_EXTERNAL_INDX), replace = TRUE, prob = c(0.5, 0.5))
set.seed(1)
sample4_oobag_indx <- sample(c(0, 1), length(SAMPLE4_EXTERNAL_INDX), replace = TRUE, prob = c(0.5, 0.5))

SAMPLE1_EXTERNAL_INDX_OOBAG <- SAMPLE1_EXTERNAL_INDX[which(sample1_oobag_indx == 1)]
SAMPLE2_EXTERNAL_INDX_OOBAG <- SAMPLE2_EXTERNAL_INDX[which(sample2_oobag_indx == 1)]
SAMPLE3_EXTERNAL_INDX_OOBAG <- SAMPLE3_EXTERNAL_INDX[which(sample3_oobag_indx == 1)]
SAMPLE4_EXTERNAL_INDX_OOBAG <- SAMPLE4_EXTERNAL_INDX[which(sample4_oobag_indx == 1)]

SAMPLE1_EXTERNAL_INDX_LOGSCORE <- SAMPLE1_EXTERNAL_INDX[which(sample1_oobag_indx == 0)]
SAMPLE2_EXTERNAL_INDX_LOGSCORE <- SAMPLE2_EXTERNAL_INDX[which(sample2_oobag_indx == 0)]
SAMPLE3_EXTERNAL_INDX_LOGSCORE <- SAMPLE3_EXTERNAL_INDX[which(sample3_oobag_indx == 0)]
SAMPLE4_EXTERNAL_INDX_LOGSCORE <- SAMPLE4_EXTERNAL_INDX[which(sample4_oobag_indx == 0)]



##################### ONLY OBSERVATIONS WHERE THE INDEX WAS  == 0:
SAMPLE1_FITTING <- SAMPLE1[-SAMPLE1_EXTERNAL_INDX, ]
SAMPLE2_FITTING <- SAMPLE2[-SAMPLE2_EXTERNAL_INDX, ] 
SAMPLE3_FITTING <- SAMPLE3[-SAMPLE3_EXTERNAL_INDX, ]
SAMPLE4_FITTING <- SAMPLE4[-SAMPLE4_EXTERNAL_INDX, ]


######################################################## ONLY OBSERVATIONS FOR EXTERNAL USE
SAMPLE1_EXTERNAL_OOBAG <- SAMPLE1[SAMPLE1_EXTERNAL_INDX_OOBAG, ]
SAMPLE2_EXTERNAL_OOBAG <- SAMPLE2[SAMPLE2_EXTERNAL_INDX_OOBAG, ]
SAMPLE3_EXTERNAL_OOBAG <- SAMPLE3[SAMPLE3_EXTERNAL_INDX_OOBAG, ]
SAMPLE4_EXTERNAL_OOBAG <- SAMPLE4[SAMPLE4_EXTERNAL_INDX_OOBAG, ]

###################################################### ONLY OBSERVATIONS FOR OOBAG RISK:
SAMPLE1_EXTERNAL <- SAMPLE1[SAMPLE1_EXTERNAL_INDX_LOGSCORE, ]
SAMPLE2_EXTERNAL <- SAMPLE2[SAMPLE2_EXTERNAL_INDX_LOGSCORE, ]
SAMPLE3_EXTERNAL <- SAMPLE3[SAMPLE3_EXTERNAL_INDX_LOGSCORE, ]
SAMPLE4_EXTERNAL <- SAMPLE4[SAMPLE4_EXTERNAL_INDX_LOGSCORE, ]


### ALL TOGETHER
SCR_DATA <- rbind(SAMPLE1, 
  SAMPLE2, 
  SAMPLE3, 
  SAMPLE4)


#### final assembly of datasets:
SCR_DATA_FITTING <- rbind(SAMPLE1_FITTING, 
  SAMPLE2_FITTING, 
  SAMPLE3_FITTING, 
  SAMPLE4_FITTING, 
  SAMPLE1_EXTERNAL_OOBAG, 
  SAMPLE2_EXTERNAL_OOBAG, 
  SAMPLE3_EXTERNAL_OOBAG, 
  SAMPLE4_EXTERNAL_OOBAG
)


SCR_DATA_EXTERNAL <- rbind(SAMPLE1_EXTERNAL,
  SAMPLE2_EXTERNAL, 
  SAMPLE3_EXTERNAL, 
  SAMPLE4_EXTERNAL)


dim(SCR_DATA)
dim(SCR_DATA_EXTERNAL)
dim(SCR_DATA_FITTING)


nrow(SCR_DATA_EXTERNAL) + nrow(SCR_DATA_FITTING)
nrow(rbind(SAMPLE1, 
  SAMPLE2, SAMPLE3, SAMPLE4))


#### creation of WEIGHTS for oobag risk during fitting: 
weights_mstop <- c(rep(0, nrow(rbind(SAMPLE1_FITTING, 
  SAMPLE2_FITTING, 
  SAMPLE3_FITTING, 
  SAMPLE4_FITTING))),
  rep(1, nrow(rbind(SAMPLE1_EXTERNAL_OOBAG, 
    SAMPLE2_EXTERNAL_OOBAG, 
    SAMPLE3_EXTERNAL_OOBAG, 
    SAMPLE4_EXTERNAL_OOBAG))))

length(weights_mstop) # should be equal to nrow(SCR_DATA_FITTING)
nrow(SCR_DATA_FITTING)
mean(weights_mstop)


################################################################### Final changes / adjustments to the data: 
### Status variables are turned into numeric variables
### Tumorstage is turned to factor with dummy.


terminal_status_numeric_FITTING     <- ifelse(SCR_DATA_FITTING$terminal_status == "deceased", 1, 0)
nonterminal_status_numeric_FITTING  <- ifelse(SCR_DATA_FITTING$nonterminal_status == "recurrence", 1, 0)

terminal_status_numeric_EXTERNAL     <- ifelse(SCR_DATA_EXTERNAL$terminal_status == "deceased", 1, 0)
nonterminal_status_numeric_EXTERNAL  <- ifelse(SCR_DATA_EXTERNAL$nonterminal_status == "recurrence", 1, 0)

library(fastDummies)

TUMOR_DUMMY_FITTING <- dummy_cols(SCR_DATA_FITTING$tumorstage)
TUMOR_DUMMY_EXTERNAL <- dummy_cols(SCR_DATA_EXTERNAL$tumorstage)

colnames(TUMOR_DUMMY_FITTING)   <- c("Original", "Reference", "TumorStage_2", "TumorStage_3", "TumorStage_4")
colnames(TUMOR_DUMMY_EXTERNAL)  <- c("Original", "TumorStage_2", "TumorStage_3", "TumorStage_4") # in case type 1 is not included in the sample here:


STUDY_DUMMY_FITTING <- dummy_cols(SCR_DATA_FITTING$StudyID)
STUDY_DUMMY_EXTERNAL <- dummy_cols(SCR_DATA_EXTERNAL$StudyID)

colnames(STUDY_DUMMY_FITTING)   <- c("Original", "Reference", "StudyID_2", "StudyID_3", "StudyID_4")
colnames(STUDY_DUMMY_EXTERNAL)  <- c("Original", "Reference", "StudyID_2", "StudyID_3", "StudyID_4")



##### DEBULKING / RESIDUAL TUMOR SIZE

residual_tumor_size_factor_FITTING  <- factor( ifelse(SCR_DATA_FITTING$residual_tumor_size == "optimal", 0, 1) )
residual_tumor_size_factor_EXTERNAL <- factor( ifelse(SCR_DATA_EXTERNAL$residual_tumor_size == "optimal", 0, 1) )



###### Replace variables with their numeric counterparts that we just created:

### Status:
SCR_DATA_FITTING$nonterminal_status   <- nonterminal_status_numeric_FITTING
SCR_DATA_EXTERNAL$nonterminal_status  <- nonterminal_status_numeric_EXTERNAL

SCR_DATA_FITTING$terminal_status      <- terminal_status_numeric_FITTING
SCR_DATA_EXTERNAL$terminal_status     <- terminal_status_numeric_EXTERNAL


### Remove first the factor variables:
SCR_DATA_FITTING_OLD <- SCR_DATA_FITTING
SCR_DATA_EXTERNAL_OLD <- SCR_DATA_EXTERNAL

SCR_DATA_FITTING <- SCR_DATA_FITTING[, -c(5, 7)]
SCR_DATA_EXTERNAL <- SCR_DATA_EXTERNAL[, -c(5, 7)]

### Re-add: Tumor Stage:
SCR_DATA_FITTING$TumorStage_2 <- TUMOR_DUMMY_FITTING$TumorStage_2
SCR_DATA_FITTING$TumorStage_3 <- TUMOR_DUMMY_FITTING$TumorStage_3
SCR_DATA_FITTING$TumorStage_4 <- TUMOR_DUMMY_FITTING$TumorStage_4

SCR_DATA_EXTERNAL$TumorStage_2 <- TUMOR_DUMMY_EXTERNAL$TumorStage_2
SCR_DATA_EXTERNAL$TumorStage_3 <- TUMOR_DUMMY_EXTERNAL$TumorStage_3
SCR_DATA_EXTERNAL$TumorStage_4 <- TUMOR_DUMMY_EXTERNAL$TumorStage_4


#### RE- ADD RESIDUAL TUMOR SIZE: 
SCR_DATA_FITTING$residual_tumor_size  <- residual_tumor_size_factor_FITTING
SCR_DATA_EXTERNAL$residual_tumor_size <- residual_tumor_size_factor_EXTERNAL



#### Add intercept column!
SCR_DATA_EXTERNAL$XInter <- as.numeric( rep(1, nrow(SCR_DATA_EXTERNAL)) )
SCR_DATA_FITTING$XInter <- as.numeric( rep(1, nrow(SCR_DATA_FITTING)) )


nrow(SCR_DATA_EXTERNAL) + nrow(SCR_DATA_FITTING)

nrow(SCR_DATA_FITTING)

length(weights_mstop)

ncol(SCR_DATA_EXTERNAL)
ncol(SCR_DATA_FITTING)


######## These observations are Y = (0, 0), this makes no sense... we remove them altogether...
SCR_DATA_FITTING <- SCR_DATA_FITTING[ -98, ]
SCR_DATA_EXTERNAL <- SCR_DATA_EXTERNAL[ -41, ]


weights_mstop <- weights_mstop[-98]



weights_mstop_original <- weights_mstop
weights_mstop <- ifelse(weights_mstop_original == 1, 0, 1)

table(weights_mstop) / nrow(SCR_DATA_FITTING)



#### SAVE: 
save( 
  SCR_DATA_EXTERNAL, 
  SCR_DATA_FITTING,
  weights_mstop,
  file = "Applications/SCR_APPLICATION/SCR_APPLICATION_DATA_FILES.RData")







