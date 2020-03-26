############################################################################
######### A Meta Analysis on Children’s Trust in Social Robots #############
############################################################################

# title: "A Meta Analysis on Children’s Trust in Social Robots"
# author: 'Rebecca Stower · Natalia Calvo-Barajas · Ginevra Castellano · Arvid Kappas'
# date: "20/03/2020"

###################### Set things up #################################### 
# To run the code in this document, you should make sure these packages are installed. 
# If they are not installed, use install.packages and the package name in ""
require(tidyverse)      ### For shaping data structures and plotting
require(RCurl)          ### For importing data via urls
require(metafor)        ### For meta-analysis
require(pwr)            ### For Power calculations
library(dplyr)
library(Hmisc)
library(clubSandwich)

source("correlation.R")           ### Correlations function
source("effectsizes_variance.R")  ### Calculate effect sizes and variance
source("MRE_model.R")             ### Multivariate Random Effect Analyses
source("visualizations.R")        ### Meta-analytic visualizations
source("heterogeneity.R")         ### Heterogeneity

# Install packages using the pacman package
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, metafor)


###################### Load Data #################################### 
# For this Meta Analysis, we have set four categories for eligibility decision: include, exclude, liking, and incompleted. 
# This part focuses on load the data and divided according with the eligibility_decision for further analysis.  
# We are directly reading our practice spreadsheet into R here
dataset_url = getURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vQeQHaYFl1a8Pm5oz-k2oYyb6IUpJ7NLeSgSo44wWSCsfYbexgxa7i7ZHha5s8wG3jCNr_dwcsoEFut/pub?output=csv")

# Now we can open the spreadsheet directly 
db_total = read.csv(textConnection(dataset_url))
db_total$mean_age_1 <- as.numeric(levels(db_total$mean_age_1))[db_total$mean_age_1]
db_total$avg_length <- as.numeric(levels(db_total$avg_length))[db_total$avg_length]

# Information that meet the inclusion criteria 
db_include <- db_total[which(db_total$eligibility_decision == 'include'),]

# Information related to liking
db_liking <- db_total[which(db_total$eligibility_decision == 'liking'),]

# Information incompleted
db_incomplete <- db_total[which(db_total$eligibility_decision == 'incomplete'),]

# Information excluded
db_exclude <- db_total[which(db_total$eligibility_decision == 'exclude'),]

cat('Total of variables: ', nrow(db_total), "\n")
cat('Variables Included for Trust (Social and Competency): ', nrow(db_include), "\n")
cat('Variables Liking: ', nrow(db_liking), "\n")
cat('Variables Incompleted', nrow(db_incomplete), "\n")
cat('Variables Excluded', nrow(db_exclude), "\n")

###################### Descriptive Statistics #################################### 
# Merge Liking + Social + Competency Trust
db_ct_st_liking <- merge(db_include, db_liking, all=TRUE, sort=TRUE)
cat('Mean Age:', mean(db_ct_st_liking$mean_age_1), '\n')
cat('SD Age:', sd(db_ct_st_liking$mean_age_1), '\n')
cat('media Interaction Length:', median(db_ct_st_liking$avg_length), '\n')
cat('studies included CT+ST+Liking:', length(unique(db_ct_st_liking$short_cite)), '\n')
cat('Media Sample Size:', median(db_ct_st_liking$N), '\n')

describe(db_ct_st_liking$IV_category_MA)
describe(db_ct_st_liking$IV_category_MA, db_ct_st_liking$dependent_measure)
cat('per measure \n')
describe(db_ct_st_liking$dependent_variable)

# Liking + Social + Competency trust only for IV: Robot_related
db_ct_st_liking <- db_ct_st_liking[which(db_ct_st_liking$IV_category_MA == 'robot'),]
cat('Total of measures for Liking + Comptency + Social Trust IV: Robot Related: ', nrow(db_ct_st_liking), "\n")
describe(db_ct_st_liking$IV_category)

##################### Imputed Correlations #######################################
db_corr_general <- compute_correlations(db_ct_st_liking)

##################### Multivariate Random Effects ###############################
### Multivariate Random Effects Models
# We are interested in four models 
# 1. Competency + Social Trust 
# 2. Social Trust + Liking 
# 3. Competency Trust 
# 4. Social Trust 

## Moderator Analysis - Multivariate Mixed-Effects Meta-Regression Model 
# We have 6 Moderators for our MA
# Moderator 1: Age 
# Moderator 2: Robot Type: Humanoid, Zoomorphic 
# Moderator 3: Interaction Type: learning_task, interview, game 
# Morerator 4: Interaction Length 
# Moderator 5: Robot Operation: WoZ, semi-autonomous, autonomous

################ 1. COMPETENCY + SOCIAL TRUST MODEL ############################
"Dependent Variables: competency_trust + social_trust 
Reported outcomes include: 
competency_trust: reliability_task, helpfulness, trust_label, following_intructions, belief_in_robot 
social_trust: self_disclosure, helpfulness(affective), honesty, keep_secrets, follow_advice"

# Competency + Social trust only for IV: Robot_related
db_Trust <- db_corr_general[which(db_corr_general$IV_category_MA == 'robot'),]
db_Trust  <- subset(db_Trust, subset = dependent_variable != 'social_liking')
cat('Total of measures for Competency + Social Trust IV: Robot Related: ', nrow(db_Trust), "\n")

# GT. Compute ES
dbES_Trust <- compute_es_and_esvar(db_Trust)

# Debug GT miss data (This is just for debugging)
db_miss_GT <- debug_incomplete_es(dbES_Trust)

# Outlier Removal and Add mean_age, mean_n, and weight_d
dbES_Trust <- remove_missing_EF(dbES_Trust)

# Variance-Covariance Matrix
dat <- data.frame(study = dbES_Trust$short_cite, 
                  yi = dbES_Trust$d_calc, 
                  vi = dbES_Trust$d_var_calc)

# r = 0.26 The modarate global effect between trust and all the factors influencing HRI
V <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.63)

########## 1.2. Multivariate Random Effects - Meta Analysis ###########
"Several studies have an overlap in the individuals used to compute multiple 
effect size estimates for both competency and social trust. We can not treat 
those studies as independent. We use an unique id per study, but also the number 
of experiments in each study as grouping variable. 
(e.g., Jones et al. (2014) reported three effect sizes for competency_trust) 
inner: outcome. Transformed to a factor. This variable indicates when one study 
share the same sample population for diferent measures and the effect sizes are dependent. 
outer: short_cite. This factor corresponds to the study identification. 
struct = 'CS': This is used to specify the variance structure corresponding to the inner. 
The random effects are allowed to have different variances for each outcome and are allowed to be correlated. 
random: inner|outer, struct='CS' "

cat('Multivariate Random Effects Model for ')
cat('Competency + Social Trust \n')
cat('inner|outer -> random = ~factor(expt_num)| short_cite \n')
cat('struct = CS \n')
cat('dependent variable as moderator \n -------------------------')

compandsoc_RE_Model <- multivariate_COR_RE_model(dbES_Trust, V)
summary(compandsoc_RE_Model, digits =3)

########## 1.3. Calculate Heterogeneity: Competency + Social Trust ########
I2_comandsoc <- calculate_I2(compandsoc_RE_Model, V)
cat('The amount of Heterogeneity I^2 for the CT+ST Model is: ', I2_comandsoc)

############ 1.4. Funel Plot: Competency + Social Trust ###########
# This is a functionality provided by the metafor package.
funnel(compandsoc_RE_Model)

########### 1.5. Forest Plot: Competency + Social Trust  ##############
plot_Forest(compandsoc_RE_Model)

######## 1.6. Moderator Analysis, Competency + Social Trust   #####
###### 1.6.1. Moderator 1: Age. DV: Competency + Social Trust ########
# Mixed effects meta-regression model with competency + social trust and mean age as predictors/covariates. 
AgeModel_comandsoc <- multivariate_COR_moderatorAge_model(dbES_Trust, V)
summary(AgeModel_comandsoc, digits=3)
# Heterogeneity Moderator 1
I2_comandsocAge <- calculate_I2(AgeModel_comandsoc, V)
cat('The amount of Heterogeneity I^2 for the CT+ST Model + Moderator:Age is: ', I2_comandsocAge)

##### 1.6.2. Moderator 2: Robot Type. DV: Competency + Social Trust ####
# Mixed effects meta-regression model with competency + social trust and Robot Type as predictors/covariates. 
RobotType_comandsoc <- multivariate_COR_moderatorRobotType_model(dbES_Trust, V)
summary(RobotType_comandsoc, digits = 3)
# Heterogeneity Moderator 2
I2_comandsocRT <- calculate_I2(RobotType_comandsoc, V)
cat('The amount of Heterogeneity I^2 for the CT+ST Model + Moderator:RobotType is: ', I2_comandsocRT)

####### 1.6.3. Moderator 3: Interaction Type. DV: Competency + Social Trust ######
# Mixed effects meta-regression model with competency + social trust and Interaction Type as predictors/covariates. 
InteractionType_comandsoc <- multivariate_COR_moderatorInteractionType_model(dbES_Trust, V)
summary(InteractionType_comandsoc, digits = 3)
# Heterogeneity Moderator 3
I2_comandsocIT <- calculate_I2(InteractionType_comandsoc, V)
cat('The amount of Heterogeneity I^2 for the CT+ST Model + Moderator:InteractionType is: ', I2_comandsocIT)

####### 1.6.4. Moderator 4: Interaction Length. DV: Competency + Social Trust #######
InteractionLength_comandsoc <- multivariate_COR_moderatorInteractionLength_model(dbES_Trust, V)
summary(InteractionLength_comandsoc, digits = 3)
# Heterogeneity Moderator 4
I2_comandsocIL <- calculate_I2(InteractionLength_comandsoc, V)
cat('The amount of Heterogeneity I^2 for the CT+ST Model + Moderator:InteractionLength is: ', I2_comandsocIL)

####### 1.6.5. Moderator 5: Robot Operation. DV: Competency + Social Trust ######
RobotOperation_comandsoc <- multivariate_COR_moderatorRobotOperation_model(dbES_Trust, V)
summary(RobotOperation_comandsoc, digits = 3)
# Heterogeneity Moderator 5
I2_comandsocRO <- calculate_I2(RobotOperation_comandsoc, V)
cat('The amount of Heterogeneity I^2 for the CT+ST Model + Moderator:RobotOperation is: ', I2_comandsocRO)

### 1.6.6. Moderator 6: Independent Variable Category. DV: Competency + Social Trust
IVCategory_comandsoc <- multivariate_COR_moderatorIV_model(dbES_Trust, V)
summary(IVCategory_comandsoc, digits = 3)
# Heterogeneity Moderator 6
I2_comandsocIVCat <- calculate_I2(IVCategory_comandsoc, V)
cat('The amount of Heterogeneity I^2 for the CT+ST Model + Moderator:IndependentVariableCategory is: ', I2_comandsocIVCat)

####### 1.7. Power Analysis Competency + Social Trust ######
# Hancock et al. (2011). Found a large global effect concerning trust and HRI d=+0.67. We used this effect size of interest.  
# We now extract the effect size from our meta-analytic model computed in Section 1.2
EffectSize = as.numeric(compandsoc_RE_Model$b)  
# Let's again assume we want to conduct a paired t-test with the effect size above, how many children do we need to reach 80% power?
comp_power <- pwr.t.test(d = EffectSize[1], power = 0.8, type = "paired", alternative = "two.sided")
cat("Power Competency Trust Experimental Design: Within \n")
comp_power
plot(comp_power)

comp_power <- pwr.t.test(d = EffectSize[1], power = 0.8, type = "two.sample", alternative = "two.sided")
cat("Power Competency Trust Experimental Design: Between \n")
comp_power
plot(comp_power)

soc_power <- pwr.t.test(d = EffectSize[2], power = 0.8, type = "paired", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Within \n")
soc_power
plot(soc_power)

soc_power <- pwr.t.test(d = EffectSize[2], power = 0.8, type = "two.sample", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Between \n")
soc_power
plot(soc_power)


###### 1.8. Current Power ####
"To calculate the real power we need to classify the papers according to its experimental desing: withing/between 
The average sample size we have for competency and social trust is for experimental desing 'between' is n_1 = 38. 
The average sample size we have for competency and social trust is for experimental desing 'within' is n_1 = 18. 
We compute the current power of the Meta Analysis Model of section 1.2"

# To compute power we need to classify the papers according to its experimental desing: withing/between
db_Trust_within  <- subset(db_Trust, subset = participant_design == 'within_two')
cat('Total of measures for Competency + Social with within as study desing ', nrow(db_Trust_within), "\n")
db_Trust_between  <- subset(db_Trust, subset = participant_design == 'between') 
cat('Total of measures for Competency + Social with between as study desing ', nrow(db_Trust_between), "\n")

comp_power_real <- pwr.t.test(n = 18 , d = EffectSize[1], type = "paired", alternative = "two.sided")
cat("Power Competency Trust Experimental Design: Within \n")
comp_power_real
plot(comp_power_real)

comp_power_real <- pwr.t.test(n = 38 , d = EffectSize[1], type = "two.sample", alternative = "two.sided")
cat("Power Competency Trust Experimental Design: Between \n")
comp_power_real
plot(comp_power_real)

soc_power_real <- pwr.t.test(n = 18, d = EffectSize[2], type = "paired", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Within \n")
soc_power_real
plot(soc_power_real)

soc_power_real <- pwr.t.test(n = 38, d = EffectSize[2], type = "two.sample", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Between\n")
soc_power_real
plot(soc_power_real)


################ 2. SOCIAL TRUST + LIKING MODEL ################################
"Dependent Variables: social_trust + social_liking 
Reported outcomes include: 
social_trust: self_disclosure, helpfulness(affective), honesty, keep_secrets, follow_advice 
social_liking: social_companion, liking, rapport, being_friends"

# Liking + Social trust only for IV: Robot_related
db_LikingTrust <- db_corr_general[which(db_corr_general$IV_category_MA == 'robot'),]
db_LikingTrust  <- subset(db_LikingTrust, subset = dependent_variable != 'competency_trust')
cat('Total of measures for Liking + Social Trust IV: Robot Related: ', nrow(db_LikingTrust), "\n")

# GT. Compute ES
dbES_LikingTrust <- compute_es_and_esvar(db_LikingTrust)

# Debug GT miss data (This is just for debugging)
db_miss_GT <- debug_incomplete_es(dbES_LikingTrust)

# Outlier Removal and Add mean_age, mean_n, and weight_d
dbES_LikingTrust <- remove_missing_EF(dbES_LikingTrust)

# Variance-Covariance Matrix
dat <- data.frame(study = dbES_LikingTrust$short_cite, 
                  yi = dbES_LikingTrust$d_calc, 
                  vi = dbES_LikingTrust$d_var_calc)

# r = 0.41 The modarate global effect between trust and all the factors influencing HRI
V_Liking <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.41)

######## 2.2. Multivariate Random Effects - Meta Analysis ##########
cat('Hierarchical Random Effects Model \n')
cat('Liking + Social Trust \n')
cat('inner|outer -> random = ~factor(expt_num)| short_cite \n')
cat('struct = CS \n')
cat('dependent variable as moderator \n -------------------------')

socandliking_RE_Model <- multivariate_COR_RE_model(dbES_LikingTrust, V_Liking)
summary(socandliking_RE_Model, digits=3)

####### 2.3. Calculate Heterogeinity: Social Trust + Liking ###############
# This statistic can be thought of as the overall I2 value that indicates how much of 
# the total variance can be attributed to the total amount of heterogeneity.
I2_socliking <- calculate_I2(socandliking_RE_Model,  V_Liking)
cat('The amount of Heterogeneity I^2 for the ST + Liking Model is: ', I2_socliking)

######### 2.4. Funel Plot: Social Trust + Liking #######
# This is a functionality provided by the metafor package.
funnel(socandliking_RE_Model)

######## 2.5. Forest Plot: Social Trust + Liking   #######
plot_Forest(socandliking_RE_Model)

## 2.6. Moderator Analysis, Social Trust + Liking  ######
### 2.6.1. Moderator 1: Age. DV: Social Trust + Liking ########
# Mixed effects meta-regression model with social trust + liking and mean age as predictors/covariates. 
AgeModel_socandLiking <- multivariate_COR_moderatorAge_model(dbES_LikingTrust, V_Liking)
summary(AgeModel_socandLiking, digits=3)
# Heterogeneity Moderator 1
I2_soclikingAge <- calculate_I2(AgeModel_socandLiking,  V_Liking)
cat('The amount of Heterogeneity I^2 for the ST + Liking Model + Moderator:Age is: ', I2_soclikingAge)

### 2.6.2. Moderator 2: Robot Type. DV: Social Trust + Liking #####
# Mixed effects meta-regression model with social trust + liking and Robot Type as predictors/covariates. 
RobotType_socandLiking <- model_x(dbES_LikingTrust, V_Liking)
summary(RobotType_socandLiking, digits = 3)
# Heterogeneity Moderator 2
I2_soclikingRT <- calculate_I2(RobotType_socandLiking,  V_Liking)
cat('The amount of Heterogeneity I^2 for the ST + Liking Model + Moderator:RobotType is: ', I2_soclikingRT)

### 2.6.3. Moderator 3: Interaction Type. DV: Social Trust + Liking #####
# Mixed effects meta-regression model with social trust + liking and Interaction Type as predictors/covariates. 
InteractionType_socandLiking <- multivariate_COR_moderatorInteractionType_model(dbES_LikingTrust, V_Liking)
summary(InteractionType_socandLiking, digits = 3)
# Heterogeneity Moderator 3
I2_soclikingIT <- calculate_I2(InteractionType_socandLiking,  V_Liking)
cat('The amount of Heterogeneity I^2 for the ST + Liking Model + Moderator:InteractionType is: ', I2_soclikingIT)

### 2.6.4. Moderator 4: Interaction Length. DV: Social Trust + Liking ####
InteractionLength_socandLiking <- multivariate_COR_moderatorInteractionLength_model(dbES_LikingTrust, V_Liking)
summary(InteractionLength_socandLiking, digits = 3)
# Heterogeneity Moderator 4
I2_soclikingIL <- calculate_I2(InteractionLength_socandLiking,  V_Liking)
cat('The amount of Heterogeneity I^2 for the ST + Liking Model + Moderator:InteractionLength is: ', I2_soclikingIL)

### 2.6.5. Moderator 5: Robot Operation. DV: Social Trust + Liking #####
RobotOperation_socandLiking <- multivariate_COR_moderatorRobotOperation_model(dbES_LikingTrust, V_Liking)
summary(RobotOperation_socandLiking, digits = 3)
# Heterogeneity Moderator 5
I2_soclikingRO <- calculate_I2(RobotOperation_socandLiking,  V_Liking) 
cat('The amount of Heterogeneity I^2 for the ST + Liking Model + Moderator:RobotOperation is: ', I2_soclikingRO)

### 2.6.6. Moderator 6: Independent Variable Category. DV: Social Trust + Liking #####
IVCategory_socandliking <- multivariate_COR_moderatorIV_model(dbES_LikingTrust, V_Liking)
summary(IVCategory_socandliking, digits = 3)
# Heterogeneity Moderator 6
I2_soclikingIVCat <- calculate_I2(IVCategory_socandliking,  V_Liking)
cat('The amount of Heterogeneity I^2 for the ST + Liking Model + Moderator:IndependentVariable_Category is: ', I2_soclikingIVCat)

## 2.7. Power Analysis Social Trust + Liking #####
# We now extract the effect size from our meta-analytic model computed in Section 2.2
EffectSize = as.numeric(socandliking_RE_Model$b)  
# how many children do we need to reach 80% power?
comp_power <- pwr.t.test(d = EffectSize[1], power = 0.8, type = "paired", alternative = "two.sided")
cat("Power Liking Experimental Design: Within \n")
comp_power
plot(comp_power)

comp_power <- pwr.t.test(d = EffectSize[1], power = 0.8, type = "two.sample", alternative = "two.sided")
cat("Power Liking Experimental Design: Between \n")
comp_power
plot(comp_power)

soc_power <- pwr.t.test(d = EffectSize[2], power = 0.8, type = "paired", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Within \n")
soc_power
plot(soc_power)

soc_power <- pwr.t.test(d = EffectSize[2], power = 0.8, type = "two.sample", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Between \n")
soc_power
plot(soc_power)

###### 2.8. Current Power ######
"To calculate the real power we need to classify the papers according to its experimental desing: withing/between 
The average sample size we have for liking and social trust is for experimental desing 'between' is n_1 = 27. 
The average sample size we have for liking and social trust is for experimental desing 'within' is n_1 = 23. 
We compute the current power of the Meta Analysis Model of section 2.2"
db_LikingTrust_within  <- subset(db_LikingTrust, subset = participant_design == 'within_two')
cat('Total of measures for Liking + Social with within as study desing ', nrow(db_LikingTrust_within), "\n")
db_LikingTrust_between  <- subset(db_LikingTrust, subset = participant_design == 'between') 
cat('Total of measures for Liking + Social with between as study desing ', nrow(db_LikingTrust_between), "\n")

comp_power_real <- pwr.t.test(n = 23 , d = EffectSize[1], type = "paired", alternative = "two.sided")
cat("Power Liking Experimental Design: Within \n")
comp_power_real
plot(comp_power_real)

comp_power_real <- pwr.t.test(n = 27 , d = EffectSize[1], type = "two.sample", alternative = "two.sided")
cat("Power Liking Experimental Design: Between \n")
comp_power_real
plot(comp_power_real)

soc_power_real <- pwr.t.test(n = 23, d = EffectSize[2], type = "paired", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Within \n")
soc_power_real
plot(soc_power_real)

soc_power_real <- pwr.t.test(n = 27, d = EffectSize[2], type = "two.sample", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Between\n")
soc_power_real
plot(soc_power_real)

########## 3. COMPETENCY TRUST MODEL ##############
"Dependent Variables: competency_trust 
Reported outcomes include: 
competency_trust: reliability_task, helpfulness, trust_label, following_intructions, belief_in_robot, seeking_information"

# We used the correlation for robot related characterisctics reported in Hancock (2011) r=0.24
# Available data for Competency Trust, only for IV: Robot_related
db_Trust <- db_corr_general[which(db_corr_general$IV_category_MA == 'robot'),]
db_CompetencyTrust  <- subset(db_Trust, subset = dependent_variable == 'competency_trust')
cat('Total of measures for Competency Trust: ', nrow(db_CompetencyTrust), "\n")

# CT. Compute Correlations
dbCorr_CompetencyTrust <- compute_correlations(db_CompetencyTrust)

# CT. Compute ES
dbES_CompetencyTrust <- compute_es_and_esvar(dbCorr_CompetencyTrust)

# Outlier Removal and Add mean_age, mean_n, and weight_d
dbES_CompetencyTrust <- remove_missing_EF(dbES_CompetencyTrust)

dat <- data.frame(study = dbES_CompetencyTrust$short_cite, 
                  yi = dbES_CompetencyTrust$d_calc, 
                  vi = dbES_CompetencyTrust$d_var_calc)

# r = 0.24 The modarate for robot-related characteristics associated with trust in  HRI
V_CT <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.24)

######## 3.2. Multivariate Random Effects - Meta Analysis #########
cat('Hierarchical Random Effects Model for one DV \n')
cat('Competency Trust \n')
cat('inner|outer -> random = ~factor(expt_num)| short_cite \n')
cat('struct = CS \n')

com_RE_Model <- multivariate_uniqueDV_RE_model(dbES_CompetencyTrust, V_CT)
summary(com_RE_Model, digits=3)

########## 3.3. Calculate Heterogeinity: Competency Trust ########
#  This statistic can be thought of as the overall I2 value that indicates how much of the total 
# variance can be attributed to the total amount of heterogeneity.
I2_com <- calculate_I2(com_RE_Model,  V_CT)
cat('The amount of Heterogeneity I^2 for the CT is: ', I2_com)

########## 3.4. Funel Plot: Competency Trust  #######
# This is a functionality provided by the metafor package.
funnel(com_RE_Model)

########## 3.5. Forest Plot: Competency Trust #########
plot_Forest(com_RE_Model)


######## 3.6. Moderator Analysis, Competency Trust  #########
######## 3.6.1. Moderator 1: Age. DV: Competency Trust #########
# Mixed effects meta-regression model with competency trust and mean age as predictors/covariates. 
AgeModel_comp <- multivariate_COR_moderatorAge_model(dbES_CompetencyTrust, V_CT)
summary(AgeModel_comp, digits=3)
# Heterogeneity Moderator 1
I2_comAge <- calculate_I2(AgeModel_comp, V_CT)
cat('The amount of Heterogeneity I^2 for the CT Model + Moderator:Age is: ', I2_comAge)

######## 3.6.2. Moderator 2: Robot Type. DV: Competency Trust #######
# Mixed effects meta-regression model with competency trust and Robot Type as predictors/covariates. 
RobotType_com <- multivariate_COR_moderatorRobotType_model(dbES_CompetencyTrust, V_CT)
summary(RobotType_com, digits = 3)
# Heterogeneity Moderator 2
I2_comRT <- calculate_I2(RobotType_com, V_CT)
cat('The amount of Heterogeneity I^2 for the CT Model + Moderator:RobotType is: ', I2_comRT)

### 3.6.3. Moderator 3: Interaction Type. DV: Competency Trust ####
# Mixed effects meta-regression model with competency trust and Interaction Type as predictors/covariates. 
InteractionType_com <- multivariate_COR_moderatorInteractionType_model(dbES_CompetencyTrust, V_CT)
summary(InteractionType_com, digits = 3)
# Heterogeneity Moderator 3
I2_comIT <- calculate_I2(InteractionType_com, V_CT)
cat('The amount of Heterogeneity I^2 for the CT Model + Moderator:InteractionType is: ', I2_comIT)

### 3.6.4. Moderator 4: Interaction Length. DV: Competency Trust #####
InteractionLength_com <- multivariate_COR_moderatorInteractionLength_model(dbES_CompetencyTrust, V_CT)
summary(InteractionLength_com, digits = 3)
# Heterogeneity Moderator 4
I2_comIL <- calculate_I2(InteractionLength_com, V_CT)
cat('The amount of Heterogeneity I^2 for the CT Model + Moderator:InteractionLength is: ', I2_comIL)

### 3.6.5. Moderator 5: Robot Operation. DV: Competency Trust #####
RobotOperation_com <- multivariate_COR_moderatorRobotOperation_model(dbES_CompetencyTrust, V_CT)
summary(RobotOperation_com, digits = 3)
# Heterogeneity Moderator 5
I2_comRO <- calculate_I2(RobotOperation_com, V_CT)
cat('The amount of Heterogeneity I^2 for the CT Model + Moderator:RobotOperation is: ', I2_comRO)

### 3.6.6. Moderator 6: Independent Variable Category. DV: Competency Trust #####
IVCategory_CT <- multivariate_COR_moderatorIV_model(dbES_CompetencyTrust, V_CT)
summary(IVCategory_CT, digits = 3)
# Heterogeneity Moderator 6
I2_comIVCat <- calculate_I2(IVCategory_CT, V_CT)

######### 3.7. Power Analysis Competency Trust #####
# We now extract the effect size from our meta-analytic model computed in Section 3.2
EffectSize = as.numeric(com_RE_Model$b)  
# How many children do we need to reach 80% power?
comp_power <- pwr.t.test(d = EffectSize, power = 0.8, type = "paired", alternative = "two.sided")
cat("Power Liking Experimental Design: Within \n")
comp_power
plot(comp_power)

comp_power <- pwr.t.test(d = EffectSize, power = 0.8, type = "two.sample", alternative = "two.sided")
cat("Power Liking Experimental Design: Between \n")
comp_power
plot(comp_power)

########### 3.8. Current Power ##########
"To calculate the real power we need to classify the papers according to its experimental desing: withing/between 
The average sample size we have for competency trust is for experimental desing 'between' is n_1 = 40. 
The average sample size we have for competency trust is for experimental desing 'within' is n_1 = 17. 
We compute the current power of the Meta Analysis Model of section 3.2"
db_CompetencyTrust_within  <- subset(db_CompetencyTrust, subset = participant_design == 'within_two')
cat('Total of measures for Competency with within as study desing ', nrow(db_CompetencyTrust_within), "\n")
db_CompetencyTrust_between  <- subset(db_CompetencyTrust, subset = participant_design == 'between') 
cat('Total of measures for Competency with between as study desing ', nrow(db_CompetencyTrust_between), "\n")

comp_power_real <- pwr.t.test(n = 17 , d = EffectSize, type = "paired", alternative = "two.sided")
cat("Power Competency Trust Experimental Design: Within \n")
comp_power_real
plot(comp_power_real)

comp_power_real <- pwr.t.test(n = 40 , d = EffectSize, type = "two.sample", alternative = "two.sided")
cat("Power Competency Trust Experimental Design: Between \n")
comp_power_real
plot(comp_power_real)
cat('The amount of Heterogeneity I^2 for the CT Model + Moderator:IndependentVariable_Category is: ', I2_comIVCat)


########### 4. SOCIAL TRUST MODEL ###########
"Dependent variable: Social trust. \
Dependent Measure: self disclosure, helpfulness, honesty, keep secrets, follow suggestions"

# We used the correlation for robot related characterisctics reported in Hancock (2011) r=0.24
# Available data for Social Trust
db_Trust <- db_corr_general[which(db_corr_general$IV_category_MA == 'robot'),]
db_SocialTrust  <- subset(db_Trust, subset = dependent_variable == 'social_trust')
cat('Total of measures for Social Trust: ', nrow(db_SocialTrust), "\n")

# ST. Compute Correlations
dbCorr_SocialTrust <- compute_correlations(db_SocialTrust)

# ST. Compute ES
dbES_SocialTrust <- compute_es_and_esvar(dbCorr_SocialTrust)

# Outlier Removal and Add mean_age, mean_n, and weight_d
dbES_SocialTrust <- remove_missing_EF(dbES_SocialTrust)
dat <- data.frame(study = dbES_SocialTrust$short_cite, 
                  yi = dbES_SocialTrust$d_calc, 
                  vi = dbES_SocialTrust$d_var_calc)

# r = 0.24 The modarate for robot-related characteristics associated with trust in  HRI
V_ST <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.24)

###### 4.2 Multivariate Random Effects - Meta Analysis ######
cat('Hierarchical Random Effects Model for one DV \n')
cat('Social Trust \n')
cat('inner|outer -> random = ~factor(expt_num)| short_cite \n')
cat('struct = CS \n')
soc_RE_Model <- multivariate_uniqueDV_RE_model(dbES_SocialTrust, V_ST)
print(soc_RE_Model, digits=3)

####### 4.3. Calculate Heterogeinity: Social Trust ########
#  This statistic can be thought of as the overall I2 value that indicates how much 
# of the total variance can be attributed to the total amount of heterogeneity.
I2_soc <- calculate_I2(soc_RE_Model,  V_ST)
cat('The amount of Heterogeneity I^2 for the ST Model is: ', I2_soc)

####### 4.4. Funel Plot: Social Trust ######
# This is a functionality provided by the metafor package.
funnel(soc_RE_Model)

####### 4.5. Forest Plot: Social Trust ######
plot_Forest(soc_RE_Model)

########## 4.6. Moderator Analysis, Social Trust ########
######## 4.6.1. Moderator 1: Age. DV: Social Trust #######
# Mixed effects meta-regression model with social trust and mean age as predictors/covariates. 
AgeModel_soc <- multivariate_COR_moderatorAge_model(dbES_SocialTrust, V_ST)
summary(AgeModel_soc, digits=3)
# Heterogeneity Moderator 1
I2_socAge <- calculate_I2(AgeModel_soc, V_ST)
cat('The amount of Heterogeneity I^2 for the ST Model + Moderator:Age is: ', I2_socAge)

### 4.6.2. Moderator 2: Robot Type. DV: Social Trust #####
# Mixed effects meta-regression model with social trust and Robot Type as predictors/covariates. 
RobotType_soc <- multivariate_COR_moderatorRobotType_model(dbES_SocialTrust, V_ST)
summary(RobotType_soc, digits = 3)
# Heterogeneity Moderator 2
I2_socRT <- calculate_I2(RobotType_soc, V_ST)
cat('The amount of Heterogeneity I^2 for the ST Model + Moderator:RobotType is: ', I2_socRT)

###### 4.6.3. Moderator 3: Interaction Type. DV: Social Trust ######
# Mixed effects meta-regression model with social trust  and Interaction Type as predictors/covariates. 
InteractionType_soc <- multivariate_COR_moderatorInteractionType_model(dbES_SocialTrust, V_ST)
summary(InteractionType_soc, digits = 3)
# Heterogeneity Moderator 3
I2_socIT <- calculate_I2(InteractionType_soc, V_ST)
cat('The amount of Heterogeneity I^2 for the ST Model + Moderator:InteractionType is: ', I2_socIT)

### 4.6.4. Moderator 4: Interaction Length. DV: Social Trust ######
InteractionLength_soc <- multivariate_COR_moderatorInteractionLength_model(dbES_SocialTrust, V_ST)
summary(InteractionLength_soc, digits = 3)
# Heterogeneity Moderator 4
I2_socIL <- calculate_I2(InteractionLength_soc, V_ST)
cat('The amount of Heterogeneity I^2 for the ST Model + Moderator:InteractionLength is: ', I2_socIL)

### 4.6.5. Moderator 5: Robot Operation. DV: Social Trust ######
RobotOperation_soc <- multivariate_COR_moderatorRobotOperation_model(dbES_SocialTrust, V_ST)
summary(RobotOperation_soc, digits = 3)
# Heterogeneity Moderator 5
I2_socRO <- calculate_I2(RobotOperation_soc, V_ST)
cat('The amount of Heterogeneity I^2 for the ST Model + Moderator:RobotOperation is: ', I2_socRO)

### 4.6.6. Moderator 6: Independent Variable Category. DV: Social Trust #####
IVCategory_ST <- multivariate_COR_moderatorIV_model(dbES_SocialTrust, V_ST)
summary(IVCategory_ST, digits = 3)
# Heterogeneity Moderator 6
I2_socIVCat <- calculate_I2(IVCategory_ST , V_ST)
cat('The amount of Heterogeneity I^2 for the ST Model + Moderator:IndenpendentVariable_Category is: ', I2_socIVCat)

### 4.7. Power Analysis Social Trust ####
# We now extract the effect size from our meta-analytic model computed in Section 2.2
EffectSize = as.numeric(soc_RE_Model$b)  
# How many children do we need to reach 80% power?
comp_power <- pwr.t.test(d = EffectSize[1], power = 0.8, type = "paired", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Within \n")
comp_power
plot(comp_power)

comp_power <- pwr.t.test(d = EffectSize[1], power = 0.8, type = "two.sample", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Between \n")
comp_power
plot(comp_power)

######### 4.8. Current Power Social Trust ########
"To calculate the real power we need to classify the papers according to its experimental desing: withing/between 
The average sample size we have for social trust is for experimental desing 'between' is n_1 = 34. 
The average sample size we have for social trust is for experimental desing 'within' is n_1 = 19. 
We compute the current power of the Meta Analysis Model of section 2.2"

db_SocialTrust_within  <- subset(db_SocialTrust, subset = participant_design == 'within_two')
cat('Total of measures for Social with within as study desing ', nrow(db_SocialTrust_within), "\n")
db_SocialTrust_between  <- subset(db_SocialTrust, subset = participant_design == 'between') 
cat('Total of measures for Social with between as study desing ', nrow(db_SocialTrust_between), "\n")

comp_power_real <- pwr.t.test(n = 19 , d = EffectSize, type = "paired", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Within \n")
comp_power_real
plot(comp_power_real)

comp_power_real <- pwr.t.test(n = 34 , d = EffectSize, type = "two.sample", alternative = "two.sided")
cat("Power Social Trust Experimental Design: Between \n")
comp_power_real
plot(comp_power_real)
