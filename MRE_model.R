# Perform Meta-analysis 
# Now that we have our effect sizes, we can put them into a meta-analytic regression model. 

# Main meta-analysis model 
# We are interested in four models 
# 1. Social Trust 
# 2. Competency Trust 
# 3. Competency + Social Trust 
# 4. Social Trust + Liking 

multivariate_COR_RE_model <- function(db, V){
  "This function computes the multivariate random effect analysis
   The model takes as input the data and the variance-covariance matrix. The model takes the DV 
   as a factor moderator, to investigate the correlation between the two dependent variables."
  
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "REML",         # restricted maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          struct = "CS",
                          mods = ~ dependent_variable - 1)

  return(Multi_REModel)  
}

multivariate_uniqueDV_RE_model <- function(db, V){
  "This function computes the multivariate random effect analysis
   The model takes as input the data and the variance-covariance matrix."
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "REML",         # restricted maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          struct = "CS")
  return(Multi_REModel)  
}

## Moderator Analysis - Multivariate Mixed-Effects Meta-Regression Model 
# We have 6 Moderators for our MA
# Moderator 1: Age 
# Moderator 2: Robot Type: Humanoid, Zoomorphic 
# Moderator 3: Interaction Type: learning_task, interview, game 
# Morerator 4: Interaction Length 
# Moderator 5: Robot Operation: WoZ, semi-autonomous, autonomous
# Moderator 6: Independent Variable Category
# Moderator 7: Type of Measure: Subjective/Objective 

multivariate_COR_moderatorAge_model <- function(db, V){
  "This function computes the multivariate mixed-effect analysis
  Moderator 1: Age"
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "REML",         # restricted maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          struct = "CS",
                          mods = ~ dependent_variable:mean_age -1)
  return(Multi_REModel)  
}

multivariate_COR_moderatorRobotType_model <- function(db, V){
  "This function computes the multivariate mixed-effect analysis
  Moderator 2: Robot Type"
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "REML",         # restricted maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          struct = "CS",
                          mods = ~ dependent_variable:robot_type_classification - 1)
  return(Multi_REModel)  
}

model_x <- function(db, V){
  "This function computes the multivariate mixed-effect analysis with a different
  method to reach convergency
  Moderator 2: Robot Type"
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "ML",           # maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          struct = "CS",
                          mods = ~ dependent_variable:robot_type_classification - 1)
  return(Multi_REModel)  
}

multivariate_COR_moderatorInteractionType_model <- function(db, V){
  "This function computes the multivariate mixed-effect analysis
  Moderator 3: Interaction Type"
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "REML",         # restricted maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          struct = "CS",
                          mods = ~ dependent_variable:interaction_type -1)
  return(Multi_REModel) 
}

multivariate_COR_moderatorInteractionLength_model <- function(db, V){
  "This function computes the multivariate mixed-effect analysis
  Moderator 4: Interaction Length"
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "REML",         # restricted maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          # struct = "CS",
                          mods = ~ dependent_variable:avg_length -1)
  return(Multi_REModel) 
}

multivariate_COR_moderatorRobotOperation_model <- function(db, V){
  "This function computes the multivariate mixed-effect analysis
  Moderator 5: Robot Operation"
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "REML",         # restricted maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          struct = "CS",
                          mods = ~ dependent_variable:robot_operation -1)
  return(Multi_REModel) 
}

multivariate_COR_moderatorIV_model <- function(db, V){
  "This function computes the multivariate mixed-effect analysis
  Moderator 6: Independent Variable Category"
  Multi_REModel <- rma.mv(db$d_calc, 
                          V,
                          data = db,
                          W = db$weights_d,
                          method = "REML",         # restricted maximum likelihood
                          level = 95,              # Confidence Internval 
                          digits = 7,              # decimal points
                          slab = db$short_cite,    # study labels
                          random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                          struct = "CS",
                          mods = ~ dependent_variable:IV_category -1)
  return(Multi_REModel) 
}

multivariate_COR_typemeasure_model <- function(db, V){
  "This function computes the multivariate mixed-effect analysis
  Moderator 7: Type of Measure"
  Multi_REModel <- rma.mv(db$d_calc, 
                    V,
                    data = db,
                    W = db$weights_d,
                    method = "REML",         # restricted maximum likelihood
                    level = 95,              # Confidence Internval 
                    digits = 7,              # decimal points
                    slab = db$short_cite,    # study labels
                    random = ~factor(outcome)| short_cite, # multilevel model to handle sample dependency
                    struct = "CS",
                    mods = ~ dependent_variable:type_measure -1)
 return(Multi_REModel) 
}
