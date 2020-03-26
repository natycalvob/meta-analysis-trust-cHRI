# This script is based on the MetaLab Repository: https://github.com/langcog/metalab2

compute_es_and_esvar <- function(db){
  "This funtion computes the effect sizes and the variance. 
   The effect sizes computed in this example reflect standardized mean differences (Cohen's d-type effect size). 
   These effect sizes can be computed based on either means/SDs, t-values, or F values. This script checks which 
   of these values are present (prioritizing computation from means/SDs), and computes effect sizes whereever possible.
   Start of decision tree where effect sizes are calculated differently based on participant design depending on which 
   data is available, effect sizes are calculated differently. "
  
  # Initiate the columns for d and d_var that we will later fill with the effect sizes and their variances. 
  # We call them d_calc and d_var_calc to distinguish them from d and d_var for cases where these values were reported in primary articles
  db$d_calc <- NA
  db$d_var_calc <- NA
  db$es_method <- "missing"
  # db$pooled_SD <- NA
  # print(db$es_method)
  
  # Start of decision tree where effect sizes are calculated differently based on participant design
  # depending on which data is available, effect sizes are calculated differently
  for(line in 1:length(db$n_1)){
    # Participant Desing: Between
    if(db[line,]$participant_design == "between"){
      db[line,]$es_method <- "between"
      # Effect size calculation 'd_calc'
      if (complete.cases(db[line,]$x_1, db[line,]$x_2, db[line,]$SD_1, db[line,]$SD_2)) {
        pooled_SD <- sqrt(((db[line,]$n_1 - 1) * db[line,]$SD_1 ^ 2 + (db[line,]$n_2 - 1) * db[line,]$SD_2 ^ 2) / (db[line,]$n_1 + db[line,]$n_2 - 2)) # Lipsey & Wilson, 3.14
        db[line,]$d_calc <- (db[line,]$x_1 - db[line,]$x_2) / pooled_SD # Lipsey & Wilson (2001)
        db[line,]$es_method <- "between_means"
      } else if (complete.cases(db[line,]$t)) {
        db[line,]$d_calc <- db[line,]$t * sqrt((db[line,]$n_1 + db[line,]$n_2) / (db[line,]$n_1 * db[line,]$n_2)) # Lipsey & Wilson, (2001)
        db[line,]$es_method <- "between_t"
      } else if (complete.cases(db[line,]$f)) {
        db[line,]$d_calc <- sqrt(db[line,]$f * (db[line,]$n_1 + db[line,]$n_2) / (db[line,]$n_1 * db[line,]$n_2)) # Lipsey & Wilson, (2001)
        db[line,]$es_method <- "between_f"
      }
      # Effect size variance calculation 'd_var_calc'
      if(complete.cases(db[line,]$n_1, db[line,]$n_2, db[line,]$d_calc)){
        db[line,]$d_var_calc <- ((db[line,]$n_1 + db[line,]$n_2) / (db[line,]$n_1 * db[line,]$n_2)) + (db[line,]$d_calc ^ 2 / (2 * (db[line,]$n_1 + db[line,]$n_2)))
      } else if (complete.cases(db[line,]$r, db[line,]$r_var)){
        # if r -Pearson's correlation- was reported instead of d -Cohen's correlation- transform it for standardization
        db[line,]$d_calc <- 2*db[line,]$r / sqrt(1 - db[line,]$r ^2)
        db[line,]$d_var_calc <- 4*db[line,]$r_var / ((1-db[line,]$r ^2)^3)
        db[line,]$es_method <- "between_r"
      } else if (complete.cases(db[line,]$d, db[line,]$d_var)){
        db[line,]$d_calc <- db[line,]$d
        db[line,]$d_var_calc <- db[line,]$d_var
        db[line,]$es_method <- "betweeb_d"
      }
    }
    # Participant Desing: Within_two
    else if(db[line,]$participant_design == "within_two") {
      db[line,]$es_method <- "within_two"
      # Effect size calculation 'd_calc'
      if(complete.cases(db[line,]$x_1, db[line,]$x_2, db[line,]$SD_1, db[line,]$SD_2 )){
        pooled_SD <- sqrt((db[line,]$SD_1^2 + db[line,]$SD_2^2) / 2)      # Lipsey & Wilson (2001)
        db[line,]$d_calc <- (db[line,]$x_1 - db[line,]$x_2) / pooled_SD   # Lipsey & Wilson (2001)
        db[line,]$es_method <- "group_means_two"
      } else if (complete.cases(db[line,]$x_1, db[line,]$x_2, db[line,]$SD_dif)) {
        within_SD <-db[line,]$SD_dif / sqrt(2 * (1 - db[line,]$corr)) # Lipsey & Wilson (2001); Morris & DeShon (2002)
        db[line,]$d_calc <- (db[line,]$x_1 - db[line,]$x_2) / within_SD # Lipsey & Wilson (2001)
        db[line,]$es_method <- "group_means_two"
      } else if(complete.cases(db[line,]$x_dif, db[line,]$SD_1, db[line,]$SD_2)){
        pooled_SD <- sqrt((db[line,]$SD_1^2 + db[line,]$SD_2^2) / 2)      # Lipsey & Wilson (2001)
        db[line,]$d_calc <- (db[line,]$x_dif / pooled_SD)
        db[line,]$es_method <- "subj_diff_two"
      } else if (complete.cases(db[line,]$x_dif, db[line,]$SD_dif)){
        wc <- sqrt(2*(1-db[line,]$corr))
        db[line,]$d_calc <- (db[line,]$x_dif / db[line,]$SD_dif) * wc    # Morris & DeShon (2002)
        db[line,]$es_method <- "subj_diff_two"
      } else if(complete.cases(db[line,]$t)){
        wc <- sqrt(2*(1-db[line,]$corr))                                 # Dunlap et al. 1996 p.171
        db[line,]$d_calc <- (db[line,]$t / sqrt(db[line,]$n_1)) * wc
        db[line,]$es_method <- "t_two"
      } else if(complete.cases(db[line,]$f)){
        wc <- sqrt(2*(1-db[line,]$corr))
        db[line,]$d_calc <- sqrt( db[line,]$f / db[line,]$n_1)*wc
        db[line,]$es_method <- "f_two"
      }
      # Effect size variance calculation 'd_var_calc'
      if(complete.cases(db[line,]$n_1, db[line,]$d_calc)){
        # We need the correlations to compute the variance of the effect size Double check with Sho
        # db[line,]$d_var_calc <- (2*(1-db[line,]$corr) / db[line,]$n_1) + (db[line,]$d_calc^2 / (2*db[line,]$n_1))  # Lipsey & Wilson (2001)
        db[line,]$d_var_calc <- (2/ db[line,]$n_1) + (db[line,]$d_calc^2 / (2*db[line,]$n_1))  # Lipsey & Wilson (2001)
      } else if(complete.cases(db[line,]$r)){
        # if r -Pearson's correlation- was reported instead of d -Cohen's correlation- transform it for standardization
        db[line,]$d_calc <- 2*db[line,]$r / sqrt(1-db[line,]$r^2)
        db[line,]$d_var_calc <- 4*db[line,]$r_var / ((1-db[line,]$r^2)^3)
        db[line,]$es_method <- "r_two"
      } else if(complete.cases(db[line,]$d, db[line,]$d_var)){
        db[line,]$d_calc <- db[line,]$d
        db[line,]$d_var_calc <- db[line,]$d_var
        db[line,]$es_method <- "d_two"
      }
    }
    # Participant Desing: Within_one
    else if(db[line,]$participant_design == "within_one"){
      # Effect size calculation 'd_calc'
      # This is super important, x2 is supposed to contain chance level where applicable, 0 where not. 
      if (complete.cases(db[line,]$x_1, db[line,]$x_2, db[line,]$SD_1)) {
        db[line,]$d_calc <- (db[line,]$x_1 -db[line,]$x_2) / db[line,]$SD_1
        db[line,]$es_method  <- "group_means_one"
      } else if (complete.cases(db[line,]$t)) {
        db[line,]$d_calc <- db[line,]$t / sqrt(db[line,]$n_1)
        db[line,]$es_method  <- "t_one"
      } else if(complete.cases(db[line,]$f)){
        db[line,]$d_calc <- sqrt(db[line,]$f / db[line,]$n_1)
        db[line,]$es_method <- "f_one"
      }
      # Effect size variance calculation 'd_var_calc'
      if (complete.cases(db[line,]$n_1, db[line,]$d_calc)) {
        # This models what is done in metafor package, escalc(measure="SMCR"() (Viechtbauer, 2010)
        db[line,]$d_var_calc <- (2 / db[line,]$n_1) + (db[line,]$d_calc ^ 2 / (2 * db[line,]$n_1))
      } else if(complete.cases(db[line,]$n_1, db[line,]$d_calc)){
        # Get variance of transformed r (z, Fisher's transformation)
        SE_z = 1 / sqrt(db[line,]$n_1 - 3)     # Howell (2010, Statistical methods for Psychology, pg 275)
        var_z = SE_z ^ 2 
        # Transform z variance to r variance 
        var_r = tanh(var_z) # from wikipedia (https://en.wikipedia.org/wiki/Fisher_transformation) for convert z -> r, consistent with Howell
        # Transform r to d 
        db[line,]$d_calc <- 2*db[line,]$r / (sqrt(1-db[line,]$r ^2))   # from Hunter and Schmidt, pg 279
        db[line,]$d_var_calc <- (4*var_r)/(1-db[line,]$r^2)^3          # from https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf (pg. 4)
        db[line,]$es_method <- "r_one"
      } else if(complete.cases(db[line,]$r)){
        # if r -Pearson's correlation- was reported instead of d -Cohen's correlation- transform it for standardization
        db[line,]$d_calc <- 2*db[line,]$r / sqrt(1 - db[line,]$r^2)
        db[line,]$d_var_calc <- 4*db[line,]$r_var / ((1 - db[line,]$r^2)^3)
        db[line,]$es_method <- "r_one"
      } else if(complete.cases(db[line,]$d))
        # if d and d_var were already reported 
        db[line,]$d_calc <- db[line,]$d
      # this models what is done in metafor package, escalc(measure="SMCR"() (Viechtbauer, 2010)
      db[line,]$d_var_calc <- (1/db[line,]$n_1) + (db[line,]$d_calc^2/(2*db[line,]$n_1))
      db[line,]$es_method <- "d_one"
    }
  }
  return(db)
}


debug_incomplete_es <- function(db){
  "This function debugs the dataframe and look for missing information.
  We need to figure out the reason why we could not compute the effect sizes"
  
  db_miss_d <- db[!complete.cases(db$d_calc),]
  db_miss_d_var <- db[!complete.cases(db$d_var_calc),]
  cat('The effect sizes could not bet computed \n')
  db_miss_d_var$short_cite
  
  return(db_miss_d)
}

remove_missing_EF <-function(db){
  "This function removes missing data and outliers
  We removed effect sizes more than 3DS away from the median effect,
  as proposed by Stanford Metalab"
  
  db = db[!is.na(db$d_calc),]
  db$nooutlier = ifelse(db$d_calc > mean(db$d_calc, na.rm = TRUE) + 3*sd(db$d_calc, na.rm = TRUE) 
                        | db$d_calc < mean(db$d_calc, na.rm = TRUE) - 3*sd(db$d_calc, na.rm = TRUE),FALSE, TRUE)
  
  db = db[db$nooutlier,]
  
  # For our MA we don't have completed data for mean_age. For some studies we have the range, we approximated the avarage. 
  # We compute weighted mean age with approximations instead. 
  
  db = db %>%
    rowwise() %>%
    mutate(mean_age = weighted.mean(c(mean_age_1, mean_age_2), c(n_1, n_2), na.rm = TRUE),
           n_mean = mean(c(n_1, n_2), na.rm = TRUE), 
           weights_d = 1/(d_var_calc)^2) %>%
    rownames_to_column("unique_row") %>%
    ungroup()

  return(db)
}