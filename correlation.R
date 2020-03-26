library(dplyr)

compute_correlations <- function(db){
  # For within_two if (x_1 && x_3 && sd_1 && sd_2 && t) then compute correlations 
  for(line in 1:length(db$n_1)){
    if(db[line,]$participant_design == "within_two") {
      if(complete.cases(db[line,]$corr)){
        db[line,]$corr = db[line,]$corr
      } else if (complete.cases(db[line,]$x_1,db[line,]$x_2,db[line,]$SD_1,db[line,]$SD_2,db[line,]$t)){
        db[line,]$corr = (db[line,]$SD_1^2 + db[line,]$SD_2^2 - (db[line,]$n_1 * (db[line,]$x_1 - db[line,]$x_2)^2 / db[line,]$t^2)) / (2 * db[line,]$SD_1 * db[line,]$SD_2)
      }
    }
  }  
  
  # We impute the remaining missing correlations.
  # First we calculate the median & variance of existing correlatio 
  median_corr = median(db$corr, na.rm = TRUE)  # remove when there is not info
  var_corr = var(db$corr, na.rm = TRUE)
  
  # We then impute missing correlations by drawing randomly from a normal distribution with that median and variance
  set.seed(111)
  db$corr_imputed = rnorm(length(db$n_1), mean = median_corr, sd = var_corr)
  db$corr = ifelse(is.na(db$corr), db$corr_imputed, db$corr)
  
  # Correlations cannot be higher than 1 or lower than 0, fixing some possible calculation and imputation issues. 
  db$corr = ifelse(db$corr > .99, .99, db$corr)
  db$corr = ifelse(db$corr < .01, .01, db$corr)
  
  return(db)
}