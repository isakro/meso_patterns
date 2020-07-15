# Function to impute data 5 times and then run each set through a nested cv
# setup using random forest
imp_nest_rf <- function(input_data, nimp, ntrees,
                        depend, nindep, indep = NULL){
  results <- list()
  
  # If island variables are specified as the predictors, pass those on with the 
  # dependent on to nested cv.  
  if(!(is.null(indep))){
    rf <- nestcv_func(input_data[, names(input_data) %in% c(depend, indep)],
                      ntrees, depend, nindep)
    return(rf)
  }
  
  # If not then 
  else{
    # Loop over each imputed dataset
    for(i in 1:nimp){
      # Computation can take a while - this is for keeping track
      print(paste("Imputation", i, "started", Sys.time()), sep = " ")
      
      # Impute
      dat <- rfImpute(as.formula(paste(depend, '~.')), input_data,
                      ntree = 5000)
      dat[, names(dat) == 'infil'] <- 
        as.numeric(as.character(dat[, names(dat) == 'infil']))
      
      # Pass to nested cv function
      rf <- nestcv_func(dat[, names(dat) != 'elev'], ntrees, depend, nindep)
      results[[i]] <- rf
  }
  return(results)
}
}

# Function to measure Brier and AUC scores accross cross-validation folds 
# by passing this to the summaryFunction parameter of carets trainControl().
# (partly repurposed from print(twoClassSummary))
summary_func <- function(data, lev = NULL, model = NULL) {
  brier <- mean((as.numeric(data$pred) - as.numeric(data$obs))^2)  
  brier_all <- sum(brier)
  
  lvls <- levels(data$obs)
  rocAUC <- auc(ifelse(data$obs == lev[2], 0,
                       1), data[, lvls[1]])
  out <- c(brier_all, rocAUC)
  
  names(out) <- c("Brier", "ROC")
  return(out)
}

# Nested 5-fold cross-validation. Tuning in inner cv cycle using AUC to find
# optimal mtry. ntree needs to be identified beforehand. Performance
# measure by Brier and AUC, as well as variable importance measures are done 
# in the outer cv cycle.
# Main sources of inspiration for syntax: 
# https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/
# https://machinelearningmastery.com/tuning-machine-learning-models-using-
# the-caret-r-package/
nestcv_func <- function(imp_data, ntrees, depend, nindep){
  # Create outer folds
  folds <- createFolds(imp_data[, depend], k = 5)
  # Loop over folds
  res <-  sapply(folds, function(fold){
    
    # Set up number of inner folds and define summary statistic
    ctrl <- trainControl(method = 'cv', number = 5,
                         classProbs = TRUE, summaryFunction = twoClassSummary)
    
    # Set up grid (vector in this case) search for optimal mtry
    tgrid <- expand.grid(.mtry = c(1:nindep)) 
    
    # Run current outer folds  through inner cv cycle
    rf_grid <- train(as.formula(paste(depend, '~.')),
                     data = imp_data[-fold,], method = 'rf', metric = 'ROC',
                     maximize = TRUE, ntree = ntrees, tuneGrid = tgrid,
                     trControl = ctrl)
    
    # Retrieve mtries accross the cv cycle, and find the best one
    mtries <- rf_grid$results[rf_grid$results$ROC == max(rf_grid$results$ROC),
                              names(rf_grid$results) %in% c('mtry', 'ROC')]
    nmtry <- mtries[[1]]
    
    # Use this mtry to run rf on the retained outer fold
    mod <- randomForest(as.formula(paste(depend, '~.')),
                        data = imp_data[fold,], ntree = ntrees, mtry = nmtry,
                        importance = TRUE)

    # Retrieve measures of performance and variable importance
    bri <- brier(mod)
    auc <- auc(mod)
    varimp <- importance(mod, scale = FALSE, type = 1)
    
    # Return result for current outer fold
    result <- list(bri, auc, varimp)
    return(result)
  })
  
  # Return results for all outer folds
  return(res)
}

# Function to scale continous variables to values between 0 and 1. 
scale_func <- function(x){
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Function for bootstrap resampling with logistic regression  
# and single imputation of null-values.
logistic_reg <- function(input_data, i, dependent, independent = NULL, ...){
  # Single imputation using all varaibles
  imp <- mice(input_data, m = 1, maxit = 20, print = FALSE)
  d <- complete(imp)
  
  # Resampling
  train_data <- d[i, ]
  valid_data <- d
  
  # If specified independent variables
  if(!(is.null(independent))){
    lrm <- glm(as.formula(c(paste(dependent, '~', paste(independent,
                                                        collapse = '+')))), 
               data = train_data, family = binomial(logit))
  }
  else{
    # If not, use all except elev, which is only used for imputation
    lrm <- glm(as.formula(paste(dependent, ' ~. -elev', sep='')), 
             data = train_data, family = binomial(logit))
  }
  
  # Find Brier and AUC score
  pred_prob <- predict(lrm, newdata = valid_data, type = 'response')
  auc_score <- auc(valid_data[,dependent], pred_prob)
  brier_score <- brier(as.numeric(valid_data[,dependent])-1, pred_prob) 
  scores <- c(brier_score, auc_score)
  names(scores) <- c('brier', 'auc')
  
  # Find VIFs
  varif <- vif(lrm)
  names(varif) <- paste('vif_', names(coef(lrm)[names(coef(lrm))
                                                != '(Intercept)']), sep = '')
  # Return data
  results <- c(coef(lrm), scores, varif)
  return(results)
}

# Function to find percentile rank for each value in df1, 
# based on the corresponding column in df2. To be used for
# comparison accross phases. Make sure the columns match before
# passing to the function.
prcr_func <- function(df1, df2, cols){
  df1 <- df1[,names(df1) %in% cols]  
  df2 <- df2[,names(df2) %in% cols]   
  
  prcr <- list()
  for(i in seq(1, ncol(df1), 1)){
    target <- df1[, i]
    base <- df2[, i]
    pr <- ecdf(base)(target) * 100
    prcr[[i]] <- pr
  }
  do.call(cbind, prcr)
  names(prcr) <- names(df1)
  return(prcr)
}
