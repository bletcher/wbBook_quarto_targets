tar_option_set(packages = c("rsample", "gbm", "xgboost", "h2o", "pdp", "lime",
                            "vtreat","tidyverse", "lubridate"))


modelXGBoost_target <-
  tar_plan(
    # http://uc-r.github.io/gbm_regression
    numModelsToSave_target = 15,
    
    ####################
    # xgboost model for growth in mass ('xxxx_W')
    ####################
    
    dML_W_target = cd1_target |> 
      filter(speciesGG == "Brook trout", !is.na(grWeight)) |> 
      mutate(age = year - cohort) |> 
      dplyr::select(riverGG, seasonGG, observedWeight, cohort, age, meanTemperature, meanFlowByRiver, grWeight), 
    
    # variable names
    features_W_target = setdiff(names(dML_W_target), "grWeight"),
    
    # Create the treatment plan from the training data
    treatplan_W_target = vtreat::designTreatmentsZ(dML_W_target, features_W_target, verbose = FALSE),
    
    # Get the "clean" variable names from the scoreFrame
    new_vars_W_target = treatplan_W_target %>%
      magrittr::use_series(scoreFrame) %>%        
      dplyr::filter(code %in% c("clean", "lev")) %>% 
      magrittr::use_series(varName) ,    
    
    # Prepare the training data
    features_train_W_target = vtreat::prepare(treatplan_W_target, dML_W_target, varRestriction = new_vars_W_target) %>% as.matrix(),
    response_train_W_target = dML_W_target$grWeight,
    
    # Prepare the test data
    features_test_W_target = vtreat::prepare(treatplan_W_target, dML_W_target, varRestriction = new_vars_W_target) %>% as.matrix(),
    response_test_W_target = dML_W_target$grWeight,
    
    # hypergrid for training model
    hyper_grid_W_target0 = expand.grid(
      eta = c(.01, .05),
      max_depth = c(3, 5, 7),
      min_child_weight = c(1, 3, 5, 7),
      subsample = c(.65, .8), 
      colsample_bytree = c(.8, .9, 1),
      optimal_trees = 0,               # a place to dump results
      min_RMSE = 0                     # a place to dump results
    ),
    
    hyper_grid_W_target = runHyperGrid(hyper_grid_W_target0, features_test_W_target, response_test_W_target),
    
    topModel_W_target = hyper_grid_W_target %>%
        dplyr::arrange(min_RMSE) %>%
        head(15), #(numModelsToSave_target), UPDATE this when ready to rerun hyperGrid
    
    finalModels_W_target = runFinalModels(
      topModel_W_target, 
      numModelsToSave_target,
      features_train_W_target,
      response_train_W_target
    )

  )


# functions
runHyperGrid <- function(d, features, response) {
  for(i in 1:nrow(d)) {
    print(i)
    # create parameter list
    params_W <- list(
      eta = d$eta[i],
      max_depth = d$max_depth[i],
      min_child_weight = d$min_child_weight[i],
      subsample = d$subsample[i],
      colsample_bytree = d$colsample_bytree[i]
    )
    
    # train model
    xgb.tune_W <- xgb.cv(
      params = params_W,
      data = features,
      label = response,
      nrounds = 5000,
      nfold = 5,
      objective = "reg:squarederror",  # for regression models
      verbose = 0,               # silent,
      early_stopping_rounds = 10 # stop if no improvement for 10 consecutive trees
    )
    
    # add min training error and trees to grid
    d$optimal_trees[i] <- which.min(xgb.tune_W$evaluation_log$test_rmse_mean)
    d$min_RMSE[i] <- min(xgb.tune_W$evaluation_log$test_rmse_mean)
  }
  return(d)
}

runFinalModels <- function(t = topModel_W_target, 
                           n = numModelsToSave_target, 
                           f = features_train_W_target,
                           r = response_train_W_target) {
  
  params <- xgb <- importanceMatrix <- list()
  # Fit the top model
  # parameter list
  for(i in 1:n) {
    print(i)
    
    params[[i]] <- list(
      eta = t$eta[i], #0.01,
      max_depth = t$max_depth[i], #7,
      min_child_weight = t$min_child_weight[i], #5,
      subsample = t$subsample[i], #0.65,
      colsample_bytree = t$colsample_bytree[i] #0.8
    )
    
    # train final model
    xgb[[i]] <- xgboost(
      params = params[[i]],
      data = f,
      label = r,
      nrounds = t$optimal_trees[i],
      objective = "reg:squarederror",
      verbose = 0
    )
    
    # create importance matrix
    importanceMatrix[[i]] <- xgb.importance(model = xgb[[i]])
  }
  return(list(params = params, xgb = xgb, importanceMatrix = importanceMatrix
              ))
  
}