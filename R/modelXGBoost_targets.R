tar_option_set(packages = c("rsample", "gbm", "xgboost", "h2o", "pdp", "lime",
                            "vtreat","tidyverse", "lubridate"))


modelXGBoost_target <-
  tar_plan(
    # http://uc-r.github.io/gbm_regression
    
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
        head(15),
    
    
    # Fit the top model
    # parameter list
    params_W_target = list(
      eta = topModel_W_target$eta[1], #0.01,
      max_depth = topModel_W_target$max_depth[1], #7,
      min_child_weight = topModel_W_target$min_child_weight[1], #5,
      subsample = topModel_W_target$subsample[1], #0.65,
      colsample_bytree = topModel_W_target$colsample_bytree[1] #0.8
    ),
    
    # train final model
    xgb.fit.final_W_target = xgboost(
      params = params_W_target,
      data = features_train_W_target,
      label = response_train_W_target,
      nrounds = topModel_W_target$optimal_trees[1],
      objective = "reg:squarederror",
      verbose = 0
    ),
    
    # create importance matrix
    importance_matrix_W_target = xgb.importance(model = xgb.fit.final_W_target)
    
    
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