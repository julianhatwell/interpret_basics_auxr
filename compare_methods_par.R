rm(list = ls())
library(psych)
library(ada)
library(randomForest)
library(gbm)
library(foreach)
library(parallel)
library(doParallel)
random_states <- 123 # :152
source("compare_methods_utils.R")

max_tests <- 1000

# model <- "rf"
model <- "samme"
# model <- "gbm"

# algorithm <- "inTrees"
algorithm <- "BRL"

results_nrows <- length(random_states) * length(datasetnames)

# adult 
# rnr <- 1; ntree_divisor <- 5; inTrees_maxdepth <- 8; lambda <- 10; eta <- 1; rule_maxlen <- 8; nchain <- 10

# bankmark
# rnr <- 2; ntree_divisor <- 20; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 4; nchain <- 10

# car
# rnr <- 3; ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10

# cardio
# rnr <- 4; ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 4; nchain <- 10

# credit
# rnr <- 5; ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 8; nchain <- 10

# german
# rnr <- 6; ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 4; nchain <- 10

# lending_tiny_samp CANNOT BE DONE IN R

# nursery
# rnr <- 8; ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10

# rcdv
rnr <- 9; ntree_divisor <- 10; inTrees_maxdepth <- 8; lambda <- 10; eta <- 1; rule_maxlen <- 4; nchain <- 10

# set counters
r <- ((rnr - 1) %/% length(datasetnames)) + 1
i <- ((rnr - 1) %% length(datasetnames)) + 1

# encapsulated set up
data_prep(i, max_tests = max_tests)

# classifier prep
fmla <- as.formula(paste(class_cols[i], "~ ."))

ntree <- fromJSON(readLines(file(paste0(
  project_dir
  , datasetnames[i]
  , pathsep
  , model
  , "_best_params_rnst_"
  , random_states[r]
  , ".json"))))$n_estimators

if (model == "rf") {
  maxnodes <- 2^(fromJSON(readLines(file(paste0(
    resfilesdirs[i]
    , model
    , "_best_params_rnst_"
    , random_states[r]
    , ".json"))))$max_depth - 1)
}

if (model %in% c("samme", "samme.r")) {
  bestparams <- fromJSON(readLines(file(paste0(
    project_dir
    , datasetnames[i]
    , pathsep
    , model
    , "_best_params_rnst_"
    , random_states[r]
    , ".json"))))
  ntree <- bestparams$n_estimators
  max_depth_pos <- regexpr("max_depth=", bestparams$base_estimator) + nchar("max_depth=")
  max_depth <- as.integer(substr(bestparams$base_estimator, max_depth_pos, max_depth_pos))
}

if (model == "gbm") {
  bestparams <- fromJSON(readLines(file(paste0(
    project_dir
    , datasetnames[i]
    , pathsep
    , model
    , "_best_params_rnst_"
    , random_states[r]
    , ".json"))))
}

# build forest
set.seed(random_states[r])
if (model == "rf") forest <- randomForest(fmla, data = dat_train
                                          , ntree = ntree
                                          , maxnodes = maxnodes)
if (model == "samme") forest <- ada(fmla, data = dat_train, iter = ntree, nu = 1
                                  , loss = "exponential", type = "discrete"
                                  , rpart.control(maxdepth = max_depth))
if (model == "samme.r") forest <- ada(fmla, data = dat_train, iter = ntree, nu = 1
                                    , loss = "exponential", type = "real"
                                    , rpart.control(maxdepth = max_depth))
if (model == "gbm") {
  if (!is.na(positive_classes[i])) {
    dat_train[, class_cols[i]] <- ifelse(dat_train[, class_cols[i]] == positive_classes[i], 1, 0)
    dat_test[, class_cols[i]] <- ifelse(dat_test[, class_cols[i]] == positive_classes[i], 1, 0)
  }

  forest <- gbm(fmla, data = dat_train
                , n.trees = base_estimator$n_estimators
                , interaction.depth = base_estimator$max_depth
                , shrinkage = base_estimator$learning_rate
                , bag.fraction = base_estimator$subsample)
  }

# run method and benchmark
if (algorithm == "inTrees") {
  benchmark <- inTrees_benchmark(forest = forest
                                 , ds_container = ds_container
                                 , ntree = ntree / ntree_divisor
                                 , maxdepth = inTrees_maxdepth
                                 , model = model
  )
} else {
  benchmark <- sbrl_benchmark(ds_container=ds_container
                              , classes = classes
                              , lambda = lambda
                              , eta = eta
                              , rule_maxlen = rule_maxlen
                              , nchain = nchain)
}

# process forest preds
if (model == "rf") {
  forest_label <- predict(forest
                          , newdata = ds_container$X_test)
  forest_preds <- which_class(as.character(forest_label))
}
if (model == "gbm") {
  forest_label <- predict(forest
                           , newdata = ds_container$X_test
                           , n.trees = base_estimator$n_estimators
                           , type = "response")
  if (!is.na(positive_classes[i])) {
    forest_label <- factor(ifelse(forest_label > 0.5, positive_classes[i], negative_classes[i]))
    forest_label <- relevel(forest_label, positive_classes[i])
  } else {
    forest_label <- apply(forest_label, 1, which.max)
  }
}
if (model %in% c("samme", "samme.r")) {
  forest_label <- predict(forest
                          , newdata = ds_container$X_test)
  forest_preds <- which_class(as.character(forest_label))
  
  # this is the best way to get forest vote share
  # Initiate cluster
  cl <- makeCluster(detectCores() - 2)
  clusterEvalQ(cl, library("ada"))
  n_est <- bestparams$n_estimators
  frst <- forest
  X_test <- ds_container$X_test[1:n_test, ]
  frst_prds <- forest_preds[1:n_test]
  clusterExport(cl, c("n_est", "frst", "X_test", "frst_prds"))
  staged_preds <- parSapply(cl, 1:n_est, function(t) {
    sp <- predict(frst, newdata=X_test
                  , n.iter = t
                  , type = "probs")[cbind(1:length(frst_prds), frst_prds)]
    return(sp)
  })
  tree_preds <- parSapply(cl, 1:n_est, function(t) {
    tp <- predict(frst$model$trees[[t]], newdata = X_test)
    return(apply(tp, 1, which.max))
  })
  stopCluster(cl)
  
  sp <- sapply(1:n_test, function(j) log(staged_preds[j, ]/(1 - staged_preds[j, ])) + log(n_classes[i] - 1))
  mx <- apply(sp, 2, function(j) max(j[j < Inf]))
  rws <- apply(sp, 2, function(j) j == Inf)
  sapply(1:n_test, function(j) sp[rws[, j]] <<- mx[j])
  sp <- sapply(1:n_test, function(j) 1 + (sp[, j] - min(sp[, j]))/ max(sp[, j]))
  sp <- sapply(1:n_test, function(j) sp[, j] / sum(sp[, j]))
  forest_vote_share <- sapply(1:n_test, function(j) sum(sp[, j][tree_preds[j, ] == forest_preds[j]]))

}

# collect     results
tpe <- time_per_explanation(benchmark$begin_time
                            , benchmark$completion_time
                            , n_test)
results_init(n_test)
for (j in 1:n_test) {
  rule <- benchmark$rule[j]
  
  covered_instances <- apply_rule(rule, ds_container$X_test)
  covered_instances <- covered_instances[-j] # drop out current instance
  loo_true_labels <- ds_container$y_test[-j]

  instance_results <- evaluate(prior_labels = loo_true_labels
                               , post_idx = covered_instances
                               , classes = classes)

  if (model == "rf") {
    forest_vote_share[j] <-
      mean(predict(forest, newdata=ds_container$X_test[j, ]
                   , predict.all = TRUE)$individual ==
             classes[forest_label[j]])
  }
  if (model == "gbm") {
    staged_preds <- numeric(bestparams$n_estimators + 1)
    for (t in 1:(bestparams$n_estimators + 1)) {
      staged_preds[t] <- predict(forest, newdata=ds_container$X_test[j, ]
                                 , n.trees = t - 1
                                 , type = "link")
    }
    staged_preds <- diff(staged_preds, 1)
    sign_class <- ifelse(forest_label[j] %in% positive_classes, 1, -1)
    agree_sign <- sum(ifelse(sign_class == sign(staged_preds), abs(staged_preds), 0))
    disagree_sign <- sum(ifelse(sign_class != sign(staged_preds), abs(staged_preds), 0))
    forest_vote_share[j] <- agree_sign / (agree_sign + disagree_sign)
  }

  prior[j] <- instance_results[["prior"]][forest_label[j]]
  coverage[j] <- instance_results[["coverage"]]
  xcoverage[j] <- instance_results[["xcoverage"]]
  proxy_precision[j] <- instance_results[["posterior"]][benchmark$label[j]]
  proxy_stability[j] <- instance_results[["stability"]][benchmark$label[j]]
  proxy_recall[j] <- instance_results[["recall"]][benchmark$label[j]]
  proxy_f1[j] <- instance_results[["f1"]][benchmark$label[j]]
  proxy_cc[j] <- instance_results[["cc"]][benchmark$label[j]]
  proxy_ci[j] <- instance_results[["ci"]][benchmark$label[j]]
  proxy_ncc[j] <- instance_results[["ncc"]][benchmark$label[j]]
  proxy_nci[j] <- instance_results[["nci"]][benchmark$label[j]]
  proxy_npv[j] <- instance_results[["npv"]][benchmark$label[j]]
  proxy_accu[j] <- instance_results[["accu"]][benchmark$label[j]]
  proxy_lift[j] <- instance_results[["lift"]][benchmark$label[j]]
  proxy_kl_div[j] <- instance_results[["kl_div"]][benchmark$label[j]]
  forest_precision[j] <- penalise_bad_prediction(forest_preds[j]
                                                 , benchmark$label[j]
                                                 , instance_results[["posterior"]][forest_preds[j]])
  forest_stability[j] <- penalise_bad_prediction(forest_preds[j]
                                                 , benchmark$label[j]
                                                 , instance_results[["stability"]][forest_preds[j]])
  forest_recall[j] <- penalise_bad_prediction(forest_preds[j]
                                              , benchmark$label[j]
                                              , instance_results[["recall"]][forest_preds[j]])
  forest_f1[j] <- penalise_bad_prediction(forest_preds[j]
                                          , benchmark$label[j]
                                          , instance_results[["f1"]][forest_preds[j]])
  forest_accu[j] <- penalise_bad_prediction(forest_preds[j]
                                            , benchmark$label[j]
                                            , instance_results[["accu"]][forest_preds[j]])
  forest_kl_div[j] <- penalise_bad_prediction(forest_preds[j]
                                            , benchmark$label[j]
                                            , instance_results[["kl_div"]][forest_preds[j]])
  forest_npv[j] <- penalise_bad_prediction(forest_preds[j]
                                              , benchmark$label[j]
                                              , instance_results[["npv"]][forest_preds[j]])
  forest_lift[j] <- ifelse(forest_precision[j] == 0, 1, proxy_lift[j])
}

# save results to file before exiting loop
write.csv(data.frame(dataset_name = rep(datasetnames[i], n_test)
                     , instance_id = test_idx[1:n_test] - 1 # conform with Python zero base
                     , algorithm = rep(algorithm, n_test)
                     , pretty_rule = benchmark$rule[1:n_test]
                     , rule_length = benchmark$rl_ln[1:n_test]
                     , true_class = as.numeric(ds_container$y_test[1:n_test]) - 1 # conform with Python zero base
                     , true_class_label = as.character(ds_container$y_test[1:n_test]) 
                     , pred_class = as.numeric(forest_label[1:n_test]) - 1 # conform with Python zero base
                     , pred_class_label = classes[forest_label][1:n_test]
                     , target_class = benchmark$label[1:n_test] - 1  # conform with Python zero base
                     , target_class_label = classes[benchmark$label][1:n_test]
                     , forest_vote_share = forest_vote_share
                     , accummulated_points = 0 # only meaningful for CHIRPS
                     , prior = prior
                     , precision_tr = proxy_precision
                     , stability_tr	= proxy_stability
                     , recall_tr = proxy_recall
                     , f1_tr = proxy_f1
                     , cc_tr = proxy_cc
                     , ci_tr = proxy_ci
                     , ncc_tr = proxy_ncc
                     , nci_tr = proxy_nci
                     , npv_tr = proxy_npv
                     , accuracy_tr = proxy_accu
                     , lift_tr = proxy_lift
                     , coverage_tr = coverage
                     , xcoverage_tr = xcoverage
                     , kl_div_tr = proxy_kl_div
                     , precision_tt = forest_precision
                     , stability_tt	= forest_stability
                     , recall_tt = forest_recall
                     , f1_tt = forest_f1
                     , cc_tt = proxy_cc # doesn't change proxy / forest
                     , ci_tt = proxy_ci
                     , ncc_tt = proxy_ncc
                     , nci_tt = proxy_nci
                     , npv_tt = forest_npv
                     , accuracy_tt = forest_accu
                     , lift_tt = forest_lift
                     , coverage_tt = coverage
                     , xcoverage_tt = xcoverage
                     , kl_div_tt = forest_kl_div
                     , elapsed_time = rep(tpe, n_test))
          , file = paste0(project_dir
                          , datasetnames[i]
                          , pathsep
                          , algorithm
                          , "_rnst_"
                          , random_states[r]
                          , ".csv")
)

# collect summary results
f_perf <- mean(forest_label == ds_container$y_test)
p_perf <- mean(benchmark$model_accurate)
fid <- mean(benchmark$label == as.numeric(forest_label))

# save results to file before exiting loop
write.csv(data.frame(dataset_name = datasetnames[i]
          , algorithm = algorithm
          , n_instances = n_test
          , n_rules = benchmark$unique_rules
          , n_rules_used = benchmark$n_rules_used
          , median_rule_cascade = benchmark$median_rule_cascade
          , mean_rule_cascade = benchmark$mean_rule_cascade
          , sd_rule_cascade = benchmark$sd_rule_cascade
          , mean_rulelen = benchmark$mean_rulelen
          , sd_rulelen = benchmark$sd_rulelen
          , begin_time = as.character(benchmark$begin_time)
          , completion_time = as.character(benchmark$completion_time)
          , forest_performance = f_perf
          , sd_forest_performance = sqrt((f_perf/(1-f_perf))/length(forest_label))
          , model_kappa = cohen.kappa(table(forest_label
                                            , factor(ds_container$y_test)))$kappa
          , proxy_performance = p_perf
          , sd_proxy_performance = sqrt((p_perf/(1-p_perf))/length(benchmark$model_accurate))
          , proxy_kappa = cohen.kappa(table(benchmark$label
                                             , forest_label))$kappa
          , fidelity = fid
          , sd_fidelity = sqrt((fid/(1-fid))/n_test)
  ), file = paste0(project_dir
                   , datasetnames[i]
                   , pathsep
                   , algorithm
                   , "_rnst_"
                   , random_states[r]
                   , "_summary.csv")
)