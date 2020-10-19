rm(list = ls())
library(psych)
library(ada)
library(randomForest)
library(gbm)
# library(foreach)
# library(parallel)
# library(doParallel)
random_states <- 123 # :152
source("compare_methods_utils.R")

max_tests <- 1000

# model <- "rf"
# model <- "samme"
model <- "samme.r"
# model <- "gbm"

# algorithm <- "inTrees"
algorithm <- "BRL"

results_nrows <- length(random_states) * length(datasetnames)

# choose a db
rnr <- 9

if (model=="rf") {
  # adult
  if (rnr == 1) ntree_divisor <- 5; inTrees_maxdepth <- 8; lambda <- 10; eta <- 1; rule_maxlen <- 8; nchain <- 10
  # bankmark
  if (rnr == 2) ntree_divisor <- 20; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 4; nchain <- 10
  # car
  if (rnr == 3) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # cardio
  if (rnr == 4) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 4; nchain <- 10
  # credit
  if (rnr == 5) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 8; nchain <- 10
  # german
  if (rnr == 6) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 4; nchain <- 10
  # lending_tiny_samp CANNOT BE DONE IN R
  
  # nursery
  if (rnr == 8) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # rcdv
  if (rnr == 9) ntree_divisor <- 10; inTrees_maxdepth <- 8; lambda <- 10; eta <- 1; rule_maxlen <- 4; nchain <- 10
}

if (model %in% c("samme", "samme.r")) {
  # adult
  if (rnr == 1) ntree_divisor <- 5; inTrees_maxdepth <- 8; lambda <- 10; eta <- 1; rule_maxlen <- 8; nchain <- 10
  # bankmark
  if (rnr == 2) ntree_divisor <- 20; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 4; nchain <- 10
  # car
  if (rnr == 3) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # cardio
  if (rnr == 4) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 4; nchain <- 10
  # credit
  if (rnr == 5) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 8; nchain <- 10
  # german
  if (rnr == 6) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 5; eta <- 1; rule_maxlen <- 4; nchain <- 10
  # lending_tiny_samp CANNOT BE DONE IN R
  
  # nursery
  if (rnr == 8) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # rcdv
  if (rnr == 9) ntree_divisor <- 10; inTrees_maxdepth <- 8; lambda <- 10; eta <- 1; rule_maxlen <- 4; nchain <- 10
}

if (model == "gbm") {
  # adult
  if (rnr == 1) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # bankmark
  if (rnr == 2) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # car
  if (rnr == 3) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # cardio
  if (rnr == 4) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # credit
  if (rnr == 5) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # german
  if (rnr == 6) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # lending_tiny_samp CANNOT BE DONE IN R
  
  # nursery
  if (rnr == 8) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 5; rule_maxlen <- 8; nchain <- 10
  # rcdv
  if (rnr == 9) ntree_divisor <- 1; inTrees_maxdepth <- 8; lambda <- 10; eta <- 1; rule_maxlen <- 4; nchain <- 10
}

# set counters
r <- ((rnr - 1) %/% length(datasetnames)) + 1
i <- ((rnr - 1) %% length(datasetnames)) + 1

# encapsulated set up
data_prep(i, max_tests = max_tests)

# classifier prep
fmla <- as.formula(paste(class_cols[i], "~ ."))

# build forest
set.seed(random_states[r])
if (model == "rf") {
  
  bestparams <- get_best_params()
  
  ntree <- bestparams$n_estimators
  
  maxnodes <- 2^bestparams$max_depth - 1
  
  forest <- randomForest(fmla
                         , data = dat_train
                         , ntree = ntree
                         , maxnodes = maxnodes)
  forest_label <- predict(forest
                          , newdata = ds_container$X_test)
  forest_preds <- as.numeric(forest_label)
  forest_vote_share <- apply(predict(forest, newdata=ds_container$X_test
                                     , type = "vote"), 1, max) # vote counts are normed by default
}
if (model %in% c("samme", "samme.r")) {
  type <- "discrete"
  if (model == "samme.r") type <- "real"
  
  bestparams <- get_best_params()
  ntree <- bestparams$n_estimators
  max_depth_pos <- regexpr("max_depth=", bestparams$base_estimator) + nchar("max_depth=")
  max_depth <- as.integer(substr(bestparams$base_estimator, max_depth_pos, max_depth_pos))
  
  if (length(classes) == 2) {
    forest <- ada(fmla, data = dat_train, iter = ntree, nu = 1
                      , loss = "exponential", type = type
                      , rpart.control(maxdepth = max_depth))
    forest_label <- predict(forest
                            , newdata = ds_container$X_test
                            , n.trees = ntree
                            , type = "vector")
    forest_label <- relevel(forest_label, classes[1]) # bug in ada? vector returns default alphabetically indexed factor, not the values used in the training.
    forest_preds <- as.integer(forest_label)
    
    tree_preds <- sapply(1:ntree, function(t) {
      tp <- predict(forest$model$trees[[t]], newdata = ds_container$X_test)
      return(apply(tp, 1, which.max))
    })
    forest_vote_share <- apply(t((forest_preds == tree_preds)) * forest$model$alpha, 2, sum) / sum(forest$model$alpha)
    
  } else {
    ds_container$y_train_bin <- mky(ds_container$y_train)
    names(ds_container$y_train_bin) <- classes
    ds_container$y_test_bin <- mky(ds_container$y_test_bin)
    names(ds_container$y_test_bin) <- classes
    forest_label <- matrix(nrow = length(ds_container$y_test)
                                         , ncol = length(classes))
    alphas <- matrix(nrow = ntree
                           , ncol = length(classes))
    tree_preds <- list()
    for(k in seq_along(classes)) {
      dat_train_bin <- dat_train
      dat_train_bin[[class_cols[i]]] <- ds_container$y_train_bin[[classes[k]]]
      if (datasetnames[i] == "nursery" && k == 1) { # trivial condition will cause an error on ada
        forest <- rpart(fmla, data = dat_train_bin, control = rpart.control(maxdepth = max_depth))
        forest_label[, k] <- ifelse(predict(forest
                                     , newdata = ds_container$X_test
                                     , n.trees = ntree
                                     , type = "vector") == 1, Inf, -Inf)
        alphas[, k] <- 1
        tree_preds[[k]] <- matrix(TRUE, nrow = length(forest_label[, k]), ncol = ntree)
      } else {
        forest <- ada(fmla, data = dat_train_bin
                      , iter = ntree, nu = 1
                      , loss = "exponential", type = type
                      , rpart.control(maxdepth = max_depth))
        forest_label[, k] <- predict(forest
                                     , newdata = ds_container$X_test
                                     , n.trees = ntree
                                     , type = "F")
        alphas[, k] <- forest$model$alpha
        tree_preds[[k]] <- sapply(1:ntree, function(t) {
          tp <- (predict(forest$model$trees[[t]], newdata = ds_container$X_test)[, "1"] > 0.5)
          return(tp)
        })
      }
    }
    forest_preds <- apply(forest_label, 1, which.max)
    forest_label <- factor(classes[forest_preds], levels = classes)
    forest_vote_share <- sapply(seq_along(forest_preds), function(j) sum(tree_preds[[forest_preds[j]]][j, ] * alphas[, forest_preds[j]]) / sum(alphas[, forest_preds[j]]))
  }
}
if (model == "gbm") {
  bestparams <- get_best_params()
  
  if (!is.na(positive_classes[i])) {
    dat_train[, class_cols[i]] <- ifelse(dat_train[, class_cols[i]] == positive_classes[i], 1, 0)
    dat_test[, class_cols[i]] <- ifelse(dat_test[, class_cols[i]] == positive_classes[i], 1, 0)
  }

  forest <- gbm(fmla, data = dat_train
                , n.trees = bestparams$n_estimators
                , interaction.depth = bestparams$max_depth
                , shrinkage = bestparams$learning_rate
                , bag.fraction = bestparams$subsample)
  
  forest_label <- predict(forest
                          , newdata = ds_container$X_test
                          , n.trees = bestparams$n_estimators
                          , type = "response")

  
  if (!is.na(positive_classes[i])) { # binary
    forest_label <- factor(ifelse(forest_label > 0.5, positive_classes[i], negative_classes[i]))
    forest_label <- relevel(forest_label, positive_classes[i])
    forest_preds <- as.numeric(forest_label)
    tree_preds <- matrix(nrow = length(forest_preds), ncol = bestparams$n_estimators + 1)
    for (t in 1:(bestparams$n_estimators + 1)) {
      single.tree = TRUE
      if (t - 1 == 0) single.tree = FALSE
      tree_preds[, t] <- predict(forest, newdata=ds_container$X_test
                                 , n.trees = t - 1
                                 , single.tree = single.tree
                                 , type = "link")
    }
    sign_class <- ifelse(forest_label %in% positive_classes, 1, -1)
    sign_tree <- sign(tree_preds) == sign_class
    agree_sign <- abs(sapply(1:length(sign_class), function(sc) {
      sum(tree_preds[sc, sign_tree[sc, ]])
    }))
    disagree_sign <- abs(sapply(1:length(sign_class), function(sc) {
      sum(tree_preds[sc, !sign_tree[sc, ]])
    }))
    forest_vote_share <- agree_sign / (agree_sign + disagree_sign)
    
  } else { # multi-nomial
    forest_preds <- apply(forest_label, 1, which.max)
    forest_label <- factor(classes[forest_preds], levels = classes)
    forest_vote_share <- rep(NA, n_test)
    # predict gbm is buggy or certainly does not work as expected for multi-class.
    # can't get forest_vote_share at all, the numbers don't make any sense
  }
}

# run method and benchmark
if (algorithm == "inTrees") {
  benchmark <- inTrees_benchmark(forest = forest
                                 , ds_container = ds_container
                                 , ntree = bestparams$n_estimators / ntree_divisor
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

# collect     results
tpe <- time_per_explanation(benchmark$begin_time
                            , benchmark$completion_time
                            , n_test)
results_init(n_test)
for (j in 1:n_test) {
  rule <- benchmark$rule[j]
  
  covered_instances <- apply_rule(rule, ds_container$X_test, algorithm)
  covered_instances <- covered_instances[-j] # drop out current instance
  loo_true_labels <- ds_container$y_test[-j]

  instance_results <- evaluate(prior_labels = loo_true_labels
                               , post_idx = covered_instances
                               , classes = classes)

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
               , pred_class = forest_preds[1:n_test] - 1 # conform with Python zero base
                     , pred_class_label = forest_label[1:n_test]
               , target_class = benchmark$label[1:n_test] - 1  # conform with Python zero base
               , target_class_label = classes[benchmark$label][1:n_test]
               , forest_vote_share = forest_vote_share[1:n_test]
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

make_sym <- function(x) {
  dm <- dim(x)
  d <- diff(dm)
  if(d == 0) return(x)
  if(d > 0) { # add rows
    x <- rbind(x, rep(rep(0, dm[2]), d))
    return(x)
  } else { # add columns
    x <- cbind(x, rep(rep(0, dm[1]), d))
    return(x)
  }
}

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
    , proxy_kappa = cohen.kappa(make_sym(table(benchmark$label
                                       , forest_label)))$kappa
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
