library(psych)
library(ada)
library(randomForest)
library(foreach)
library(doParallel)
n_cores <- 1 # detectCores() - 2
random_states <- 123 # :152
source("compare_methods_utils.R")

max_tests <- 1000

model <- "rf"
# model <- "ada"
algorithm <- "inTrees"
# algorithm <- "BRL"

results_nrows <- length(random_states) * length(datasetnames)

cl <- makeCluster(n_cores)
registerDoParallel(cl)

gres <- foreach(rnr = 1:results_nrows
        , .packages = c("jsonlite"
                        , "randomForest"
                        , "ada"
                        , "rattle"
                        , "inTrees"
                        , "sbrl"))  %dopar% {

for (rnr in 1:results_nrows) {
          
          # set counters
          r <- ((rnr - 1) %/% length(datasetnames)) + 1
          i <- ((rnr - 1) %% length(datasetnames)) + 1
          
          # encapsulated set up
          data_prep(i, max_tests = 1000)
          
          # classifier prep
          fmla <- as.formula(paste(class_cols[i], "~ ."))
          
          ntree <- fromJSON(readLines(file(paste0(
            resfilesdirs[i]
            , model
            , "_best_params_rnst_"
            , random_states[r]
            , ".json"))))$n_estimators
          
          if (model == "rf") {
            maxnodes <- 2^fromJSON(readLines(file(paste0(
              resfilesdirs[i]
              , model
              , "_best_params_rnst_"
              , random_states[r]
              , ".json"))))$max_depth
          }
          
          if (model == "ada") {
            base_estimator <- fromJSON(readLines(file(paste0(
              resfilesdirs[i]
              , model
              , "_best_params_rnst_"
              , random_states[r]
              , ".json"))))$base_estimator
            max_depth_pos <- regexpr("max_depth=", base_estimator) + nchar("max_depth=")
            max_depth <- as.integer(substr(base_estimator, max_depth_pos, max_depth_pos))
            }
          
          # build forest
          set.seed(random_states[r])
          if (model == "rf") forest <- randomForest(fmla, data = dat_train
                                                    , ntree = ntree
                                                    , maxnodes = maxnodes)
          if (model == "ada") forest <- ada(fmla, data = dat_train, iter = ntree
                                            , rpart.control(maxdepth = max_depth))
          forest_label <- which_class(as.character(predict(forest
                                                           , newdata = ds_container$X_test)))
          # run method and benchmark
          if (algorithm == "inTrees") {
            benchmark <- inTrees_benchmark(forest = forest
                                           , ds_container = ds_container
                                           , ntree = ntree
                                           , maxdepth = 1000
            )
          } else {
            benchmark <- sbrl_benchmark(ds_container=ds_container
                                        , classes = classes)
          }
          
          # collect results
          tpe <- time_per_explanation(benchmark$begin_time
                                      , benchmark$completion_time
                                      , n_test)
          results_init(n_test)
          for (j in 1:n_test) {
            rule <- benchmark$rule_idx[j]
            
            covered_instances <- benchmark$rule_idx == rule
            covered_instances <- covered_instances[-j] # drop out current instance
            loo_true_labels <- as.character(ds_container$y_test)[-j]

            instance_results <- evaluate(prior_labels = loo_true_labels
                                         , post_idx = covered_instances
                                         , classes = classes)

            forest_vote_share[j] <-
              mean(predict(forest, newdata=ds_container$X_test[1, ]
                           , predict.all = TRUE)$individual ==
                     classes[forest_label[j]])
            prior[j] <- instance_results[["prior"]][forest_label[j]]
            coverage[j] <- instance_results[["coverage"]][1]
            xcoverage[j] <- instance_results[["xcoverage"]][1]
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
            forest_precision[j] <- penalise_bad_prediction(forest_label[j]
                                                           , benchmark$label[j]
                                                           , instance_results[["posterior"]][forest_label[j]])
            forest_stability[j] <- penalise_bad_prediction(forest_label[j]
                                                           , benchmark$label[j]
                                                           , instance_results[["stability"]][forest_label[j]])
            forest_recall[j] <- penalise_bad_prediction(forest_label[j]
                                                        , benchmark$label[j]
                                                        , instance_results[["recall"]][forest_label[j]])
            forest_f1[j] <- penalise_bad_prediction(forest_label[j]
                                                    , benchmark$label[j]
                                                    , instance_results[["f1"]][forest_label[j]])
            forest_accu[j] <- penalise_bad_prediction(forest_label[j]
                                                      , benchmark$label[j]
                                                      , instance_results[["accu"]][forest_label[j]])
            forest_kl_div[j] <- penalise_bad_prediction(forest_label[j]
                                                      , benchmark$label[j]
                                                      , instance_results[["kl_div"]][forest_label[j]])
            forest_npv[j] <- penalise_bad_prediction(forest_label[j]
                                                        , benchmark$label[j]
                                                        , instance_results[["npv"]][forest_label[j]])
            forest_lift[j] <- ifelse(forest_precision[j] == 0, 1, proxy_lift[j])
          }
          
          # save results to file before exiting loop
          write.csv(data.frame(dataset_name = rep(datasetnames[i], n_test)
                               , instance_id = test_idx - 1 # conform with Python zero base
                               , algorithm = rep(algorithm, n_test)
                               , pretty_rule = benchmark$rule
                               , rule_length = benchmark$rl_ln
                               , pred_class = forest_label - 1 # conform with Python zero base
                               , pred_class_label = classes[forest_label]
                               , target_class = benchmark$label - 1  # conform with Python zero base
                               , target_class_label = classes[benchmark$label]
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
                    , file = paste0(resfilesdirs[i], algorithm, "_rnst_", random_states[r], ".csv")
          )
          
          # collect summary results
          f_perf <- mean(forest_label == as.numeric(ds_container$y_test))
          p_perf <- mean(benchmark$model_accurate)
          fid <- mean(benchmark$label == forest_label)
          
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
                    , model_kappa = cohen.kappa(table(factor(names(forest_label)
                                                       , levels = levels(factor(ds_container$y_test)))
                                                , factor(ds_container$y_test)))$kappa
                    , proxy_performance = p_perf
                    , sd_proxy_performance = sqrt((p_perf/(1-p_perf))/length(benchmark$model_accurate))
                    , proxy_kappa = cohen.kappa(table(factor(benchmark$label
                                                       , levels = levels(factor(forest_label)))
                                                , factor(forest_label)))$kappa
                    , fidelity = fid
                    , sd_fidelity = sqrt((fid/(1-fid))/n_test)
            ), file = paste0(resfilesdirs[i], algorithm, "_rnst_", random_states[r], "_summary.csv")
          )
          
          # don't return anything to foreach
          benchmark
}

} 
stopCluster(cl)
