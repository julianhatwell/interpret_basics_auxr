library(randomForest)
library(jsonlite)
random_states <- 123:(123 + 30)
source("compare_methods_utils.R")

algorithm <- "inTrees"
# algorithm <- "BRL"

for (i in seq_along(data_files)) {
  
  data_prep(i)
  
  for (r in seq_along(random_states)) {
    
    set.seed(random_states[r])
    randfor <- randomForest(fmla, data=dat_train, ntree=ntree)
    forest_label <- which_class(as.character(predict(randfor
                                                     , newdata = ds_container$X_test)))   
    
    if (algorithm == "inTrees") {
      benchmark <- inTrees_benchmark(forest=randfor, ds_container=ds_container)
    } else {
      benchmark <- sbrl_benchmark(ds_container=ds_container, classes)
    }
    
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
        mean(predict(randfor, newdata=ds_container$X_test[1, ]
                     , predict.all = TRUE)$individual ==
               classes[forest_label[j]])
      prior[j] <- instance_results[["prior"]][forest_label[j]]
      coverage[j] <- instance_results[["coverage"]][1]
      xcoverage[j] <- instance_results[["xcoverage"]][1]
      proxy_precision[j] <- instance_results[["posterior"]][benchmark$label[j]]
      proxy_stability[j] <- instance_results[["stability"]][benchmark$label[j]]
      proxy_recall[j] <- instance_results[["recall"]][benchmark$label[j]]
      proxy_f1[j] <- instance_results[["f1"]][benchmark$label[j]]
      proxy_accu[j] <- instance_results[["accu"]][benchmark$label[j]]
      proxy_lift[j] <- instance_results[["lift"]][benchmark$label[j]]
      proxy_kl_div[j] <- entropy_corrected(instance_results[["posterior"]], instance_results[["prior"]])
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
      forest_lift[j] <- penalise_bad_prediction(forest_label[j]
                                                , benchmark$label[j]
                                                , instance_results[["lift"]][forest_label[j]])
      forest_kl_div[j] <- entropy_corrected(instance_results[["posterior"]], instance_results[["prior"]])
    }
    
    # collect results
    # save each run for comparable analysis
    write.csv(data.frame(dataset_name = rep(datasets[i], n_test)
                         , instance_id = test_idx
                         , algorithm = rep(algorithm, n_test)
                         , pretty_rule = benchmark$rule
                         , rule_length = benchmark$rl_ln
                         , pred_class = forest_label - 1 # conform with Python zero base
                         , pred_class_label = classes[forest_label]
                         , target_class = benchmark$label - 1  # conform with Python zero base
                         , target_class_label = classes[benchmark$label]
                         , forest_vote_share = forest_vote_share
                         , prior = prior
                         , precision_tr = proxy_precision
                         , stability_tr	= proxy_stability
                         , recall_tr = proxy_recall
                         , f1_tr = proxy_f1
                         , accuracy_tr = proxy_accu
                         , lift_tr = proxy_lift
                         , kl_div_tr = proxy_kl_div
                         , precision_tt = forest_precision
                         , stability_tt	= forest_stability
                         , recall_tt = forest_recall
                         , f1_tt = forest_f1
                         , accuracy_tt = forest_accu
                         , lift_tt = forest_lift
                         , kl_div_tt = forest_kl_div)
              , file = paste0(output_dirs[i], algorithm, "_rndst_", random_states[r], ".csv")
    )
    
    this_run <- length(random_states) * (i - 1) + r
    dataset[this_run] <- datasets[i]
    n_instances[this_run] <- n_test
    random_state[this_run] <- random_states[r]
    n_rules[this_run] <- benchmark$unique_rules
    n_rules_used[this_run] <- length(unique(benchmark$rule_idx))
    f_perf <- mean(forest_label == as.numeric(ds_container$y_test))
    forest_performance[this_run] <- f_perf
    sd_forest_performance[this_run] <- (f_perf/(1-f_perf))/length(forest_label)
    p_perf <- mean(benchmark$model_accurate)
    proxy_performance[this_run] <- p_perf
    sd_proxy_performance[this_run] <- (p_perf/(1-p_perf))/length(benchmark$model_accurate) # binomial sd
    fid <- mean(benchmark$label == forest_label)
    fidelity[this_run] <- fid
    sd_fidelity[this_run] <- (fid/(1-fid))/n_test # binomial sd
    mean_rule_cascade[this_run] <- mean(benchmark$rule_idx)
    sd_rule_cascade[this_run] <- sd(benchmark$rule_idx)
    mean_rulelen[this_run] <- mean(benchmark$rl_ln)
    sd_rulelen[this_run] <- sd(benchmark$rl_ln)
    mean_coverage[this_run] <- mean(coverage)
    sd_coverage[this_run] <- sd(coverage)
    mean_xcoverage[this_run] <- mean(xcoverage)
    sd_xcoverage[this_run] <- sd(xcoverage)
    mean_proxy_precision[this_run] <- mean(proxy_precision)
    sd_proxy_precision[this_run] <- sd(proxy_precision)
    mean_proxy_stability[this_run] <- mean(proxy_stability)
    sd_proxy_stability[this_run] <- sd(proxy_stability)
    mean_proxy_recall[this_run] <- mean(proxy_recall)
    sd_proxy_recall[this_run] <- sd(proxy_recall)
    mean_proxy_f1[this_run] <- mean(proxy_f1)
    sd_proxy_f1[this_run] <- sd(proxy_f1)
    mean_proxy_accu[this_run] <- mean(proxy_accu)
    sd_proxy_accu[this_run] <- sd(proxy_accu)
    mean_proxy_lift[this_run] <- mean(proxy_lift)
    sd_proxy_lift[this_run] <- sd(proxy_lift)
    mean_forest_precision[this_run] <- mean(forest_precision)
    sd_forest_precision[this_run] <- sd(forest_precision)
    mean_forest_stability[this_run] <- mean(forest_stability)
    sd_forest_stability[this_run] <- sd(forest_stability)
    mean_forest_recall[this_run] <- mean(forest_recall)
    sd_forest_recall[this_run] <- sd(forest_recall)
    mean_forest_f1[this_run] <- mean(forest_f1)
    sd_forest_f1[this_run] <- sd(forest_f1)
    mean_forest_accu[this_run] <- mean(forest_accu)
    sd_forest_accu[this_run] <- sd(forest_accu)
    mean_forest_lift[this_run] <- mean(forest_lift)
    sd_forest_lift[this_run] <- sd(forest_lift)
    
    print(c(datasets[i], random_states[r]))
  }
}

grand_results <- data.frame(dataset = dataset
                            , n_instances = n_instances
                            , random_state = random_state
                            , forest_performance = forest_performance
                            , sd_forest_performance = sd_forest_performance
                            , proxy_performance = proxy_performance
                            , sd_proxy_performance = sd_proxy_performance
                            , fidelity = fidelity
                            , sd_fidelity = sd_fidelity
                            , mean_rule_cascade = mean_rule_cascade
                            , sd_rule_cascade = sd_rule_cascade
                            , mean_rulelen = mean_rulelen
                            , sd_rulelen = sd_rulelen
                            , mean_coverage = mean_coverage
                            , sd_coverage = sd_coverage
                            , mean_xcoverage = mean_xcoverage
                            , sd_xcoverage = sd_xcoverage
                            , mean_proxy_precision = mean_proxy_precision
                            , sd_proxy_precision = sd_proxy_precision
                            , mean_proxy_stability = mean_proxy_stability
                            , sd_proxy_stability = sd_proxy_stability
                            , mean_proxy_recall = mean_proxy_recall
                            , sd_proxy_recall = sd_proxy_recall
                            , mean_proxy_f1 = mean_proxy_f1
                            , sd_proxy_f1 = sd_proxy_f1
                            , mean_proxy_accu = mean_proxy_accu
                            , sd_proxy_accu = sd_proxy_accu
                            , mean_proxy_lift = mean_proxy_lift
                            , sd_proxy_lift = sd_proxy_lift
                            , mean_forest_precision = mean_forest_precision
                            , sd_forest_precision = sd_forest_precision
                            , mean_forest_stability = mean_forest_stability
                            , sd_forest_stability = sd_forest_stability
                            , mean_forest_recall = mean_forest_recall
                            , sd_forest_recall = sd_forest_recall
                            , mean_forest_f1 = mean_forest_f1
                            , sd_forest_f1 = sd_forest_f1
                            , mean_forest_accu = mean_forest_accu
                            , sd_forest_accu = sd_forest_accu
                            , mean_forest_lift = mean_forest_lift
                            , sd_forest_lift = sd_forest_lift)

View(grand_results)
write.csv(grand_results
          , file=paste0(project_dir
                        , "grand_results_"
                        , algorithm
                        , "_"
                        , Sys.Date()))
