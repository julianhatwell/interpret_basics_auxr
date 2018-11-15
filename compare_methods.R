library(randomForest)
library(jsonlite)
source("compare_methods_utils.R")

# algorithm <- "inTrees"
algorithm <- "BRL"

random_states <- 123:(123 + 0)

# importing data
data_dir <- "C:\\Users\\id126493\\Documents\\GitHub\\explain_te\\CHIRPS\\datafiles\\"
project_dir <- "V:\\whiteboxing\\"

class_cols <- c("income", "y"
                , "acceptability"
                , "NSP", "A16"
                , "rating", "loan_status"
                , "decision", "recid")

data_files <- c("adult_small_samp.csv.gz", "bankmark_samp.csv.gz"
, "car.csv.gz", "cardio.csv.gz", "credit.csv.gz"
, "german.csv.gz" , "lending_tiny_samp.csv.gz"
, "nursery_samp.csv.gz", "rcdv_samp.csv.gz")

datasets <- sapply(data_files, get_dataset_name)
output_dirs <- paste0(project_dir, datasets, "\\")

results_nrows <- length(random_states) * length(datasets)

dataset <- character(results_nrows)
n_instances <- integer(results_nrows)
random_state <- integer(results_nrows)
n_rules <- integer(results_nrows)
faith <- numeric(results_nrows)
sd_faith <- numeric(results_nrows)
mean_rule_cascade <- numeric(results_nrows)
mean_rulelen <- numeric(results_nrows)
mean_pred_acc <- numeric(results_nrows)
sd_rule_cascade <- numeric(results_nrows)
sd_rulelen <- numeric(results_nrows)
sd_pred_acc <- numeric(results_nrows)

mean_coverage <- numeric(results_nrows)
sd_coverage <- numeric(results_nrows)
mean_xcoverage <- numeric(results_nrows)
sd_xcoverage <- numeric(results_nrows)
mean_learner_precision <- numeric(results_nrows)
sd_learner_precision <- numeric(results_nrows)
mean_learner_stability <- numeric(results_nrows)
sd_learner_stability <- numeric(results_nrows)
mean_learner_recall <- numeric(results_nrows)
sd_learner_recall <- numeric(results_nrows)
mean_learner_f1 <- numeric(results_nrows)
sd_learner_f1 <- numeric(results_nrows)
mean_learner_accu <- numeric(results_nrows)
sd_learner_accu <- numeric(results_nrows)
mean_learner_lift <- numeric(results_nrows)
sd_learner_lift <- numeric(results_nrows)
mean_forest_precision <- numeric(results_nrows)
sd_forest_precision <- numeric(results_nrows)
mean_forest_stability <- numeric(results_nrows)
sd_forest_stability <- numeric(results_nrows)
mean_forest_recall <- numeric(results_nrows)
sd_forest_recall <- numeric(results_nrows)
mean_forest_f1 <- numeric(results_nrows)
sd_forest_f1 <- numeric(results_nrows)
mean_forest_accu <- numeric(results_nrows)
sd_forest_accu <- numeric(results_nrows)
mean_forest_lift <- numeric(results_nrows)
sd_forest_lift <- numeric(results_nrows)

for (i in seq_along(data_files)) {
  dat <- read.csv(gzfile(paste0(data_dir, data_files[i])))
  train_idx <- read.csv(paste0(output_dirs[i], "train_index.csv"), header = FALSE)$V1
  test_idx <- read.csv(paste0(output_dirs[i], "test_index.csv"), header = FALSE)$V1
  dat_train <- dat[train_idx, ]
  dat_test <- dat[test_idx, ]
  
  n_test <- length(test_idx)
  
  ds_container <- list(
    X_train = dat_train[, names(dat) != class_cols[i]],
    y_train = dat_train[, class_cols[i]],
    X_test = dat_test[, names(dat) != class_cols[i]],
    y_test = dat_test[, class_cols[i]]
  )
  
  classes <- as.character(unique(ds_container$y_train))
  
  fmla <- as.formula(paste(class_cols[i], "~ ."))
  
  ntree <- fromJSON(readLines(file(paste0(
    output_dirs[i]
    , "best_params_rndst_"
    , random_states[i]
    , ".json"))))$n_estimators
  for (r in seq_along(random_states)) {
    
    set.seed(random_states[r])
    rf <- randomForest(fmla, data=dat_train, ntree=ntree)
    forest_label <- which_class(as.character(predict(rf
                                , newdata = ds_container$X_test)))
    
    if (algorithm == "inTrees") {
      benchmark <- inTrees_benchmark(forest=rf, ds_container=ds_container)
    } else {
      benchmark <- sbrl_benchmark(ds_container=ds_container)
    }
    
    forest_vote_share <- numeric(n_test)
    prior <- numeric(n_test)
    coverage <- numeric(n_test)
    xcoverage <- numeric(n_test)
    learner_precision <- numeric(n_test)
    learner_stability <- numeric(n_test)
    learner_counts <- integer(n_test)
    learner_recall <- numeric(n_test)
    learner_f1 <- numeric(n_test)
    learner_accu <- numeric(n_test)
    learner_lift <- numeric(n_test)
    forest_precision <- numeric(n_test)
    forest_stability <- numeric(n_test)
    forest_counts <- integer(n_test)
    forest_recall <- numeric(n_test)
    forest_f1 <- numeric(n_test)
    forest_accu <- numeric(n_test)
    forest_lift <- numeric(n_test)
    for (j in 1:n_test) {
      rule <- benchmark$rule_idx[j]
      covered_instances <- benchmark$rule_idx == rule
      covered_instances <- covered_instances[-j] # drop out current instance
      loo_true_labels <- as.character(ds_container$y_test)[-j]
      
      instance_results <- evaluate(prior_labels = loo_true_labels
               , post_idx = covered_instances
               , classes = classes)
      
      forest_vote_share[j] <-
        mean(predict(rf, newdata=ds_container$X_test[1, ]
                     , predict.all = TRUE)$individual ==
               classes[forest_label[j]])
      prior[j] <- instance_results[["prior"]][forest_label[j]]
      coverage[j] <- instance_results[["coverage"]][1]
      xcoverage[j] <- instance_results[["xcoverage"]][1]
      learner_precision[j] <- instance_results[["posterior"]][benchmark$label[j]]
      learner_stability[j] <- instance_results[["stability"]][benchmark$label[j]]
      learner_recall[j] <- instance_results[["recall"]][benchmark$label[j]]
      learner_f1[j] <- instance_results[["f1"]][benchmark$label[j]]
      learner_accu[j] <- instance_results[["accu"]][benchmark$label[j]]
      learner_lift[j] <- instance_results[["lift"]][benchmark$label[j]]
      forest_precision[j] <- instance_results[["posterior"]][forest_label[j]]
      forest_stability[j] <- instance_results[["stability"]][forest_label[j]]
      forest_recall[j] <- instance_results[["recall"]][forest_label[j]]
      forest_f1[j] <- instance_results[["f1"]][forest_label[j]]
      forest_accu[j] <- instance_results[["accu"]][forest_label[j]]
      forest_lift[j] <- instance_results[["lift"]][forest_label[j]]
    }
    
    # collect results
    # save each run for comparable analysis
    rl_ln <- sapply(gregexpr("&", benchmark$model[benchmark$rule_idx,4]), function(x) {length(x)[[1]] + 1})
    write.csv(data.frame(dataset_name = rep(datasets[i], n_test)
                         , instance_id = test_idx
                         , algorithm = rep(algorithm, n_test)
                         , pretty_rule = benchmark$model[benchmark$rule_idx, 4]
                         , rule_length = rl_ln
                         , pred_class = forest_label
                         , pred_class_label = classes[forest_label]
                         , target_class = benchmark$label
                         , target_class_label = classes[benchmark$label]
                         , forest_vote_share = forest_vote_share
                         , prior = prior
                         , precision_tr = learner_precision
                         , stability_tr	= learner_stability
                         , recall_tr = learner_recall
                         , f1_tr = learner_f1
                         , accuracy_tr = learner_accu
                         , lift_tr = learner_lift
                         , precision_tt = forest_precision
                         , stability_tt	= forest_stability
                         , recall_tt = forest_recall
                         , f1_tt = forest_f1
                         , accuracy_tt = forest_accu
                         , lift_tt = forest_lift)
              , file = paste0(output_dirs[i], algorithm, "_rndst_", random_states[r], ".csv")
    )
    
    
    this_run <- length(random_states) * (i - 1) + r
    dataset[this_run] <- datasets[i]
    n_instances[this_run] <- n_test
    random_state[this_run] <- random_states[r]
    n_rules[this_run] <- nrow(learner)
    fth <- mean(learner_label == forest_label)
    faith[this_run] <- fth
    sd_faith[this_run] <- (fth/(1-fth))/n_test
    mean_rule_cascade[this_run] <- mean(rule_idx)
    sd_rule_cascade[this_run] <- sd(rule_idx)
    mean_rulelen[this_run] <- mean(rl_ln)
    sd_rulelen[this_run] <- sd(rl_ln)
    mean_coverage[this_run] <- mean(coverage)
    sd_coverage[this_run] <- sd(coverage)
    mean_xcoverage[this_run] <- mean(xcoverage)
    sd_xcoverage[this_run] <- sd(xcoverage)
    mean_learner_precision[this_run] <- mean(learner_precision)
    sd_learner_precision[this_run] <- sd(learner_precision)
    mean_learner_stability[this_run] <- mean(learner_stability)
    sd_learner_stability[this_run] <- sd(learner_stability)
    mean_learner_recall[this_run] <- mean(learner_recall)
    sd_learner_recall[this_run] <- sd(learner_recall)
    mean_learner_f1[this_run] <- mean(learner_f1)
    sd_learner_f1[this_run] <- sd(learner_f1)
    mean_learner_accu[this_run] <- mean(learner_accu)
    sd_learner_accu[this_run] <- sd(learner_accu)
    mean_learner_lift[this_run] <- mean(learner_lift)
    sd_learner_lift[this_run] <- sd(learner_lift)
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
  }
}

grand_results <- data.frame(dataset = dataset
                      , n_instances = n_instances
                      , random_state = random_state
                      , faith = faith
                      , sd_faith = sd_faith
                      , mean_rule_cascade = mean_rule_cascade
                      , sd_rule_cascade = sd_rule_cascade
                      , mean_rulelen = mean_rulelen
                      , sd_rulelen = sd_rulelen
                      , mean_precision = mean_precision
                      , sd_precision = sd_precision
                      , mean_coverage = mean_coverage
                      , sd_coverage = sd_coverage
                      , mean_xcoverage = mean_xcoverage
                      , sd_xcoverage = sd_xcoverage
                      , mean_learner_precision = mean_learner_precision
                      , sd_learner_precision = sd_learner_precision
                      , mean_learner_stability = mean_learner_stability
                      , sd_learner_stability = sd_learner_stability
                      , mean_learner_recall = mean_learner_recall
                      , sd_learner_recall = sd_learner_recall
                      , mean_learner_f1 = mean_learner_f1
                      , sd_learner_f1 = sd_learner_f1
                      , mean_learner_accu = mean_learner_accu
                      , sd_learner_accu = sd_learner_accu
                      , mean_learner_lift = mean_learner_lift
                      , sd_learner_lift = sd_learner_lift
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

