library(sbrl)
library(inTrees)
library(rattle)
library(jsonlite)

source("data_files_mgmt.R")

# utility functions
p_count <- function(numvec) {
  tab <- table(numvec)
  counts <- as.vector(tab)
  return(list(
    labels = names(tab)
    , counts = counts
    , p_counts = counts / length(numvec)
    , s_counts = counts / (length(numvec) + 1)
  ))
  
}

# insert zeros for unrepresented classes
p_count_corrected <- function(arr, classes) {
  n_classes <- length(classes)
  p_counts <- p_count(arr)
  sc <- numeric(n_classes)
  pc <- numeric(n_classes)
  c <- integer(n_classes)
  
  for (i in seq_along(classes)) {
    if (classes[i] %in% p_counts[["labels"]]) {
      sc[i] <- p_counts[["s_counts"]][which(p_counts[["labels"]] == classes[i])]
      pc[i] <- p_counts[["p_counts"]][which(p_counts[["labels"]] == classes[i])]
      c[i] <- p_counts[["counts"]][which(p_counts[["labels"]] == classes[i])]
    }
  }
  
  return(list(labels = classes
              , counts = c
              , p_counts = pc
              , s_counts = sc))
}

entropy_corrected <- function(p, q) {
  n <- length(p)
  p_norm <- p/sum(p)
  q_norm <- q/sum(q)
  p_smooth <- runif(n)
  q_smooth <- runif(n)
  p_smoothed <- 0.01 * p_smooth + 0.99 * p_norm
  q_smoothed <- 0.01 * q_smooth + 0.99 * q_norm
  return(-sum(p_smoothed * log(q_smoothed/p_smoothed)))
}

evaluate <- function(prior_labels, post_idx, classes) {
  prior <- p_count_corrected(prior_labels, classes)
  
  coverage = mean(post_idx) # tp + fp / tp + fp + tn + fn  
  xcoverage = sum(post_idx) / (length(post_idx) + 1) # tp + fp / tp + fp + tn + fn + current instance
  
  p_counts = p_count_corrected(prior_labels[post_idx], classes)
  
  posterior = p_counts[["p_counts"]]
  stability = p_counts[["s_counts"]]
  counts = p_counts[["counts"]]
  labels = p_counts[["labels"]]
  
  recall = counts / prior[["counts"]] # TPR (recall) TP / (TP + FN)
  
  # F1
  p_corrected = ifelse(posterior > 0.0, posterior, 1.0) # to avoid div by zeros
  r_corrected = ifelse(recall > 0.0, recall, 1.0) # to avoid div by zeros
  f1 = 2 * ((posterior * recall) / (p_corrected + r_corrected))
  
  not_covered_counts = counts + sum(prior[["counts"]]) - prior[["counts"]] - (sum(counts) - counts)
  # accuracy = (TP + TN) / num_instances
  accu = not_covered_counts/sum(prior[["counts"]])
  
  # to avoid div by zeros
  pri_corrected = ifelse(prior[["p_counts"]] > 0.0, prior[["p_counts"]], 1.0)
  pos_corrected = ifelse(prior[["p_counts"]] > 0.0, posterior, 1.0)
  
  if (sum(counts) == 0) {
    rec_corrected <- numeric(length(pos_corrected))
    cov_corrected <- rep(1, length(pos_corrected))
  } else {
    rec_corrected <- counts / sum(counts)
    cov_corrected <- sum(counts) / sum(prior[["counts"]])
  }
  
  # lift = precis / (total_cover * prior)
  lift = pos_corrected / ( cov_corrected * pri_corrected )
  
  chisq <- chisq.test(rbind(counts, prior[["counts"]]))[["p.value"]]
  
  return(list(coverage = coverage
                    , xcoverage = xcoverage
                    , stability = stability
                    , prior = prior[["p_counts"]]
                    , posterior = p_counts[["p_counts"]]
                    , counts = counts
                    , labels = labels
                    , recall = recall
                    , f1 = f1
                    , accu = accu
                    , lift = lift
                    , chisq = chisq))
}

data_prep <- function(i) {
  dat <- read.csv(gzfile(paste0(data_dir, data_files[i])))
  # ensure y is a factor
  if (class(dat[, class_cols[i]]) != "factor") dat[, class_cols[i]] <- factor(dat[, class_cols[i]])
  classes <<- levels(dat[, class_cols[i]])
  
  # test train split according to indices exported from Python
  train_idx <<- read.csv(paste0(resfilesdirs[i], "train_index.csv"), header = FALSE)$V1 + 1
  test_idx <<- read.csv(paste0(resfilesdirs[i], "test_index.csv"), header = FALSE)$V1 + 1
  dat_train <<- dat[train_idx, ]
  dat_test <<- dat[test_idx, ]
  
  n_test <<- length(test_idx)
  
  ds_container <<- list(
    X_train = dat_train[, names(dat) != class_cols[i]],
    y_train = dat_train[, class_cols[i]],
    X_test = dat_test[, names(dat) != class_cols[i]],
    y_test = dat_test[, class_cols[i]]
  )
  
  fmla <<- as.formula(paste(class_cols[i], "~ ."))
  
  ntree <<- fromJSON(readLines(file(paste0(
    resfilesdirs[i]
    , "best_params_rnst_"
    , random_states[r]
    , ".json"))))$n_estimators
}

results_init <<- function(n_test) {
  forest_vote_share <<- numeric(n_test)
  prior <<- numeric(n_test)
  coverage <<- numeric(n_test)
  xcoverage <<- numeric(n_test)
  proxy_precision <<- numeric(n_test)
  proxy_stability <<- numeric(n_test)
  proxy_counts <<- integer(n_test)
  proxy_recall <<- numeric(n_test)
  proxy_f1 <<- numeric(n_test)
  proxy_accu <<- numeric(n_test)
  proxy_lift <<- numeric(n_test)
  proxy_kl_div <<- numeric(n_test)
  forest_precision <<- numeric(n_test)
  forest_stability <<- numeric(n_test)
  forest_counts <<- integer(n_test)
  forest_recall <<- numeric(n_test)
  forest_f1 <<- numeric(n_test)
  forest_accu <<- numeric(n_test)
  forest_lift <<- numeric(n_test)
  forest_kl_div <<- numeric(n_test)
}

whichRule_inTrees <- function(learner, X) {
  left_idx <- 1:nrow(X)
  rules_idx <- numeric(length(left_idx))
  for (i in 1:nrow(learner)) {
    match_idx <- eval(parse(text = paste("which(", learner[i, 
                                                           "condition"], ")")))
    match_idx <- intersect(left_idx, match_idx)
    rules_idx[match_idx] <- i
    left_idx <- setdiff(left_idx, match_idx)
  }
  return(rules_idx)
}

whichRule_sbrl <- function (model, tdata) 
{
  mat_data_feature <- get_data_feature_mat(tdata, model$featurenames)
  mat_data_rules <- mat_data_feature %*% model$mat_feature_rule
  mat_data_rules <- t(t(mat_data_rules) >= c(colSums(model$mat_feature_rule))) + 
    0
  nrules <- ncol(model$mat_feature_rule)
  nsamples <- nrow(tdata)
  mat_idx <- matrix(0, nrow = nsamples, ncol = nrules)
  for (i in 1:(nrow(model$rs) - 1)) {
    mat_idx[, model$rs$V1[i]] = i
  }
  mat_satisfy <- mat_data_rules * mat_idx
  mat_caps <- as.matrix(apply(mat_satisfy, 1, function(x) ifelse(!identical(x[x > 
                                                                                0], numeric(0)), min(x[x > 0]), NaN)))
  mat_caps[is.na(mat_caps)] = nrow(model$rs)
  
  return(model$rs$V1[mat_caps])
}

which_class <- function(preds) {
  sapply(preds, function(p) {
    which(classes == p)
  })
}

penalise_bad_prediction <- function(mc, tc, value) {
  return(ifelse(mc == tc, value, 0))
}

inTrees_benchmark <- function(forest, ds_container, ntree, maxdepth) {
  begin_time <- as.character(Sys.time())
  treeList <- RF2List(forest) # transform rf object to an inTrees" format
  extract <- extractRules(treeList, ds_container$X_train, ntree = ntree, maxdepth = maxdepth) # R-executable conditions
  ruleMetric <- getRuleMetric(extract, ds_container$X_train, ds_container$y_train) # get rule metrics
  ruleMetric <- pruneRule(ruleMetric, ds_container$X_train, ds_container$y_train)
  ruleMetric <- selectRuleRRF(ruleMetric, ds_container$X_train, ds_container$y_train)
  learner <- buildLearner(ruleMetric, ds_container$X_train, ds_container$y_train)
  
  learner_preds <- applyLearner(learner, ds_container$X_test)
  learner_label <- which_class(applyLearner(learner, ds_container$X_test))
  model_accurate <- ifelse(learner_label == as.numeric(ds_container$y_test), 1, 0)
  
  rule_idx <- whichRule_inTrees(learner, ds_container$X_test)
  rl_ln <- sapply(gregexpr("&", learner[rule_idx,4]), function(x) {length(x)[[1]] + 1})
  
  return(list(this_i = i
              , this_r = r
              , this_run = length(random_states) * (i - 1) + r
              , random_state = random_states[r]
              , datasetname = datasetnames[i]
              , label = learner_label
              , rule_idx = rule_idx # which rule applies to which instance
              , rl_ln = rl_ln
              , unique_rules = nrow(learner)
              , n_rules_used = length(unique(rule_idx))
              , rule = learner[rule_idx, 4]
              , mean_rule_cascade = mean(rule_idx)
              , sd_rule_cascade = sd(rule_idx)
              , mean_rulelen = mean(rl_ln)
              , sd_rulelen = sd(rl_ln)
              , model = learner
              , model_accurate = model_accurate
              , model_type="inTrees"
              , begin_time = begin_time
              , completion_time = as.character(Sys.time())
              ))
  
}

sbrl_data_prep <- function(dat) {
  for(cl in names(dat)) {
    if(class(dat[[cl]]) != "factor") {
      dat[[cl]] <- binning(dat[[cl]], method="quantile", ordered = FALSE)
    }
  }
  return(dat)
}

set_labels <- function(y_tr, zvl) {
  factor(ifelse(y_tr == zvl, 0, 1)) # error if not in this format  
}

sbrl_benchmark <- function(ds_container, classes) {
  begin_time <- as.character(Sys.time())
  # transform for sbrl
  if (datasetnames[i] == "adult_small_samp") zero_val_label <- "<=50K"
  if (datasetnames[i] == "bankmark_samp") zero_val_label <- "no"
  if (datasetnames[i] == "car") zero_val_label <- "acc"
  if (datasetnames[i] == "credit") zero_val_label <- "minus"
  if (datasetnames[i] == "german") zero_val_label <- "bad"
  if (datasetnames[i] == "lending_tiny_samp") zero_val_label <- "Charged Off"
  if (datasetnames[i] == "rcdv_samp") zero_val_label <- "N"
  
  train_label <- set_labels(ds_container$y_train, zero_val_label)
  test_label <- set_labels(ds_container$y_test, zero_val_label)
  
  train_data <- sbrl_data_prep(ds_container$X_train)
  train_data <- cbind(train_data, label = train_label)
  
  test_data <- sbrl_data_prep(ds_container$X_test)
  
  model <- sbrl(tdata=train_data, rule_minlen=1, rule_maxlen=5, 
                minsupport_pos=0.05, minsupport_neg=0.05,
                lambda=10.0, eta=2.5, nchain=10)
  
  sbrl_label <- ifelse(predict(model, tdata=test_data)$V1 > 0.5, 1, 2)
  model_accurate <- ifelse(sbrl_label == as.numeric(test_label), 1, 0)
  rules_plus_default <- c(model$rulenames, "{default}") # format of returned rules. The = sign will help counting 
  n_rules <- length(model$rulenames) + 1
  rule_idx <- whichRule_sbrl(model, test_data)
  rule_pos <- sapply(rule_idx, function(x) {which(model$rs$V1 == x)}) # create a positional index
  
  rule_idx <- ifelse(rule_idx == 0, n_rules, rule_idx) # makes rule extract easier
  rl_ln <- sapply(gregexpr("=", rules_plus_default[rule_idx]), function(x) {ifelse(x[1] == -1, 0, length(x)[1])})
  rule <- rules_plus_default[rule_idx]
  
  return(list(this_i = i
              , this_r = r
              , this_run = length(random_states) * (i - 1) + r
              , random_state = random_states[r]
              , datasetname = datasetnames[i]
              , label=sbrl_label
              , rule_idx = rule_pos
              , rl_ln = rl_ln
              , unique_rules = length(model$rs$V1)
              , n_rules_used = length(unique(rule_idx))
              , rule = rule
              , mean_rule_cascade = mean(rule_pos)
              , sd_rule_cascade = sd(rule_pos)
              , mean_rulelen = mean(rl_ln)
              , sd_rulelen = sd(rl_ln)
              , model = model
              , model_accurate = model_accurate
              , model_type="sbrl"
              , begin_time = begin_time
              , completion_time = as.character(Sys.time())
              ))
}
