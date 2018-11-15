library(sbrl)
library(inTrees)
library(rattle)

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

get_dataset_name <- function(x) {
  sub(".csv.gz", "", x)
}

whichRule <- function(learner, X) {
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

which_class <- function(preds) {
  sapply(preds, function(p) {
    which(classes == p)
  })
}

inTrees_benchmark <- function(forest, ds_container) {
  
  treeList <- RF2List(forest) # transform rf object to an inTrees" format
  extract <- extractRules(treeList, ds_container$X_train, ntree = ntree, maxdepth = 1000) # R-executable conditions
  ruleMetric <- getRuleMetric(extract, ds_container$X_train, ds_container$y_train) # get rule metrics
  ruleMetric <- pruneRule(ruleMetric, ds_container$X_train, ds_container$y_train)
  ruleMetric <- selectRuleRRF(ruleMetric, ds_container$X_train, ds_container$y_train)
  learner <- buildLearner(ruleMetric, ds_container$X_train, ds_container$y_train)
  
  learner_preds <- applyLearner(learner, ds_container$X_test)
  learner_label <- which_class(applyLearner(learner, ds_container$X_test))
  
  rule_idx <- whichRule(learner, ds_container$X_test)
  return(list(label=learner_label
              , rule_idx=rule_idx
              , model=learner
              , model_type="inTrees")) # which rule applies to which instance
  
}

sbrl_benchmark <- function(ds_container) {
  # transform for sbrl
  label <- factor(ifelse(ds_container$y_train == "<=50K", 0, 1))
  tdata <- ds_container$X_train
  # discretise where necessary
  for(cl in names(tdata)) {
    if(class(tdata[[cl]]) != "factor") {
      tdata[[cl]] <- binning(tdata[[cl]], method="quantile", ordered = FALSE)
    }
  }
  tdata <- cbind(tdata, label)
  model <- sbrl(tdata, rule_minlen=1, rule_maxlen=1000, 
                minsupport_pos=0.02, minsupport_neg=0.02, 
                lambda=10.0, eta=2.5, nchain=10
                )
  return(model)
}
