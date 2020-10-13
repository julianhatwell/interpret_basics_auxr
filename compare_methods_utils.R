library(sbrl)
library(inTrees)
library(jsonlite)
library(rattle)
library(rlang)
library(dplyr)

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
  
  all_c <- length(prior_labels)
  prior <- p_count_corrected(prior_labels, classes)
  
  # basic results
  if (sum(post_idx) == 0) {
    # no coverage
    zeros <- rep(0, length(classes))
    counts <- zeros
    covered <- zeros
    ci <- zeros
    labels <- factor(character(0), levels = classes)
    posterior <- zeros
  } else {
    p_counts = p_count_corrected(prior_labels[post_idx], classes)
    counts = p_counts[["counts"]]
    covered <- sum(counts)
    ci <- covered - counts
    labels <- p_counts[["labels"]]
    posterior <- p_counts[["p_counts"]]
  }

  # coverage
  coverage <- covered / all_c # tp + fp / tp + fp + tn + fn  
  xcoverage <- (covered + 1) / (all_c + length(classes) + 1)  # tp + fp / tp + fp + tn + fn + current instance, laplace corrected
  
  # stab = tp / tp + fp + current instance. laplace corrected
  stability <- (counts + 1) / (sum(counts) + length(classes) + 1)
  
  # negative results
  np_counts <- p_count_corrected(prior_labels[!post_idx], classes)
  ncounts <- np_counts[["counts"]]
  nc <- all_c - covered
  nci <- sum(ncounts) - ncounts
  nposterior <- np_counts[["p_counts"]]
  # negative predictive value = tn / tn + fn
  npv <- nci / (nci + ci)
  
  chisq <- chisq.test(rbind(counts, prior[["counts"]]))[["p.value"]]
  kl_div <- entropy_corrected(posterior, prior[["p_counts"]])
  
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
  
  return(list(count_all = all_c
              , covered = covered
              , not_covered = nc
              , cc = counts
              , ci = ci
              , ncc = ncounts
              , nci = nci
              , coverage = coverage
              , xcoverage = xcoverage
              , npv = npv
              , stability = stability
              , prior = prior[["p_counts"]]
              , posterior = posterior
              , counts = counts
              , labels = labels
              , recall = recall
              , f1 = f1
              , accu = accu
              , lift = lift
              , chisq = chisq
              , kl_div = kl_div))
}

get_train_test_sizes <- function(i) {
  train_idx <- read.csv(paste0(resfilesdirs[i], "train_index.csv"), header = FALSE)$V1 + 1
  test_idx <- read.csv(paste0(resfilesdirs[i], "test_index.csv"), header = FALSE)$V1 + 1
  return(c(length(train_idx), length(test_idx)))
}

data_prep <- function(i, max_tests) {
  dat <- read.csv(gzfile(paste0(data_dir, data_files[i])))
  # ensure y is a factor
  if (class(dat[, class_cols[i]]) != "factor") dat[, class_cols[i]] <- factor(dat[, class_cols[i]])
  if (!(is.na(datasets_master$positive_classes[i]))) dat[, class_cols[i]] <- relevel(dat[, class_cols[i]], datasets_master$positive_classes[i])
  classes <<- levels(dat[, class_cols[i]])
  
  # test train split according to indices exported from Python
  train_idx <<- read.csv(paste0(project_dir, datasetnames[i], pathsep, "train_index.csv"), header = FALSE)$V1 + 1
  test_idx <<- read.csv(paste0(project_dir, datasetnames[i], pathsep, "test_index.csv"), header = FALSE)$V1 + 1
  dat_train <<- dat[train_idx, ]
  dat_test <<- dat[test_idx, ]
  
  n_test <<- min(length(test_idx), max_tests)
  
  ds_container <<- list(
    X_train = dat_train[, names(dat) != class_cols[i]],
    y_train = dat_train[, class_cols[i]],
    X_test = dat_test[, names(dat) != class_cols[i]],
    y_test = dat_test[, class_cols[i]]
  )
}

mky <- function(y) {
  lapply(classes, function(k, y) {
    factor(2 * (y == k) - 1, levels = c("1", "-1"))}, y = y)
}

results_init <<- function(n_test) {
  prior <<- numeric(n_test)
  coverage <<- numeric(n_test)
  xcoverage <<- numeric(n_test)
  proxy_precision <<- numeric(n_test)
  proxy_stability <<- numeric(n_test)
  proxy_counts <<- integer(n_test)
  proxy_recall <<- numeric(n_test)
  proxy_f1 <<- numeric(n_test)
  proxy_cc <<- numeric(n_test)
  proxy_ci <<- numeric(n_test)
  proxy_ncc <<- numeric(n_test)
  proxy_nci <<- numeric(n_test)
  proxy_npv <<- numeric(n_test)
  proxy_accu <<- numeric(n_test)
  proxy_lift <<- numeric(n_test)
  proxy_kl_div <<- numeric(n_test)
  forest_precision <<- numeric(n_test)
  forest_stability <<- numeric(n_test)
  forest_counts <<- integer(n_test)
  forest_recall <<- numeric(n_test)
  forest_f1 <<- numeric(n_test)
  forest_cc <<- numeric(n_test)
  forest_ci <<- numeric(n_test)
  forest_ncc <<- numeric(n_test)
  forest_nci <<- numeric(n_test)
  forest_npv <<- numeric(n_test)
  forest_accu <<- numeric(n_test)
  forest_lift <<- numeric(n_test)
  forest_kl_div <<- numeric(n_test)
}

intTrees_whichRule <- function(learner, X) {
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

sbrl_whichRule <- function (model, tdata) 
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

inTrees_benchmark <- function(forest, ds_container, ntree, maxdepth, model) {
  begin_time <- Sys.time()
  if(model != "gbm") {
    treeList <- RF2List(forest) # transform rf object to an inTrees" format
  } else {
    treeList <- GBM2List(forest, ds_container$X_train) # transform rf object to an inTrees" format
  }
  
  extract <- extractRules(treeList, ds_container$X_train, ntree = ntree, maxdepth = maxdepth) # R-executable conditions
  ruleMetric <- getRuleMetric(extract, ds_container$X_train, ds_container$y_train) # get rule metrics
  ruleMetric <- pruneRule(ruleMetric, ds_container$X_train, ds_container$y_train)
  ruleMetric <- selectRuleRRF(ruleMetric, ds_container$X_train, ds_container$y_train)
  learner <- buildLearner(ruleMetric, ds_container$X_train, ds_container$y_train)
  
  learner_preds <- applyLearner(learner, ds_container$X_test)
  learner_label <- which_class(applyLearner(learner, ds_container$X_test))
  model_accurate <- ifelse(learner_label == as.numeric(ds_container$y_test), 1, 0)
  
  rule_idx <- intTrees_whichRule(learner, ds_container$X_test)
  rl_ln <- sapply(gregexpr("&", learner[rule_idx,4])
                  , function(x) {
                    if (x == -1) {
                      return(0)
                    } else {
                      length(x)[[1]][1] + 1
                    }
                  })
  
  rl_ln <- rl_ln + sapply(gregexpr("%in%", learner[rule_idx,4])
                          , function(x) {
                            if (x == -1) {
                              return(0)
                            } else {
                              length(x)[[1]][1]
                            }
                          })
  
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
              , median_rule_cascade = median(rule_idx)
              , mean_rule_cascade = mean(rule_idx)
              , sd_rule_cascade = sd(rule_idx)
              , mean_rulelen = mean(rl_ln)
              , sd_rulelen = sd(rl_ln)
              , model = learner
              , model_accurate = model_accurate
              , model_type="inTrees"
              , begin_time = begin_time
              , completion_time = Sys.time()
  ))
  
}

sbrl_data_prep <- function(dat) {
  discrete <- character(0)
  continuous <- character(0)
  for(cl in names(dat)) {
    if(class(dat[[cl]]) != "factor") {
      continuous <- c(continuous, cl)
      dat[[cl]] <- binning(dat[[cl]], method="quantile", ordered = FALSE)
      levs <- levels(dat[[cl]])
      levs <- sub(",", ";", levs)
      levels(dat[[cl]]) <- levs
    } else {
      discrete <- c(discrete, cl)
    }
  }
  discrete <<- discrete
  continuous <<- continuous
  return(dat)
}

set_labels <- function(y_tr, pc) {
  factor(ifelse(y_tr != pc, 0, 1)) # error if not in this format  
}

time_per_explanation <- function(b_time, c_time, n_test) {
  difference <- c_time - b_time
  tpe <- as.numeric(c_time - b_time) / n_test
  if (attr(difference, "units") == "minutes") tpe <- tpe * 60
  if (attr(difference, "units") == "hours") tpe <- tpe * 3600
  return(tpe)
}

sbrl_generate_rule <- function(rule, reverse = FALSE) {
  if (rule == "{default}") return(rule)
  rules <- gsub("\\{|\\}", "", rule)
  rules <- strsplit(rules, ",")[[1]]
  var_names <- regmatches(rules, regexpr("[^=]*", rules, perl = TRUE))
  attrib_values <- mapply(sub, paste0(var_names, "="), "", rules)
  lwrs <- sapply(attrib_values, function(av) {
    rgx <- gregexpr("[\\(\\[][0-9.]*[0-9e\\+]*;", av)
    if (rgx < 0) return(NA)
    av <- regmatches(av, rgx)
    av <- gsub("\\[", "> ", gsub("\\(", ">= ", gsub(";", "", av)))
  })
  uprs <- sapply(attrib_values, function(av) {
    rgx <- gregexpr(";[0-9.]*[0-9e\\+]*", av)
    if (rgx < 0) return(NA)
    av <- regmatches(av, rgx)
    av <- gsub(";", "< ", av)
  })
  
  rule <- paste(
    apply(matrix(c(ifelse(is.na(lwrs) & is.na(uprs), paste0(var_names, " == '", attrib_values, "'"), NA)
                   , ifelse(is.na(lwrs), NA, paste(var_names, lwrs))
                   , ifelse(is.na(uprs), NA, paste(var_names, uprs))), ncol = 3)
          , 1, function(x) paste(x[!is.na(x)], collapse = " & "))
    , collapse = " & "
  )
  
  if (reverse) {
    return(paste("!(", rule, ")"))
  } else {
    return(rule)
  }
}

apply_rule <- function(rule, instances, algorithm) {
  if (rule %in% c("{default}", "X[,1]==X[,1]")) return(rep(TRUE, nrow(instances)))
  instances$idx <- rownames(instances)
  if (algorithm == "inTrees") {
    X <- instances
    covered <- tryCatch(filter(X, eval(parse_expr(rule))), error = function(e) select_all(X[0, ]) )
  } else {
    covered <- tryCatch(filter(instances, eval(parse_expr(rule))), error = function(e) select_all(instances[0, ]) )
  }
  return(ifelse(rownames(instances) %in% covered$idx, TRUE, FALSE))
}

sbrl_benchmark <- function(ds_container, classes, lambda, eta, rule_maxlen, nchain) {
  begin_time <- Sys.time()
  
  # transform for sbrl
  train_data <- sbrl_data_prep(ds_container$X_train)
  # using the quantile break points from the training set
  test_data <- as.data.frame(sapply(names(ds_container$X_test), function(nm) {
    x <- ds_container$X_test[[nm]]
    if(class(x) != "factor") {
      levs <- levels(train_data[[nm]])
      br <- attr(train_data[[nm]], "breaks")
      if (length(levs) == 4) {
        factor(ifelse(x <= br[2], levs[1]
                      , ifelse(x <= br[3], levs[2]
                               , ifelse(x <= br[4], levs[3], levs[4])))
               , levels = levs)
      } else {
        if (length(levs) == 3) {
          factor(ifelse(x <= br[2], levs[1]
                        , ifelse(x <= br[3], levs[2], levs[3]))
                 , levels = levs)
        } else {
          factor(ifelse(x <= br[2], levs[1], levs[2])
                 , levels = levs)
        }
      }
    } else {
      x
    }
  }))
  
  if (length(classes) > 2) {
    class_iter <- classes
  } else {
    class_iter <- positive_classes[i]
  }
  
  sbrls <- list()
  for (pc in class_iter) {
    # label for training
    train_label <- set_labels(ds_container$y_train, pc)
    # reset
    train_data$label <- NULL
    train_data <- cbind(train_data, label = train_label)
    test_label <- set_labels(ds_container$y_test, pc)
    
    sbrls[[pc]] <- list()
    sbrls[[pc]]$model <- sbrl(tdata=train_data, rule_minlen=1, rule_maxlen=rule_maxlen, 
                              minsupport_pos=0.05, minsupport_neg=0.05,
                              lambda=lambda, eta=eta, nchain=nchain)
    
    sbrls[[pc]]$preds <- predict(sbrls[[pc]]$model, tdata=test_data)$V2
    sbrls[[pc]]$rules_plus_default <- c(sbrls[[pc]]$model$rulenames, "{default}")
    sbrls[[pc]]$n_rules <- length(sbrls[[pc]]$model$rulenames) + 1
    sbrls[[pc]]$rule_idx <- sbrl_whichRule(sbrls[[pc]]$model, test_data)
    sbrls[[pc]]$rule_pos <- sapply(sbrls[[pc]]$rule_idx, function(x) {which(sbrls[[pc]]$model$rs$V1 == x)}) # create a positional index
    sbrls[[pc]]$rule_idx <- ifelse(sbrls[[pc]]$rule_idx == 0, sbrls[[pc]]$n_rules, sbrls[[pc]]$rule_idx) # makes rule extract easier
    sbrls[[pc]]$rule <- sbrls[[pc]]$rules_plus_default[sbrls[[pc]]$rule_idx]
    
    default_rule_pos <- max(sbrls[[pc]]$rule_pos)
    concatenate_rule <- function(rule_pos) {
      
      if (rule_pos == default_rule_pos) return("{default}")
      if (rule_pos == 1) return(sbrl_generate_rule(sbrls[[pc]]$rules_plus_default[sbrls[[pc]]$model$rs$V1[rule_pos]]))
      
      paste(
        paste(
          sapply(sapply(1:(rule_pos - 1), function(pos) {
            sbrls[[pc]]$rules_plus_default[sbrls[[pc]]$model$rs$V1[pos]]
          }), sbrl_generate_rule, reverse = TRUE)
          , collapse = " & "
          )
        , "&", sbrl_generate_rule(sbrls[[pc]]$rules_plus_default[sbrls[[pc]]$model$rs$V1[rule_pos]])
      )
    }
    sbrls[[pc]]$concatenated_rule <- sapply(sbrls[[pc]]$rule_pos, concatenate_rule)

    bangs <- gregexpr("!", sbrls[[pc]]$concatenated_rule)
    last_bang_pos <- sapply(bangs, function(b) b[length(b)])
    last_paren_pos <- unlist(sapply(1:length(last_bang_pos)
                                    , function(b) gregexpr(")"
                                                           , substr(sbrls[[pc]]$concatenated_rule[b]
                                                                    , start = last_bang_pos[b][[1]]
                                                                    , stop = nchar(sbrls[[pc]]$concatenated_rule[b])))))
    sbrls[[pc]]$rl_ln <- sapply(1:length(last_bang_pos), function(b) {
      bang <- 0
      if(bangs[[b]][1] != -1) bang <- length(bangs[[b]])
      last_pos <- ifelse(last_bang_pos[b] + last_paren_pos[b] <= -1, 0
                         , last_bang_pos[b] + last_paren_pos[b])
      and <- 0
      ands <- unlist(gregexpr("&", substr(sbrls[[pc]]$concatenated_rule[b]
                                          , start = last_pos
                                          , stop = nchar(sbrls[[pc]]$concatenated_rule[b]))))
      if (ands[1] != -1) and <- length(ands)
      rul_len <- 0
      if (bang + and > 1) rul_len <- bang + and
      return(rul_len)
    })
  }
  
  out <- list(this_i = i
              , this_r = r
              , this_run = length(random_states) * (i - 1) + r
              , random_state = random_states[r]
              , datasetname = datasetnames[i]
              , model_type = "sbrl"
              , begin_time = begin_time
              , completion_time = Sys.time())
  
  if (length(class_iter) == 1) { # binary classification
    out$label <- ifelse(sbrls[[pc]]$preds > 0.5, 1, 2)
    out$model_accurate <- ifelse(out$label == as.numeric(ds_container$y_test), 1, 0)  
    out$rule_idx = sbrls[[pc]]$rule_pos
    out$rl_ln <- sbrls[[pc]]$rl_ln
    out$unique_rules <- length(sbrls[[pc]]$model$rs$V1)
    out$n_rules_used <- length(unique(sbrls[[pc]]$rule_idx))
    out$base_rule <- sbrls[[pc]]$rule
    out$rule <- sbrls[[pc]]$concatenated_rule
    out$model <- sbrls[[pc]]$model
  } else { # multi-class
    out$label <- factor(classes[
      apply(sapply(sbrls, function(x) {x$preds}), 1, which.max)
      ]
      , levels = levels(ds_container$y_test))
    out$model_accurate <- ifelse(out$label == ds_container$y_test, 1, 0)
    
    get_multi_class_item <- function(item) {
      sapply(seq_along(out$label), function(x, item) {
        sapply(classes, function(pc, it) {
          sbrls[[pc]][[it]]}, it = item)[x, out$label[x]]
      }, item = item)  
    }
    
    out$label <- as.numeric(out$label)
    out$rule_idx <- get_multi_class_item("rule_pos")
    out$rl_ln <- as.vector(get_multi_class_item("rl_ln"))
    out$unique_rules <- sum(sapply(classes, function(pc) {
      length(sbrls[[pc]]$model$rs$V1)
    }))
    out$n_rules_used <-length(unique(apply(cbind(names(out$rule_idx), out$rule_idx), 1, paste0, collapse = "")))
    out$rule_idx <- as.vector(out$rule_idx)
    out$base_rule <- as.vector(get_multi_class_item("rule"))
    out$rule <- as.vector(get_multi_class_item("concatenated_rule"))
    out$model <- list()
    for (pc in classes) {
      out$model[[pc]] <- sbrls[[pc]]$model
    }
  }
  out$median_rule_cascade <- median(out$rule_idx)
  out$mean_rule_cascade <- mean(out$rule_idx)
  out$sd_rule_cascade <- sd(out$rule_idx)
  out$mean_rulelen <- mean(out$rl_ln)
  out$sd_rulelen <- sd(out$rl_ln)
  
  return(out)
}
