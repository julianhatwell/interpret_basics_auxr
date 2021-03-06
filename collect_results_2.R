# collect_results
library(dplyr)
source("compare_methods_utils.R")

# sensitivity analysis
count <- 1
for (i in seq_along(resfilesdirs)) {
  datasetname <- datasetnames[i]
  filepath <- normalizePath(file.path(resfilesdirs[i], "sensitivity\\"))
  filenames <- dir(filepath)
  for (fn in filenames) {
    
    # loading the results sheets
    if (!grepl("summary", fn)) {
      # get the meta data from the file name
      pattern <- gregexpr("rnst\\_[0-9]{3}", fn)
      pattern <- regmatches(fn, pattern)[[1]]
      rnst <- as.numeric(gsub("rnst\\_", "", pattern)) 
      pattern <- gregexpr("sp\\_[0-9]{1}\\.[0-9]{1,2}", fn)
      pattern <- regmatches(fn, pattern)[[1]]
      sp <- as.numeric(gsub("sp\\_", "", pattern))
      pattern <- gregexpr("ap\\_[0-9]{1}\\.[0-9]{1}", fn)
      pattern <- regmatches(fn, pattern)[[1]]
      ap <- as.numeric(gsub("ap\\_", "", pattern))
      pattern <- gregexpr("dpb\\_[0-9]{1,2}", fn)
      pattern <- regmatches(fn, pattern)[[1]]
      dpb <- as.numeric(gsub("dpb\\_", "", pattern))
      pattern <- gregexpr("sf\\_[0-9]", fn)
      pattern <- regmatches(fn, pattern)[[1]]
      sf <- as.numeric(gsub("sf\\_", "", pattern))
      pattern <- gregexpr("w\\_((chisq)|(nothing))", fn)
      pattern <- regmatches(fn, pattern)[[1]]
      w <- gsub("w\\_", "", pattern)
      
      # load a sheet
      results <- read.csv(normalizePath(file.path(filepath, fn))
                          , stringsAsFactors = FALSE)
      
      dataset <- rep(datasetname, nrow(results))
      support <- rep(sp, nrow(results))
      alpha_paths <- rep(ap, nrow(results))
      disc_path_bins <- rep(dpb, nrow(results))
      score_func <- rep(sf, nrow(results))
      weighting <- rep(w, nrow(results))
      
      results <- cbind(results
                       , dataset
                       , as.data.frame(cbind(support
                                             , alpha_paths
                                             , disc_path_bins
                                             , score_func
                                             , weighting)))
      
      # seem to be messing something up down the track
      dataset <- NULL
      support <- NULL
      alpha_paths <- NULL
      disc_path_bins <- NULL
      score_func <- NULL
      weighting <- NULL
      
      if (count == 1) {
        sens_results <- results
      } else {
        sens_results <- rbind(sens_results, results)
      }
      count <- count + 1
    }
  }
}

# comparative analysis
algorithms <- c("Anchors", "BRL", "defragTrees", "inTrees")
patt <- paste0("(", paste(algorithms, collapse = ")|("), ")")

first_comp <- TRUE
first_comp_summ <- TRUE
for (i in seq_along(resfilesdirs)) {
  datasetname <- datasetnames[i]
  filepath <- normalizePath(file.path(resfilesdirs[i]))
  filenames <- dir(filepath)
  for (fn in filenames) {
    if (grepl(patt, fn)) {
      # load a sheet
      results <- read.csv(normalizePath(file.path(filepath, fn))
                          , stringsAsFactors = FALSE)
      
      if (!grepl("summary", fn)) {
        if (first_comp == TRUE) {
          comp_results <- results
          first_comp <- FALSE
        } else {
          names(results) <- names(comp_results)
          comp_results <- rbind(comp_results, results)
        }
      } else { # summary
        if (!any(grepl("model_kappa", names(results)))) { # no kappa model column
          results <- results %>%
            mutate(model_kappa = NA)
        }
        if (!any(grepl("proxy_kappa", names(results)))) { # no kappa model column
          results <- results %>%
            mutate(proxy_kappa = NA)
        }
        if (!any(grepl("median_rule_cascade", names(results)))) { # no med rc column
          results <- results %>%
            mutate(median_rule_cascade = NA)
        }
        if (!any(grepl("^fidelity", names(results)))) { # no kappa model column
          results <- results %>%
            mutate(fidelity = NA)
        }
        if (!any(grepl("^proxy_performance", names(results)))) { # no med rc column
          results <- results %>%
            mutate(proxy_performance = NA)
        }
        # if (!any(grepl("sd_proxy_performance", names(results)))) { # no med rc column
        #   results <- results %>%
        #     mutate(median_rule_cascade = NA)
        # }
        results <- results %>%
          dplyr::select(X, dataset_name, algorithm, n_instances, n_rules, n_rules_used
                        , median_rule_cascade, mean_rule_cascade, sd_rule_cascade
                        , mean_rulelen, sd_rulelen, begin_time, completion_time
                        , forest_performance, sd_forest_performance
                        , model_kappa, proxy_performance, sd_proxy_performance
                        , proxy_kappa, fidelity, sd_fidelity)
        if (first_comp_summ == TRUE) {
          comp_summ_results <- results
          first_comp_summ <- FALSE
        } else {
          comp_summ_results <- rbind(comp_summ_results, results)
        }
      }
    }
  }
}
names(comp_results) <- sub("_name", "", names(comp_results))

comp_results$pretty.rule[is.na(comp_results$pretty.rule)] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == ""] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == "[]"] <- "{default}"

# adjust 0-length rules to be valued as 1, {default}
sens_results$rule.length <- ifelse(sens_results$pretty.rule == "{default}"
                                   , 0
                                   , sens_results$rule.length)

comp_results$rule.length <- ifelse(comp_results$pretty.rule == "{default}"
                                   , 0
                                   , comp_results$rule.length)

