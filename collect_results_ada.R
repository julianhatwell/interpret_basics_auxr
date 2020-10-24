rm(list=ls())
# data settings
# data_dir <- "C:\\Users\\id126493\\Documents\\GitHub\\explain_te\\CHIRPS\\datafiles\\"
# data_dir <- "C:\\Users\\Crutt\\Documents\\GitHub\\explain_te\\CHIRPS\\datafiles\\"
data_dir <- "~/Documents/github/explain_te/CHIRPS/datafiles/"
# datafilesdir <- "c:\\Dev\\Study\\Python\\explain_te\\forest_surveyor\\datafiles\\"

# project_dir <- "V:\\whiteboxing\\tests\\"
# project_dir <- "V:\\whiteboxing\\"
# project_dir <- "C:\\Users\\Crutt\\Documents\\whiteboxing\\tests\\"
project_dir <- "/datadisk/whiteboxing/2020Ada2/"  # folder Ada1 or Ada2?

options(max.print=20*72)
algo_variant <- "Ada-WHIPS" # there is a hardcoded version that must be changed. Find and replace.
algo_input_csv <- c("Anchors", "BRL", "CHIRPS", "defragTrees", "lore")
patt_input_csv <- paste0("(", paste(algo_input_csv, collapse = ")|("), ")")
algorithms <- c("Anchors", "BRL", algo_variant, "defragTrees", "LORE")
algorithms <- sort(algorithms)

source("data_files_mgmt.R")

# comparative analysis
library(dplyr)
library(tidyr)
library(PMCMRplus)
library(cowplot)
library(rlang)

standard_rounding <- 2

resfilesdirs <- paste0(project_dir, datasetnames, pathsep)

first_comp <- TRUE
first_comp_summ <- TRUE

for (i in seq_along(resfilesdirs)) {
  datasetname <- datasetnames[i]
  filepath <- normalizePath(file.path(resfilesdirs[i]))
  filenames <- dir(filepath)
  for (fn in filenames) {
    if (fn == "y_test.csv") {
      test_set_size[datasetname] <- nrow(read.csv(normalizePath(file.path(filepath, fn))
                                                  , stringsAsFactors = FALSE))
    }
    if (fn == "y_train.csv") {
      train_set_size[datasetname] <- nrow(read.csv(normalizePath(file.path(filepath, fn))
                                                   , stringsAsFactors = FALSE))
    }
    if (grepl(patt_input_csv, fn)) {
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
comp_results$algorithm <- ifelse(comp_results$algorithm == "greedy_stab", algo_variant, comp_results$algorithm)
comp_results$algorithm <- ifelse(comp_results$algorithm == "lore", "LORE", comp_results$algorithm)
# comp_results$algorithm <- ifelse(comp_results$algorithm == "GTB_HIPS", "gbt-HIPS", comp_results$algorithm)
comp_results$pretty.rule[is.na(comp_results$pretty.rule)] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == ""] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == "[]"] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == "X[,1]==X[,1]"] <- "{default}"
comp_results$dataset[comp_results$dataset == "lending_tiny_samp"] <- "lending"
comp_results$dataset[comp_results$dataset == "bankmark"] <- "bank"

# sort out the lengths of defragTrees rules
relength_pretty_rule <- function(pretty_rule) {
  length(unique(sapply(strsplit(pretty_rule, " AND ")[[1]]
                       , function(x) {
                         string <- substr(x, 1, regexpr(" ", x) - 1)
                         underscore <- regexpr("_", x)
                         if (underscore < 0) underscore <- nchar(string) + 1
                         string <- substr(x, 1, underscore - 1)
                         string
                       })))
}

comp_results[comp_results$algorithm == "defragTrees", "rule.length"] <- 
  sapply(comp_results[comp_results$algorithm == "defragTrees", "pretty.rule"]
         , relength_pretty_rule)

# ensure {default} is valued at zero
comp_results$rule.length <- ifelse(comp_results$pretty.rule %in% c("{default}", "")
                                   , 0
                                   , comp_results$rule.length)

# calculate weighted excl.cov
comp_results <- within(comp_results, {
  wxcoverage.tr. <- xcoverage.tr. * (nci.tr./(ci.tr. + nci.tr.))
  wxcoverage.tt. <- xcoverage.tt. * (nci.tt./(ci.tt. + nci.tt.))
})

get_func_of <- function(func, metric, ...) {
  tapply(comp_results[[metric]]
         , list(comp_results$dataset
                , comp_results$algorithm)
         , func, ...)
}

get_mean_of <- function(metric) {
  get_func_of(mean, metric, na.rm = TRUE)
}
get_sd_of <- function(metric) {
  get_func_of(sd, metric, na.rm = TRUE)
}
get_median_of <- function(metric) {
  get_func_of(median, metric, na.rm = TRUE)
}

lwrq <- function(x) {
  quantile(x, probs = 0.25, na.rm = TRUE)
}
uprq <- function(x) {
  quantile(x, probs = 0.75, na.rm = TRUE)
}

get_lwrq_of <- function(metric) {
  get_func_of(lwrq, metric)
}
get_uprq_of <- function(metric) {
  get_func_of(uprq, metric)
}

above.75 <- function(x) {
  mean(x > 0.75) # shouldn't be any NA to rm.
}

get_above.75_of <- function(metric) {
  get_func_of(above.75, metric)
}

equal.zero <- function(x) {
  mean(x > 0) # shouldn't be any NA to rm.
}

get_equal.zero_of <- function(metric) {
  get_func_of(equal.zero, metric)
}

test_set_size_sqrt <- sapply(test_set_size, function(x) {
  sqrt(min(1000, x))
})

get_mean_ranks_of <- function(meas) {
  qmeas <- quo(meas)
  cres <- dplyr::select(comp_results, instance_id, dataset, !! qmeas, algorithm) %>%
    group_by(dataset) %>%
    spread(algorithm, !! qmeas)
  cres[, algorithms] <- t(apply(-cres[, algorithms], 1, rank))
  summarise(cres
            , Anchors = mean(Anchors)
            , BRL = mean(BRL)
            , `Ada-WHIPS` = mean(`Ada-WHIPS`)
            , defragTrees = mean(defragTrees)
            , LORE = mean(LORE))
}

post_hoc_ztest <- function(meas) {
  mr <- select(get_mean_ranks_of(meas), -dataset)
  
  top_rank <- as.vector(as.matrix(mr)[matrix(c(1:nrow(mr), apply(mr, 1, which.min)), ncol = 2)])
  scnd_rank <- c()
  for (i in 1:nrow(mr)) {
    scnd_rank[i] <- mr[i, ][, as.vector(mr[i, ] != top_rank[i])][which.min(mr[i, ][, as.vector(mr[i, ] != top_rank[i])])]
  }
  scnd_rank <- unlist(scnd_rank)
  ranks <- as.data.frame(t(matrix(c(top_rank, scnd_rank), ncol = 2)))
  names(ranks) <- datasetnames
  print(ranks)
  md <- scnd_rank - top_rank
  k <- c(6, 6, 6, 5, 6, 6, 4, 5, 6)
  z <- (md) / sqrt((k * (k + 1)) / (6 * test_set_size))
  print(z)
  ztest <- pnorm(z, lower.tail = FALSE)
  print(ztest)
  print("reject null bonferroni")
  print(ztest < 0.025/k)
}

get_fidelity <- function() {
  tapply(comp_results$predicted.class == comp_results$target.class
         , list(comp_results$dataset
                , comp_results$algorithm)
         , mean)
}

round(get_fidelity(), standard_rounding)

meas <- "stability.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), standard_rounding)
round(st_err, standard_rounding)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)
round(get_above.75_of(meas), standard_rounding)

meas <- "precision.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), standard_rounding)
round(st_err, standard_rounding)
get_mean_ranks_of(meas)

meas <- "coverage.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), standard_rounding)
round(st_err, standard_rounding)
get_mean_ranks_of(meas)

meas <- "wxcoverage.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), standard_rounding)
round(st_err, standard_rounding)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)

meas <- "rule.length"
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), standard_rounding)
round(st_err, standard_rounding)
round(get_equal.zero_of("rule.length"), standard_rounding)

meas <- "elapsed_time"
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), standard_rounding)
round(st_err, standard_rounding)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)
