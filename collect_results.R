source("dirs.R")
source("collect_results_utils.R")

first <- TRUE
for (i in seq_along(folders)) {
  fld <- folders[i]
  normpath <- np(fld)
  fls <- dir(normpath(""))
  for (fl in fls) {
    meth <- is_which_method(fl)
    if (!is.na(meth)) {
      print(fl)
      rnst <- substr(fl
                     , start = regexpr("rnst_", fl)[1] + 5
                     , stop = nchar(fl) - 4)
      results <- read.csv(normpath(fl))
      nr <- nrow(results)
      random_state <- rep(rnst, nr)
      
      # temp fix code because forget cov and xcov
      if (meth %in% c("BRL", "inTrees")) {
        results <- cbind(results[, c("X"                 
                               , "dataset_name"
                               , "instance_id"
                               , "algorithm"
                               , "pretty_rule"
                               , "rule_length"
                               , "pred_class"
                               , "pred_class_label"
                               , "target_class"
                               , "target_class_label"
                               , "forest_vote_share"
                               , "prior"
                               , "precision_tr"
                               , "stability_tr"
                               , "recall_tr"
                               , "f1_tr"
                               , "accuracy_tr"
                               , "lift_tr")]
                         , coverage.tr. = rep(NA, nr)
                         , xcoverage.tr. = rep(NA, nr)
                         , results[, c("kl_div_tr"
                                       , "precision_tt"
                                       , "stability_tt"
                                       , "recall_tt"
                                       , "f1_tt"
                                       , "accuracy_tt"
                                       , "lift_tt")]
                         , coverage.tt. = rep(NA, nr)
                         , xcoverage.tt. = rep(NA, nr)
                         , results[, c("kl_div_tt")])
      }
      
      # temp fix because I stupidly named the dataset containers in Python
      results$dataset_name <- sub("_data", "", results$dataset_name)
      
      # don't forget this bit - not to be deleted when temp fix is removed
      results <- cbind(results, random_state)
      
      if (first) {
        main_results <- results
      } else {
        names(results) <- names(main_results)
        main_results <- rbind(main_results, results)
      }
      first <- FALSE
    }
  }
}

dfs <- list()
dps <- list()
for (i in seq_along(metrics)) {
  md <- mean_df(results_mat(main_results[, metrics_names[i]]))
  dfs[[metrics[i]]] <- md
  dps[[metrics[i]]] <- my_dotplot(dat = md, metric = metrics[i])
}
dps$stability
dps$xcoverage
dps$recall
dps$f1
