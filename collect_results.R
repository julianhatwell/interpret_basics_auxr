source("dirs.R")

# main results collection
folders <- c("adult_small_samp"
          , "bankmark_samp"
          , "car"
          , "cardio"
          , "credit"
          , "german"
          , "lending_tiny_samp"
          , "nursery_samp"
          , "rcdv_samp")

meths <- c("Anchors", "BRL", "CHIRPS", "defragTrees", "inTrees")

np <- function(fld) {
  return(function(fl) {
    normalizePath(file.path(resfilesdir, fld, fl))  
  })
}

is_which_method <- function(x) {
  for (meth in meths) {
    res <- grepl(meth, x)
    if (res) return(meth)
  }
  return(NA)
}

# temp - correcting for mismatch in file names
for(fld in folders) {
  normpath <- np(fld)
  
  fls <- dir(normpath(""))
  for (fl in fls) {
    if (grepl("results_", fl)) {
      file.rename(normpath(fl), normpath(sub("results_", "", fl)))
    }
  }
}

# temp - correcting for mismatch in file names
for(fld in folders) {
  normpath <- np(fld)
  
  fls <- dir(normpath(""))
  for (fl in fls) {
    if (grepl("rndst_", fl)) {
      file.rename(normpath(fl), normpath(sub("rndst_", "rnst_", fl)))
    }
  }
}

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


# example results
precision_mat <- tapply(main_results$precision.tt.
                        , list(main_results$dataset_name
                               , main_results$algorithm
                               , main_results$random_state)
                        , mean)


mean_narm <- function(x) {
  mean(x, na.rm = TRUE)
}
apply(precision_mat, c(1, 2), mean_narm)
