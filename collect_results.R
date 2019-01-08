source("dirs.R")
source("collect_results_utils.R")
library(tikzDevice)
library(stargazer)

random_states <- 123:(123 + 29)
oob <- numeric(30)
accu <- numeric(30)
coka <- numeric(30)
ntree <- integer(30)

first <- TRUE
for (i in seq_along(folders)) {
  fld <- folders[i]
  normpath <- np(fld)
  fls <- dir(normpath(""))
  
  for (fl in fls) {
    if (is_forestperf_file(fl)) {
      rnst <- get_rnst_from_filename(fl)
      if (rnst < 153) {# silly me, did one too many
        js <- fromJSON(file=normpath(fl))
        oob[which(random_states == rnst)] <- js$oobe_score
        accu[which(random_states == rnst)] <- js$main$test_accuracy
        coka[which(random_states == rnst)] <- js$main$test_kappa
      }
    }
    if (is_bestparams_file(fl)) {
      rnst <- get_rnst_from_filename(fl)
      if (rnst < 153) {# silly me, did one too many
        jp <- fromJSON(file=normpath(fl))
        ntree[which(random_states == rnst)] <- jp$n_estimator
      }
    }
  }
  rftun <- data.frame(dataset_name = rep(fld, 30)
                      , random_state = random_states
                      , ntree = ntree
                      , oob = oob
                      , accuracy = accu
                      , kappa = coka)
  if (first) {
    
    rfperf_results <- rftun
  } else {
    rfperf_results <- rbind(rfperf_results, rftun)
  }
  first <- FALSE
}

rfperf_mat(rfperf_results, "ntree", mean)
rfperf_mat(rfperf_results, "ntree", Mode)
rfperf_mat(rfperf_results, "ntree", median)
rfperf_mat(rfperf_results, "ntree", sd)

oob_mn <- rfperf_mat(rfperf_results, "oob", mean)
oob_sd <- rfperf_mat(rfperf_results, "oob", sd)


acc_mn <- rfperf_mat(rfperf_results, "accuracy", mean)
acc_sd <- rfperf_mat(rfperf_results, "accuracy", sd)

coka_mn <- rfperf_mat(rfperf_results, "kappa", mean)
coka_sd <- rfperf_mat(rfperf_results, "kappa", sd)

rfperf_summary <- data.frame(OOB = oob_mn
                             , oob_sd
                             , Accuracy = acc_mn
                             , acc_sd
                             , Cohensk = coka_mn
                             , coka_sd)

stargazer(rfperf_summary
          , summary = FALSE
          , rownames = FALSE)

for (fld in folders) {
  tikz(file=paste0("C:\\Users\\id126493\\OneDrive\\Documents\\PhD\\", fld, "_rfntree.tikz")
       , width = 0.75, height=0.25)
  print(my_histo(subset(rfperf_results, dataset_name==fld)))
  dev.off()
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
      rnst <- get_rnst_from_filename(fl)
      if (rnst < 153) {# silly me, did one too many
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
main_results$algorithm <- ifelse(main_results$algorithm == "defragTrees"
                                 , "dfrgTrs"
                                 , main_results$algorithm)

main_results$algorithm <- ifelse(main_results$algorithm == "greedy_stab"
                                 , "CHIRPS"
                                 , main_results$algorithm)


for (i in seq_along(metrics)) {
  md <- mean_df(results_mat(main_results, metrics_names[i]))
  md$qmeasure <- metrics[i]
  if (i == 1) {
    results_long <- md
  } else {
    results_long <- rbind(results_long, md)
  }
}

my_dotplot(subset(results_long
                  , qmeasure %in% c("precision"
                                    , "recall"
                                    , "stability"
                                    , "xcoverage")))

dps$stability
dps$xcoverage
dps$recall
dps$f1
