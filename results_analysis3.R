library(dplyr)
library(ggplot2)

# fudge until added absolute coverage, train and test set sizes to results
train_set_size <- integer(length(datasetnames))
test_set_size <- integer(length(datasetnames))

for (i in seq_along(datasetnames)) {
  train_set_size[i] <- get_train_test_sizes(i)[1]
  test_set_size[i] <- get_train_test_sizes(i)[2]
}

names(train_set_size) <- datasetnames
names(test_set_size) <- datasetnames

laplace_corrections <- function(results_set, datasetname, sensitivity = TRUE) {
  # filter all results to one dataset
  analysis_out <- results_set %>% 
    filter(dataset == datasetname)
  
  if (sensitivity == TRUE) {
    analysis_out <- analysis_out %>% mutate(support = factor(support)
                            , alpha_paths = factor(alpha_paths)
                            , disc_path_bins = factor(disc_path_bins)
                            , score_func = factor(score_func)
                            , weighting = weighting)
  }
  
  analysis_out <- analysis_out %>% mutate(
           # covered is train set size * coverage
           covered.tr. = (train_set_size[datasetname] * coverage.tr.)
           # covered correct is (covered + 1 (explanandum)) * stability
           # , cc.tr. = (covered.tr. + 1) * stability.tr.
           , cc.tr. = covered.tr. * precision.tr.
           # covered incorrect is covered - cc
           , ci.tr. = covered.tr. - cc.tr.
           , stability.tr.lapl. = (cc.tr. + 1) / (covered.tr. + n_classes[datasetname] + 1) # laplace and stability
           , xcoverage.tr.lapl. = (covered.tr. + 1) / (train_set_size[datasetname] + n_classes[datasetname] + 1) # laplace and xcoverage
           , odds.tr.lapl. = (cc.tr. + 1 + 1/n_classes[datasetname])/(ci.tr. + (n_classes[datasetname] - 1) + (n_classes[datasetname] - 1)/n_classes[datasetname])
           , lodds.tr.lapl. = log(odds.tr.lapl.)
           # covered is test set size - 1  * coverage (because of LOO)
           , covered.tt. = ((test_set_size[datasetname] - 1) * coverage.tt.)
           # covered correct is (covered + 1 (explanandum)) * stability
           # , cc.tt. = (covered.tt. + 1) * stability.tt.
           , cc.tt. = covered.tt. * precision.tt.
           # covered incorrect is covered - cc
           , ci.tt. = covered.tt. - cc.tt.
           , stability.tt.lapl. = (cc.tt. + 1) / (covered.tt. + n_classes[datasetname] + 1) # laplace and stability
           , xcoverage.tt.lapl. = (covered.tt. + 1) / (test_set_size[datasetname] + n_classes[datasetname] + 1) # laplace and xcoverage
           , odds.tt.lapl. = (cc.tt. + 1 + 1/n_classes[datasetname])/(ci.tt. + (n_classes[datasetname] - 1) + (n_classes[datasetname] - 1)/n_classes[datasetname])
           , lodds.tt.lapl. = log(odds.tt.lapl.))
  return(analysis_out)
}

get_CHIRPS_analysis <- function(measure, results_set
                                , top_mean_block = NA
                                , top_ranksum_block = NA) {
  # sensitivity analysis
  analysis_out <- list()

  for (ds in datasetnames) {
  
  # results collection
  analysis_out[[ds]] <- list()
  
  # filter all results to one dataset
  analysis <- laplace_corrections(results_set, ds)
  
  # get the listing of all the sensitivity analyses by grid
  analysis_groups <- with(analysis, expand.grid(
    support = unique(support)
    , alpha = unique(alpha_paths)
    , bins = unique(disc_path_bins)
    , func = unique(score_func)
    , weights = unique(weighting)))
  
  # provide id numbers
  analysis_groups$id <- as.numeric(rownames(analysis_groups))
  
  # extract the values of "measure"
  analysis_values <- tapply(analysis[, measure]
                            , list(analysis$support
                                   , analysis$alpha_paths
                                   , analysis$disc_path_bins
                                   , analysis$score_func
                                   , analysis$weighting)
                            , identity)
  
  analysis_values <- matrix(unlist(analysis_values), ncol = nrow(analysis_groups))
  analysis_out[[ds]][["analysis_values"]] <- analysis_values
  
  # get a first taste of results
  analysis_groups$mean <- apply(analysis_values, 2, mean)
  analysis_groups$rank_mean <- colMeans(t(apply(analysis_values, 1, rank)))
  analysis_groups$rank_sum <- colSums(t(apply(analysis_values, 1, rank)))
  
  # plot the facet, hopefully showing lack of bins effect and use of chisq
  bins_plot <- ggplot(
    data = analysis_groups
    , aes(y = mean
          , x = id
          , shape = alpha
          , colour = func
          , size = support)) +
    geom_point() +
    scale_size_discrete(range = c(2, 4)) +
    theme_bw() +
    facet_grid(weights~bins)
  
  # collect plot
  analysis_out[[ds]][["bins_weights_facet"]] <- bins_plot
  
  # collect groups results
  analysis_out[[ds]][["analysis_groups"]] <- analysis_groups
  
  # which is the best mean measure discover or feed in from best stability
  if (is.na(top_mean_block)) {
    ds_top_mean_block <- which.max(analysis_out[[ds]]$analysis_groups$mean)
  } else {
    ds_top_mean_block <- top_mean_block
  }
  analysis_out[[ds]][["top_mean_block"]] <- ds_top_mean_block
  
  if (is.na(top_ranksum_block)) {
    ds_top_ranksum_block <- which.max(analysis_out[[ds]]$analysis_groups$rank_sum)
  } else {
    ds_top_ranksum_block <- top_ranksum_block
  }
  analysis_out[[ds]][["top_ranksum_block"]] <- ds_top_ranksum_block
  
  top_group_stats <- filter(analysis_groups, id == ds_top_ranksum_block) %>%
    dplyr::select(support, alpha, bins, func, weights)
  
  analysis_out[[ds]][["ds_bestmean_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out[[ds]][["top_mean_block"]]]
  analysis_out[[ds]][["ds_bestranksum_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out[[ds]][["top_ranksum_block"]]]
  
  sens_raw <- analysis %>% filter(support == top_group_stats$support
                                  , alpha_paths == top_group_stats$alpha
                                  , disc_path_bins == top_group_stats$bins
                                  , score_func == top_group_stats$func
                                  , weighting == top_group_stats$weights)
  
  analysis_out[[ds]][["sens_raw"]] <- sens_raw
  }
  
  mean_all_ds <- matrix(NA
                        , ncol = length(datasetnames)
                        , nrow = ncol(analysis_out$adult_small_samp$analysis_values))
  ranksum_all_ds <- matrix(NA
                           , ncol = length(datasetnames)
                           , nrow = ncol(analysis_out$adult_small_samp$analysis_values))
  
  for (j in seq_along(datasetnames)) {
    mean_all_ds[, j] <- analysis_out[[datasetnames[j]]]$analysis_groups$mean
    ranksum_all_ds[, j] <- analysis_out[[datasetnames[j]]]$analysis_groups$rank_sum
  }
  
  analysis_out$overall_bestmean_block <- which.max(apply(mean_all_ds, 1, mean)) # 5
  analysis_out$overall_bestranksum_block <- which.max(apply(ranksum_all_ds, 1, mean)) # 5
  
  for (ds in datasetnames) {
    analysis_out[[ds]][["overall_bestmean_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out$overall_bestmean_block]
    analysis_out[[ds]][["overall_bestranksum_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out$overall_bestranksum_block]
  }
  
  return(analysis_out)
}

# comparative analysis
get_comparative_analysis <- function(measure, results_set, CHIRPS_analysis) {
  
  for (ds in datasetnames) {
    
    # isolate and transform for one dataset
    analysis <- laplace_corrections(results_set, ds, sensitivity = FALSE)

    # transpose the raw values for analysis
    analysis_values <- tapply(analysis[, measure]
                              , analysis$algorithm
                              , identity)
    
    algos <- names(analysis_values) # no brl for cardio or nursery, no lending for brl or intrees

    analysis_values <- matrix(unlist(analysis_values), ncol = length(algos))
    analysis_values <- cbind(analysis_values, CHIRPS_analysis[[ds]]$ds_bestranksum_values)
    
    # collect
    algos <- c(algos, "CHIRPS")
    colnames(analysis_values) <- algos
    CHIRPS_analysis[[ds]][["comp_raw"]] <- analysis
    CHIRPS_analysis[[ds]][["comp_values"]] <- analysis_values
    
    # friedman test
    CHIRPS_analysis[[ds]][["comp_frd.tt"]] <- friedman.test(analysis_values)
    chisqstat <- CHIRPS_analysis[[ds]][["comp_frd.tt"]]$statistic
    df1 <- CHIRPS_analysis[[ds]][["comp_frd.tt"]]$parameter
    N <- nrow(analysis_values)
    fstat <- (chisqstat * (N - 1)) / (N * df1 - chisqstat)
    names(fstat) <- "Friedman F"
    df2 <- df1 * (N - 1)
    fpvalue <- pf(fstat, df1=df1, df2=df2, lower.tail = FALSE)
    
    # collect
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.statistic"]] <- fstat
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.p.value"]] <- fpvalue
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.N"]] <- N
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.df1"]] <- df1
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.df2"]] <- df2
  }
  return(CHIRPS_analysis)
}

results_in_detail <- function(analysis_in, rounding = 2, sgn = -1) {
  for (ds in datasetnames) {
    print(ds)
    print("mean qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, mean), rounding))
    print("sd qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, sd), rounding))
    print("min qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, min), rounding))
    print("max qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, max), rounding))
    print("med qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, median), rounding))
    mr <- apply(apply(sgn * analysis_in[[ds]]$comp_values, 1, rank), 1, mean)
    print("mean ranks")
    print(round(mr, rounding))
    print("Fried test")
    fstat <- analysis_in[[ds]][["comp_frd.tt"]][["F.statistic"]]
    df1 <- analysis_in[[ds]][["comp_frd.tt"]][["F.df1"]]
    df2 <- analysis_in[[ds]][["comp_frd.tt"]][["F.df2"]]
    N <- analysis_in[[ds]][["comp_frd.tt"]][["F.N"]]
    f.p.value <- analysis_in[[ds]][["comp_frd.tt"]][["F.p.value"]]
    print(c(fstat, df1, df2, N, f.p.value))
    print("post hoc tests")
    algs <- names(mr)[names(mr) != "CHIRPS"]
    k <- analysis_in[[ds]]$comp_frd.tt$parameter + 1
    md <- mr[algs] - mr["CHIRPS"]
    print("rank diff")
    print(round(md, rounding))
    z <- (md) / sqrt((k * (k + 1)) / (6 * N))
    ztest <- pnorm(z, lower.tail = FALSE)
    print("z stat")
    print(z)
    print("post-hoc z-test")
    print(ztest)
    print("reject null bonferroni")
    print(ztest < 0.025/df1)
    print(0.025/df1)
  }
}

measures <- c("stability.tr.lapl.", "precision.tt.", "stability.tt.lapl.", "xcoverage.tt.", "rule.length")
measure <- "stability.tr.lapl."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results)

for (ds in datasetnames) {
  print(c(CHIRPS_analysis[[ds]][["top_mean_block"]]
          , CHIRPS_analysis[[ds]][["top_ranksum_block"]])
  )
}

measure <- "stability.tt.lapl."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results)
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)

measure <- "precision.tt."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = CHIRPS_analysis[[ds]][["top_mean_block"]]
                                       , top_ranksum_block = CHIRPS_analysis[[ds]][["top_ranksum_block"]])
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)

measure <- "xcoverage.tt.lapl."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = CHIRPS_analysis[[ds]][["top_mean_block"]]
                                       , top_ranksum_block = CHIRPS_analysis[[ds]][["top_ranksum_block"]])
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)

measure <- "rule.length"
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = CHIRPS_analysis[[ds]][["top_mean_block"]]
                                       , top_ranksum_block = CHIRPS_analysis[[ds]][["top_ranksum_block"]])
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4, sgn = 1) # positive ranking. shortest wins

measure <- "recall.tt."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = CHIRPS_analysis[[ds]][["top_mean_block"]]
                                       , top_ranksum_block = CHIRPS_analysis[[ds]][["top_ranksum_block"]])
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)

measure <- "f1.tt."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = CHIRPS_analysis[[ds]][["top_mean_block"]]
                                       , top_ranksum_block = CHIRPS_analysis[[ds]][["top_ranksum_block"]])
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)


measure <- "rule.length"
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = CHIRPS_analysis[[ds]][["top_mean_block"]]
                                       , top_ranksum_block = CHIRPS_analysis[[ds]][["top_ranksum_block"]])
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4, sgn = 1) # positive ranking. shortest wins

for (ds in datasetnames) {
  A_rules <- with(tt_analysis[[ds]]$comp_raw
                  , tt_analysis[[ds]]$comp_raw[algorithm == "Anchors"
                                              , c("pretty.rule", "rule.length")])

  B_rules <- with(tt_analysis[[ds]]$comp_raw
                  , tt_analysis[[ds]]$comp_raw[algorithm == "BRL"
                                              , c("pretty.rule", "rule.length")])

  D_rules <- with(tt_analysis[[ds]]$comp_raw
                  , tt_analysis[[ds]]$comp_raw[algorithm == "defragTrees"
                                              , c("pretty.rule", "rule.length")])
  
  I_rules <- with(tt_analysis[[ds]]$comp_raw
                  , tt_analysis[[ds]]$comp_raw[algorithm == "inTrees"
                                              , c("pretty.rule", "rule.length")])
  
  C_rules <- with(tt_analysis[[ds]]$sens_raw
                  , tt_analysis[[ds]]$sens_raw[, c("pretty.rule", "rule.length")])
  
  uniques <- c(length(unique(A_rules$pretty.rule))
               , length(unique(B_rules$pretty.rule))
               , length(unique(D_rules$pretty.rule))
               , length(unique(I_rules$pretty.rule))
               , length(unique(C_rules$pretty.rule)))
  
  names(uniques) <- c("A", "B", "D", "I", "C")
  print(ds)
  print(uniques)
}

measure <- "mean_rule_cascade"
analysis <- comp_summ_results
# transpose the raw values for analysis
analysis_values <- tapply(analysis[, measure]
                          , list(analysis$dataset_name
                                 , analysis$algorithm)
                          , identity)

print(analysis_values)

analysis <- comp_summ_results
# transpose the raw values for analysis
analysis_values <- tapply(analysis[, measure]
                          , list(analysis$dataset_name
                                 , analysis$algorithm)
                          , identity)

print(analysis_values)
