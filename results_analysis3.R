library(dplyr)
library(ggplot2)
library(Rmpfr)
library(gridExtra)
library(tikzDevice)
library(stargazer)

# pexact function for post hoc wilcox tests
pexactfrsd <- function(d,k,n,option) {
  if (any(n < 1))                     stop("n out-of-bounds: min = 1")
  if (any(k < 2))                     stop("k out-of-bounds: min = 2")
  if (any(d < 0) || any(d > n*(k-1))) stop("d out-of-bounds: min,max = 0,n(k-1)")
  if (missing(option)) {option = "pvalue"}
  result <- 0
  for (h in 0:n) {
    sum1 <- chooseZ(n,h)/mpfr((pow.bigz(k,h) * pow.bigz(1-k,n)), precBits = 2048)
    sum2 <- 0
    for (s in 0:h) {
      if (any(k*s-d-h >= 0)) {
        if (option == "pvalue" || option == "cumulative no of compositions"){
          sum2 <- sum2 + (-1)^s * chooseZ(2*h,h+s) * chooseZ(k*s-d+h,k*s-d-h)}
        if (option == "probability" || option == "no of compositions"){
          sum2 <- sum2 + (-1)^s * chooseZ(2*h,h+s) * chooseZ(k*s-d+h-1,k*s-d-h)}
      }
    }
    result <- result + sum1 * sum2
  }
  if (any(d == 0) & option== "pvalue") return(1)
  if (any(d != 0) & option== "pvalue") return(as.numeric(2*result))
  if (any(d == 0) & option== "probability") return(as.numeric(result))
  if (any(d != 0) & option== "probability") return(as.numeric(2*result))
  if (any(d != 0) & 
      (option== "no of compositions" || option== "cumulative no of compositions") )            
    return(round(2*result*mpfr(pow.bigz(k*(k-1),n), precBits = 2048)))
  if (any(d == 0) & option== "no of compositions")            
    return(round(result*mpfr(pow.bigz(k*(k-1),n), precBits = 2048)))
  if (any(d == 0) & option== "cumulative no of compositions") 
    return(round(mpfr(pow.bigz(k*(k-1),n), precBits = 2048)))
}

# fudge until added absolute coverage, train and test set sizes to results
train_set_size <- integer(length(datasetnames))
test_set_size <- integer(length(datasetnames))

for (i in seq_along(datasetnames)) {
  train_set_size[i] <- get_train_test_sizes(i)[1]
  test_set_size[i] <- get_train_test_sizes(i)[2]
}

names(train_set_size) <- datasetnames
names(test_set_size) <- datasetnames

# sensitivity analysis
analysis_out <- list()

get_CHIRPS_analysis <- function(measure) {
  for (ds in datasetnames) {
    
    # results collection
    analysis_out[[ds]] <- list()
    
    # filter all results to one dataset
    analysis <- sens_results %>% 
      filter(dataset == ds) %>%
      mutate(support = factor(support)
             , alpha_paths = factor(alpha_paths)
             , disc_path_bins = factor(disc_path_bins)
             , score_func = factor(score_func)
             , weighting = weighting
             # covered is train set size * coverage
             , covered.tr. = (train_set_size[ds] * coverage.tr.)
             # covered correct is (covered + 1 (explanandum)) * stability
             # , cc.tr. = (covered.tr. + 1) * stability.tr.
             , cc.tr. = covered.tr. * precision.tr.
             # covered incorrect is covered - cc
             , ci.tr. = covered.tr. - cc.tr.
             , stability.tr.lapl. = (cc.tr. + 1) / (covered.tr. + n_classes[ds] + 1) # laplace and stability
             , odds.tr.lapl. = (cc.tr. + 1 + 1/n_classes[ds])/(ci.tr. + (n_classes[ds] - 1) + (n_classes[ds] - 1)/n_classes[ds])
             , lodds.tr.lapl. = log(odds.tr.lapl.)
             # covered is test set size - 1  * coverage (because of LOO)
             , covered.tt. = ((test_set_size[ds] - 1) * coverage.tt.)
             # covered correct is (covered + 1 (explanandum)) * stability
             # , cc.tt. = (covered.tt. + 1) * stability.tt.
             , cc.tt. = covered.tt. * precision.tt.
             # covered incorrect is covered - cc
             , ci.tt. = covered.tt. - cc.tt.
             , stability.tt.lapl. = (cc.tt. + 1) / (covered.tt. + n_classes[ds] + 1) # laplace and stability
             , odds.tt.lapl. = (cc.tt. + 1 + 1/n_classes[ds])/(ci.tt. + (n_classes[ds] - 1) + (n_classes[ds] - 1)/n_classes[ds])
             , lodds.tt.lapl. = log(odds.tt.lapl.))
    
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
    
    # which is the best mean measure
    analysis_out[[ds]][["top_mean_block"]] <- which.max(analysis_out[[ds]]$analysis_groups$mean)
    analysis_out[[ds]][["top_ranksum_block"]] <- which.max(analysis_out[[ds]]$analysis_groups$rank_sum)
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
    analysis_out[[ds]][["ds_bestmean_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out[[ds]][["top_mean_block"]]]
    analysis_out[[ds]][["ds_bestranksum_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out[[ds]][["top_ranksum_block"]]]
  }
  
  return(analysis_out)
}



# comparative analysis
get_comparative_analysis <- function(measure, analysis_in) {
  
  for (ds in datasetnames) {
    
    # isolate and transform for one dataset
    analysis <- comp_results %>% 
      filter(dataset_name == ds) %>%
      mutate(# covered is test set size - 1  * coverage (because of LOO)
        covered.tt. = ((test_set_size[ds] - 1) * coverage.tt.)
        # covered correct is (covered + 1 (explanandum)) * stability
        , cc.tt. = covered.tt. * precision.tt.
        # covered incorrect is covered - cc
        , ci.tt. = covered.tt. - cc.tt.
        , stability.tt.lapl. = (cc.tt. + 1) / (covered.tt. + n_classes[ds] + 1) # laplace and stability
        , odds.tt.lapl. = (cc.tt. + 1 + 1/n_classes[ds])/(ci.tt. + (n_classes[ds] - 1) + (n_classes[ds] - 1)/n_classes[ds])
        , lodds.tt.lapl. = log(odds.tt.lapl.))
    
    # transpose the raw values for analysis
    analysis_values <- tapply(analysis[, measure]
                              , analysis$algorithm
                              , identity)
    
    algos <- names(analysis_values) # no brl for cardio or nursery, no lending for brl or intrees

    analysis_values <- matrix(unlist(analysis_values), ncol = length(algos))
    analysis_values <- cbind(analysis_values, analysis_out[[ds]]$ds_bestranksum_values)
    
    # collect
    algos <- c(algos, "CHIRPS")
    colnames(analysis_values) <- algos
    analysis_in[[ds]][["comp_raw"]] <- analysis
    analysis_in[[ds]][["comp_values"]] <- analysis_values
    
    # friedman test
    analysis_in[[ds]][["comp_frd.tt"]] <- friedman.test(analysis_values)
    chisqstat <- analysis_in[[ds]][["comp_frd.tt"]]$statistic
    df1 <- analysis_in[[ds]][["comp_frd.tt"]]$parameter
    N <- nrow(analysis_values)
    fstat <- (chisqstat * (N - 1)) / (N * df1 - chisqstat)
    names(fstat) <- "Friedman F"
    df2 <- df1 * (N - 1)
    fpvalue <- pf(fstat, df1=df1, df2=df2, lower.tail = FALSE)
    
    # collect
    analysis_in[[ds]][["comp_frd.tt"]][["F.statistic"]] <- fstat
    analysis_in[[ds]][["comp_frd.tt"]][["F.p.value"]] <- fpvalue
    analysis_in[[ds]][["comp_frd.tt"]][["F.N"]] <- N
    analysis_in[[ds]][["comp_frd.tt"]][["F.df1"]] <- df1
    analysis_in[[ds]][["comp_frd.tt"]][["F.df2"]] <- df2
  }
  return(analysis_in)
}

results_in_detail <- function(analysis_in, rounding = 2) {
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
    mr <- apply(apply(-analysis_in[[ds]]$comp_values, 1, rank), 1, mean)
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
analysis_out <- get_CHIRPS_analysis(measure)

for (ds in datasetnames) {
  print(c(analysis_out[[ds]][["top_mean_block"]]
          , analysis_out[[ds]][["top_ranksum_block"]])
  )
}

measure <- "stability.tt.lapl."
stab_tt_analysis_out <- get_comparative_analysis(measure, analysis_out)
measure <- "precision.tt."
prec_tt_analysis_out <- get_comparative_analysis(measure, analysis_out)
measure <- "xcoverage.tt."
xcov_tt_analysis_out <- get_comparative_analysis(measure, analysis_out)
measure <- "rule.length"
rlen_tt_analysis_out <- get_comparative_analysis(measure, analysis_out)

results_in_detail(stab_tt_analysis_out, rounding = 4)


