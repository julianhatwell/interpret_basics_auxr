library(dplyr)
library(ggplot2)
library(Rmpfr)
library(gridExtra)
library(tikzDevice)

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
    
    analysis_out <- analysis_in
    
    analysis <- comp_results %>% 
      filter(dataset_name == ds) %>%
      mutate(# covered is test set size - 1  * coverage (because of LOO)
        covered.tt. = ((test_set_size[ds] - 1) * coverage.tt.)
        # covered correct is (covered + 1 (explanandum)) * stability
        # , cc.tt. = (covered.tt. + 1) * stability.tt.
        , cc.tt. = covered.tt. * precision.tt.
        # covered incorrect is covered - cc
        , ci.tt. = covered.tt. - cc.tt.
        , stability.tt.lapl. = (cc.tt. + 1) / (covered.tt. + n_classes[ds] + 1) # laplace and stability
        , odds.tt.lapl. = (cc.tt. + 1 + 1/n_classes[ds])/(ci.tt. + (n_classes[ds] - 1) + (n_classes[ds] - 1)/n_classes[ds])
        , lodds.tt.lapl. = log(odds.tt.lapl.))
    
    
    analysis_values <- tapply(analysis[, measure]
                              , analysis$algorithm
                              , identity)
    
    algos <- names(analysis_values) # no brl for cardio or nursery, no lending for brl or intrees
    bonf.conf.level <- 1 - 0.05/length(algos)
    
    
    analysis_values <- matrix(unlist(analysis_values), ncol = length(algos))
    analysis_values <- cbind(analysis_values, analysis_out[[ds]]$ds_bestranksum_values)
    
    algos <- c(algos, "CHIRPS")
    colnames(analysis_values) <- algos
    
    analysis_out[[ds]][["comp_frd.tt"]] <- friedman.test(analysis_values)
    # CHIRPS post hoc tests here
    wilcos <- list()
    for (j in 1:(length(algos)-1)) {
      wilcos[[algos[j]]] <- wilcox.test(analysis_values[,length(algos)], analysis_values[, j]
                                        , paired = TRUE, conf.int = TRUE, conf.level = bonf.conf.level)
      
      wilcos[[algos[j]]]$data.name <- paste(algos[j], "and CHIRPS")
      
      rank_sum_total <- tail(cumsum(1:nrow(analysis_values)), 1)
      
      wilcos[[algos[j]]]$effect_size <- wilcos[[algos[j]]]$statistic/rank_sum_total
      
      # exact p values
      k <- length(algos)
      n <- nrow(analysis_values)
      rank_sums <- apply(apply(-analysis_values, 1, rank), 1, sum)
      d <- abs(rank_sums[j] - rank_sums["CHIRPS"])
      wilcos[[algos[j]]]$p.value <- pexactfrsd(d,k,n, "pvalue")
    }
    
    analysis_out[[ds]][["comp_wilcox.tt"]] <- wilcos
    analysis_out[[ds]][["comp_raw"]] <- analysis
    analysis_out[[ds]][["comp_values"]] <- analysis_values
  }
  return(analysis_out)
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
analysis_out <- get_comparative_analysis(measure, analysis_out)
