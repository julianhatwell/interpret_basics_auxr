library(dplyr)
library(ggplot2)
library(Rmpfr)

# sensitivity analysis
analysis_out <- list()

measures <- c("precision.tt.", "stability.tt.", "xcoverage.tt.", "rule.length")
measure <- measures[2]

for (ds in datasetnames) {
  analysis <- sens_results %>% 
    filter(dataset == ds) %>%
    mutate(support = factor(support)
           , alpha_paths = factor(alpha_paths)
           , disc_path_bins = factor(disc_path_bins)
           , score_func = factor(score_func)
           , weighting = weighting)
  
  analysis_groups <- with(analysis, expand.grid(
    support = unique(support)
    , alpha = unique(alpha_paths)
    , bins = unique(disc_path_bins)
    , func = unique(score_func)
    , weights = unique(weighting)))
  
  analysis_groups$id <- as.numeric(rownames(analysis_groups))
  
  analysis_values <- tapply(analysis[, measure]
                            , list(analysis$support
                                   , analysis$alpha_paths
                                   , analysis$disc_path_bins
                                   , analysis$score_func
                                   , analysis$weighting)
                            , identity)
  
  # n_test <- length(analysis_values[[1]])
  analysis_values <- matrix(unlist(analysis_values), ncol = nrow(analysis_groups))
  
  analysis_groups$mean <- apply(analysis_values, 2, mean)
  analysis_groups$rank_mean <- colMeans(t(apply(-analysis_values, 1, rank)))
  analysis_groups$rank_sum <- colSums(t(apply(-analysis_values, 1, rank)))
  
  bins_plot <- ggplot(
    data = analysis_groups
    , aes(y = rank_sum # I(1-mean)
          , x = id
          , shape = alpha
          , colour = func
          , size = support)) +
    geom_point() + 
    # scale_color_grey() +
    scale_size_discrete(range = c(2, 4)) +
    theme_bw() +
    facet_grid(weights~bins)

  analysis_groups_4 <- analysis_groups %>% filter(bins == 4)
  
  analysis_groups_8 <- analysis_groups %>% filter(bins == 8)
  
  
  analysis_bins_check <- inner_join(analysis_groups_4, analysis_groups_8
                                    , by = c("support", "alpha", "func", "weights")
                                    , suffix = c("_4", "_8")) %>%
    select(-bins_4, -bins_8)
  
  # bins_wilcox <- wilcox.test(analysis_bins_check$mean_4
  #                           , analysis_bins_check$mean_8
  #                           , paired = TRUE)
  bins_wilcox <- wilcox.test(analysis_values[, analysis_bins_check$id_4]
                             , analysis_values[, analysis_bins_check$id_8]
                             , paired = TRUE
                             , conf.int = TRUE
                             , conf.level = 0.99)
  
  bins_effect_size <- mean(analysis_values[, analysis_bins_check$id_4] -
         analysis_values[, analysis_bins_check$id_8]) /
    sd(analysis_values[, analysis_bins_check$id_4] -
         analysis_values[, analysis_bins_check$id_8])
  
  # we know we only need to keep bins_4
  analysis_groups_chisq <- analysis_groups_4 %>% filter(weights == "chisq")
  analysis_groups_nothing <- analysis_groups_4 %>% filter(weights == "nothing")
  
  analysis_weights_check <- inner_join(analysis_groups_chisq, analysis_groups_nothing
                                    , by = c("support", "alpha", "func")
                                    , suffix = c("_chisq", "_nothing")) %>%
    select(-weights_chisq, -weights_nothing)
  
  # weights_wilcox <- wilcox.test(analysis_weights_check$mean_chisq
  #                              , analysis_weights_check$mean_nothing
  #                              , paired = TRUE)
  
  weights_wilcox <- wilcox.test(analysis_values[, analysis_weights_check$id_chisq]
              , analysis_values[, analysis_weights_check$id_nothing]
              , paired = TRUE
              , conf.int = TRUE
              , conf.level = 0.99)
  
  weights_effect_size <- mean(analysis_values[, analysis_weights_check$id_chisq] -
         analysis_values[, analysis_weights_check$id_nothing]) /
    sd(analysis_values[, analysis_weights_check$id_chisq] -
         analysis_values[, analysis_weights_check$id_nothing])
  
  weights_plot <- ggplot(
    data = analysis_groups_4
    , aes(y = I(1-mean)
          , x = id
          , shape = alpha
          , colour = func
          , size = support)) +
    geom_point() + 
    scale_size_discrete(range = c(2, 4)) +
    theme_bw() +
    facet_grid(.~weights)
  
  analysis_included_values <- analysis_values[, analysis_groups_chisq$id]
  
  analysis_groups_chisq$mean <- apply(analysis_included_values, 2, mean)
  analysis_groups_chisq$rank_mean <- colMeans(t(apply(-analysis_included_values, 1, rank)))
  
  best_blocks_plot <- ggplot(
    data = analysis_groups_chisq
    , aes(y = mean_rank
          , x = id
          , shape = alpha
          , colour = func
          , size = support)) +
    geom_point() + 
    scale_size_discrete(range = c(2, 4)) +
    theme_bw()
  
  analysis_out[[ds]] <- list()
  analysis_out[[ds]][["all_values"]] <- analysis_values
  analysis_out[[ds]][["bins_plot"]] <- bins_plot
  analysis_out[[ds]][["bins_check"]] <- analysis_bins_check
  analysis_out[[ds]][["bins_wilcox"]] <- bins_wilcox
  analysis_out[[ds]][["bins_effect_size"]] <- bins_effect_size
  analysis_out[[ds]][["weights_plot"]] <- weights_plot
  analysis_out[[ds]][["weights_check"]] <- analysis_weights_check
  analysis_out[[ds]][["weights_wilcox"]] <- weights_wilcox
  analysis_out[[ds]][["weights_effect_size"]] <- weights_effect_size
  analysis_out[[ds]][["best_blocks_plot"]] <- best_blocks_plot
  analysis_out[[ds]][["sens_blocks"]] <- analysis_groups_chisq
  analysis_out[[ds]][["sens_values"]] <- analysis_included_values
  analysis_out[[ds]][["sens_raw"]] <- sens_results %>% filter(
    dataset_name == ds & disc_path_bins == 4 & weighting == "chisq"
  )
  analysis_out[[ds]][["sens_frd.tt"]] <- friedman.test(analysis_included_values)
}

# comparative analysis
for (ds in datasetnames) {
  
  analysis <- comp_results %>% 
    filter(dataset_name == ds)

  analysis_values <- tapply(analysis[, measure]
                      , analysis$algorithm
                      , identity)
  analysis_values <- matrix(unlist(analysis_values), ncol = length(algorithms))
  
  analysis_out[[ds]][["comp_frd.tt"]] <- friedman.test(analysis_values)
  # CHIRPS post hoc tests here
  
  analysis_out[[ds]][["comp_raw"]] <- analysis
  analysis_out[[ds]][["comp_values"]] <- analysis_values
  
}

for (ds in datasetnames) {
  print(ds)
  print(analysis_out[[ds]]$bins_wilcox)
  print(analysis_out[[ds]]$bins_effect_size)
  print("\n")
  }

for (ds in datasetnames) {
  print(ds)
  print(analysis_out[[ds]]$weights_wilcox)
  print(analysis_out[[ds]]$weights_effect_size)
  print("\n")}

# exact p-values for sensitivity
for (ds in datasetnames) {
  rs <- with(analysis_out[[ds]]
             , sens_blocks[order(-sens_blocks$rank_sum),])[, "rank_sum"]
  k <- length(rs)
  d <- (rs[1] - rs)[-1]
  n <- nrow(analysis_out[[ds]]$sens_values)
  analysis_out[[ds]]$sens_posthoc_pvalues # p value for best is signif better than others
}

# p-value to beat for sens post hoc
0.05 / 17

# p-value to beat for CHIRPS post - hoc
0.05 / 4

with(analysis_out$adult_small_samp
     , sens_blocks[order(-sens_blocks$mean),])
with(analysis_out$adult_small_samp
     , sens_blocks[order(-sens_blocks$rank_sum),])

analysis_out$adult_small_samp$sens_frd.tt

rs <- with(analysis_out$adult_small_samp
     , sens_blocks[order(-sens_blocks$rank_sum),])[, "rank_sum"]
k <- length(rs)
d <- (rs[1] - rs)[-1]
n <- nrow(analysis_out$adult_small_samp$sens_values)

# pexactfrsd(d,k,n, "pvalue")



# comp_results_means <- apply(analysis_values, 2, mean)
# com_results_mean_ranks <- colMeans(t(apply(-analysis_values, 1, rank)))



  