library(dplyr)

plots_b <- list()
bins_check <- list()
plots_w <- list()
weights_check <- list()
plots_r <- list()

for (ds in datasetnames) {
  analysis <- main_results %>% 
    filter(dataset == ds) %>%
    mutate(support = factor(support)
           , alpha_paths = factor(alpha_paths)
           , disc_path_bins = factor(disc_path_bins)
           , score_func = factor(score_func)
           , weighting = weighting)
  measures <- c("stability.tt.", "xcoverage.tt.", "rule.length")
  measure <- measures[1]
  mean_stab <- tapply(analysis[, measure]
                      , list(support = analysis$support
                             , alpha = analysis$alpha_paths
                             , bins = analysis$disc_path_bins
                             , func = analysis$score_func
                             , weights = analysis$weighting)
                      , mean)
  
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
  analysis_groups$mean_rank <- colMeans(t(apply(-analysis_values, 1, rank)))
  
  plots_b[[ds]] <- ggplot(
    data = analysis_groups
    , aes(y = mean_rank # I(1-mean)
          , x = id
          , shape = alpha
          , colour = func
          , size = support)) +
    geom_point() + 
    scale_color_grey() +
    scale_size_discrete(range = c(2, 4)) +
    theme_bw() +
    facet_grid(weights~bins)

  analysis_groups_4 <- analysis_groups %>% filter(bins == 4)
  
  analysis_groups_8 <- analysis_groups %>% filter(bins == 8)
  
  
  analysis_bins_check <- inner_join(analysis_groups_4, analysis_groups_8
                                    , by = c("support", "alpha", "func", "weights")
                                    , suffix = c("_4", "_8")) %>%
    select(-bins_4, -bins_8)
  
  bins_check[[ds]] <- wilcox.test(analysis_bins_check$mean_4
                                  , analysis_bins_check$mean_8
                                  , paired = TRUE)
  
  # we know we only need to keep bins_4
  
  analysis_groups_chisq <- analysis_groups_4 %>% filter(weights == "chisq")
  analysis_groups_nothing <- analysis_groups_4 %>% filter(weights == "nothing")
  
  analysis_weights_check <- inner_join(analysis_groups_chisq, analysis_groups_nothing
                                    , by = c("support", "alpha", "func")
                                    , suffix = c("_nothing", "_chisq")) %>%
    select(-weights_chisq, -weights_nothing)
  
  weights_check[[ds]] <- wilcox.test(analysis_weights_check$mean_chisq
                                     , analysis_weights_check$mean_nothing
                                     , paired = TRUE)
  
  plots_w[[ds]] <- ggplot(
    data = analysis_groups_4
    , aes(y = I(1-mean)
          , x = id
          , shape = alpha
          , colour = func
          , size = support)) +
    geom_point() + 
    scale_color_grey() +
    scale_size_discrete(range = c(2, 4)) +
    theme_bw() +
    facet_grid(.~weights)
  
  analysis_included_values <- analysis_values[, analysis_groups_chisq$id]
  
  analysis_groups_chisq$mean <- apply(analysis_included_values, 2, mean)
  analysis_groups_chisq$mean_rank <- colMeans(t(apply(-analysis_included_values, 1, rank)))
  
  plots_r[[ds]] <- ggplot(
    data = analysis_groups_chisq
    , aes(y = I(1-mean)
          , x = id
          , shape = alpha
          , colour = func
          , size = support)) +
    geom_point() + 
    scale_color_grey() +
    scale_size_discrete(range = c(2, 4)) +
    theme_bw()
  
  
}

analysis_groups_chisq

bins_check
weights_check





rank(-analysis_values[1, ])

analysis_simp <- analysis %>% select(support, alpha_paths, disc_path_bins, score_func, weighting, stability.tt.)