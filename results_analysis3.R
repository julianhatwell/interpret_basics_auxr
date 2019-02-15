library(dplyr)


analysis <- main_results %>% 
  filter(dataset == "nursery_samp") %>%
  mutate(support = factor(support)
            , alpha_paths = factor(alpha_paths)
            , disc_path_bins = factor(disc_path_bins)
            , score_func = factor(score_func)
            , weighting = weighting)
measures <- c("stability.tt.", "xcoverage.tt.", "rule.length")
measure <- measures[1]
tapply(analysis[, measure]
       , list(analysis$support
              , analysis$alpha_paths
              , analysis$disc_path_bins
              , analysis$score_func
              , analysis$weighting)
       , mean)

xtabs(~weighting+score_func, data = analysis)
