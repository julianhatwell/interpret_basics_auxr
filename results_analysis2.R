source("collect_results_from_files.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tikzDevice)
library(lme4)
library(stargazer)

analysis_allrand <- with(main_results
                         , main_results[
                           target_prec == 0.95 &
                             (
                               result_set == "anchors" | 
                                 (
                                   support %in% c(0.02, 0.05) & # check when comparing support
                                     alpha_paths == 0.5 & 
                                     alpha_scores == 0.75 
                                 )
                             )
                           , ])

datasets <- unique(analysis_allrand$datasetname)

# test stats
tt.stats <- data.frame(
  dataset = rep(NA, length(datasets))
  , ma = rep(NA, length(datasets))
  , mc2 = rep(NA, length(datasets))
  , mdiffac2 = rep(NA, length(datasets))
  , wlxac2 = rep(NA, length(datasets))
  , wlxac2p = rep(NA, length(datasets))
  
  , mc5 = rep(NA, length(datasets))
  , mdiffac5 = rep(NA, length(datasets))
  , wlxac5 = rep(NA, length(datasets))
  , wlxac5p = rep(NA, length(datasets))
  , N = rep(NA, length(datasets))
)

measures <- c("stability.tt.", "excl.cov.tt.", "rule.length")
measure <- measures[3]

for (d in seq_along(datasets)) {
  
  tt.wide <- mutate(analysis_allrand
                              , instance_id = as.numeric(instance_id)
                              , random_state = as.numeric(as.vector(random_state))) %>% # don't want a factor
    filter(datasetname == datasets[d]) %>%
    rename(measure = !! measure) %>%
    transmute(instance_id = factor(random_state * 100000 + instance_id)
              , rset_supp = rset_supp
              , measure = measure) %>%
    spread(rset_supp, measure)
  
  wlxac2.test <- wilcox.test(tt.wide[, "CHIRPS_0.02"]
                             , tt.wide[, "anchors"]
                             , data=, paired=TRUE)
  wlxac5.test <- wilcox.test(tt.wide[, "CHIRPS_0.05"]
                             , tt.wide[, "anchors"]
                             , paired=TRUE)
  
  tt.stats <- within(tt.stats, {
    dataset[d] <- as.character(datasets[d])
    ma[d] <- mean(tt.wide[["anchors"]])
    mc2[d] <- mean(tt.wide[["CHIRPS_0.02"]])
    mdiffac2[d] <- t.test(tt.wide[, "CHIRPS_0.02"]
                          , tt.wide[, "anchors"]
                          , data=, paired=TRUE)$estimate
    wlxac2[d] <- wlxac2.test$statistic
    wlxac2p[d] <- wlxac2.test$p.value
    
    mc5[d] <- mean(tt.wide[["CHIRPS_0.05"]])
    mdiffac5[d] <- t.test(tt.wide[, "CHIRPS_0.05"]
                          , tt.wide[, "anchors"]
                          , data=, paired=TRUE)$estimate
    wlxac5[d] <- wlxac5.test$statistic
    wlxac5p[d] <- wlxac5.test$p.value
    N[d] <- nrow(tt.wide)
  })
}

methodcols <- c(1,2,3,4,11,6)
methodcols <- c(1,2,7,8,11,10)

method.stats <- tt.stats[,methodcols]
names(method.stats) <- c("dataset"
                        , "anchors"
                        , "CHIRPS"
                        , "diff"
                        , "N"
                        , "p.val")

stargazer(method.stats
          , summary = FALSE
          , rownames = FALSE)

# aggregated
tt.means <- data.frame(randst = 123:127)
tt.means <- cbind(tt.means, matrix(NA, nrow = 5, ncol = 9))
names(tt.means) = c("fold", as.character(datasets))

tt.means <- rbind(tt.means, tt.means) # double up
tt.means <- rbind(tt.means, tt.means) # and again

tt.means$measure <- rep(measures, each = 10)
methods <- c("CHIRPS_0.02", "CHIRPS_0.05")
tt.means$method <- rep(rep(methods, each = 5), 2)

for (iid in 1:5) {
  for (d in seq_along(datasets)) {
    for (measure in measures) {
      for (method in methods) {
        tt.wide <- filter(analysis_allrand, datasetname == datasets[d] &
                            random_state == iid + 122) %>%
          rename(measure = !! measure) %>%
          select(instance_id, rset_supp, measure) %>%
          spread(rset_supp, measure)
        tt.means[tt.means$fold == 122 + iid &
                   tt.means$measure == measure &
                   tt.means$method == method
                 , as.character(datasets[d])] <- 
          t.test(
            tt.wide[, method], tt.wide[, "anchors"]
            , paired=TRUE)$estimate    
      }  
    }
  }
}

# do one stargazer per measure
measure <- measures[2]
tt.ttest <- select(tt.means, 1, 11:12, 2:10)
t.test2 <- sapply(filter(tt.ttest, measure == !! measure
                         , method == methods[1]) %>%
                    select(-(1:3))
                  , t.test)

t.test5 <- sapply(filter(tt.ttest, measure == !! measure
                         , method == methods[2]) %>%
                    select(-(1:3))
                  , t.test)


diff2 <- as.numeric(unlist(t.test2["estimate", ]))
p.value2 <- as.numeric(unlist(t.test2["p.value", ]))
diff5 <- as.numeric(unlist(t.test5["estimate", ]))
p.value5 <- as.numeric(unlist(t.test5["p.value", ]))

df = rep(4, 9)
agg.method.stats <- data.frame(datasets, diff2, p.value2, diff5, p.value5, df)
stargazer(agg.method.stats
          , summary = FALSE
          , rownames = FALSE)


# for plotting
tt.means.plotting <- gather(tt.means, key = dataset, value = paired.mean.diff, -fold, -method, -measure)
tt.means.plotting$method <- gsub("CHIRPS_", "Support = ", tt.means.plotting$method)

tikz(file = "stability_diff_gg.tikz", width = 3, height = 3)
plot <- ggplot(data = tt.means.plotting[tt.means.plotting$measure == measures[1], ]
               , aes(y = paired.mean.diff
                     , x = dataset
                     , colour = factor(fold))) +
  geom_bar(stat = "identity"
           , position = position_dodge(width = 0.5)
           , fill = "#999999"
           , width = 0.05) +
  geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  scale_color_grey(start = 0.2, end = 0.3) +
  facet_grid(method~.) +
  labs(y = "Mean paired differences") +
  theme_bw() +
  theme(legend.position = "none"
        , axis.title = element_text(size = rel(0.75))
        , axis.text.x = element_text(angle = -20, vjust = 0.85)
        , strip.background = element_rect(colour = "black", fill = "white"))

print(plot)
dev.off()

tikz(file = "excl_cov_diff_gg.tikz", width = 3, height = 3)
plot <- ggplot(data = tt.means.plotting[tt.means.plotting$measure == measures[2], ]
               , aes(y = paired.mean.diff
                     , x = dataset
                     , colour = factor(fold))) +
  geom_bar(stat = "identity"
           , position = position_dodge(width = 0.5)
           , fill = "#999999"
           , width = 0.05) +
  geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  scale_color_grey(start = 0.2, end = 0.3) +
  facet_grid(method~.) +
  labs(y = "Mean paired differences") +
  theme_bw() +
  theme(legend.position = "none"
        , axis.title = element_text(size = rel(0.75))
        , axis.text.x = element_text(angle = -20, vjust = 0.85)
        , strip.background = element_rect(colour = "black", fill = "white"))

print(plot)
dev.off()


# support results
analysis_supp <- support_results %>%
  select(datasetname, time.per.ex.0.01
         , time.per.ex.0.02, time.per.ex.0.05
         , time.per.ex.0.1, time.per.ex.0.2) %>%
  rename(x0.01 = time.per.ex.0.01
         , x0.02 = time.per.ex.0.02
         , x0.05 = time.per.ex.0.05
         , x0.1 = time.per.ex.0.1
         , x0.2 = time.per.ex.0.2) %>%
  gather(key = support, value = time_per_exp
         , x0.01, x0.02, x0.05, x0.1, x0.2) %>%
  mutate(support = as.numeric(gsub("x", "", support)))

analysis_supp <- within(analysis_supp, {
  time_per_exp[time_per_exp == 0] <- NA  
})

tikz(file = "time_support_gg.tikz", width = 3, height = 2)
plot <- ggplot(aes(y = time_per_exp
                   , x = support)
               , data = analysis_supp) +
  geom_point() +
  geom_smooth(colour = "#777777"
              , method = "gam"
              , formula = y~x
              , method.args = list(family = "Gamma")
  ) +
  labs(y = "Time per explanation\n(seconds)") +
  theme_bw() +
  theme(axis.title = element_text(size = rel(0.75)))
print(plot)
dev.off()

# rearrange the data
analysis_time <- time_results %>%
  select(-mtry, -sp.0.02.timing, precis.thresh
         , -sp.0.05.timing, -time..anch.
         , -acc..tt., -acc..anch.
         , - coka..tt., - coka..anch.) %>%
  rename(CHIRPS_0.02 = time.per.ex.0.02
         , CHIRPS_0.05 = time.per.ex.0.05
         , anchors = time.per.ex..anch.
         , oobe = acc..oobe.
         , max_depth = max.depth
         , min_samples = min.samples) %>%
  gather(key = rset_supp, value = time_per_exp
         , anchors, CHIRPS_0.02, CHIRPS_0.05) %>%
  mutate(datasetname = factor(tolower(datasetname))
         , rset_supp = factor(rset_supp)
         , log_time_per_exp = log(time_per_exp)) %>%
  filter(precis.thresh == 0.95) %>%
  select(-precis.thresh)

# hparams
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

rf_hparams <- sapply(analysis_time[, c("max_depth", "min_samples", "n_trees")]
                     , function(x) {tapply(x, analysis_time$datasetname, Mode)})

stargazer(rf_hparams)

# timing stats
timings <- select(analysis_time, datasetname, rset_supp, randst, log_time_per_exp) %>%
  spread(key = datasetname, value = log_time_per_exp)

timing_means <- t(with(timings
                       , sapply(select(timings, -(1:2))
                                , function(x) {tapply(x, list(rset_supp), FUN = mean)})))

timing_sds <- t(with(timings
                     , sapply(select(timings, -(1:2))
                              , function(x) {tapply(x, list(rset_supp), FUN = sd)})))

timing_lower <- data.frame(exp(timing_means - timing_sds)) %>%
  mutate(dataset = as.character(datasets)) %>%
  gather(key = rset_supp, value = time_per_exp_lwr, -dataset)

timing_upper <- data.frame(exp(timing_means + timing_sds)) %>%
  mutate(dataset = as.character(datasets)) %>%
  gather(key = rset_supp, value = time_per_exp_upr, -dataset)

timing_means <- data.frame(exp(timing_means)) %>%
  mutate(dataset = as.character(datasets)) %>%
  gather(key = rset_supp, value = time_per_exp, -dataset)

timing_means <- merge(timing_means, merge(timing_lower, timing_upper))
# plotting correction
timing_means$rset_supp <- gsub("a", "A", gsub("_", " ", timing_means$rset_supp))

tikz(file = "time_compare_gg.tikz", width = 3, height = 3)
plot <- ggplot(data = timing_means
       , aes(y = time_per_exp
             , ymin = time_per_exp_lwr
             , ymax = time_per_exp_upr
             , x = dataset
             , colour = rset_supp
             , linetype = rset_supp)) +
  geom_pointrange(position = position_dodge(width = 0.25)
                  , stat = "identity"
                  , size = 0.25) +
  labs(y = "Time per explanation\n(seconds)"
       , colour = "Method"
       , linetype = "Method") +
  scale_color_grey(start = 0.7, end = 0.3) +
  theme_bw() +
  theme(axis.title = element_text(size = rel(0.75))
      , axis.text.x = element_text(angle = -20, vjust = 0.85)
      , legend.title = element_blank()
      , legend.text = element_text(size = rel(0.5))
      , legend.position = "top")

print(plot)
dev.off()

time.model2 <- glmer(time_per_exp ~ 
                       rset_supp + 
                       (1 + rset_supp|datasetname)
                     
                     , data = analysis_time
                     , family = "Gamma"
                     , control = glmerControl(optimizer="bobyqa")
                     , start = list(fixef = c(0.37467
                                              , -0.02748
                                              , 0.54133))
)

plot(time.model2)
summary(time.model2)

y_haty <- data.frame(actual = analysis_time$log_time_per_exp
                     , fitted = log(predict(time.model2, type = "response"))
                     , rset_supp = analysis_time$rset_supp
                     , dataset = analysis_time$datasetname)

tikz(file = "time_model_gg.tikz", width = 3, height = 2)
plot <- ggplot(data = y_haty
               , aes(x = actual, y = fitted
                     , colour = rset_supp
                     , shape = rset_supp)
               ) +
  scale_shape_manual(values = c(1, 3, 4)) +
  scale_color_grey(start = 0.5, end = 0.3) +
  geom_point(size = 2.5) +
  theme_bw() +
  labs(x = "Actual time per explanation\n(log seconds)"
       , y = "Fitted time per explanation\n(log seconds)"
       , colour = "Method"
       , shape = "Method") +
  theme(axis.title = element_text(size = rel(0.75)))

print(plot)
dev.off()
