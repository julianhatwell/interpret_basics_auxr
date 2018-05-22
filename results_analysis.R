source("collect_results_from_files.R")
library(ggplot2)
library(lattice)
library(dplyr)
library(tidyr)
library(car)
library(lme4)
library(stargazer)


analysis <- with(main_results
                , main_results[
                  target_prec == 0.95 &
                    random_state == 123 & # :127 # == 123
                    (
                    result_set == "anchors" | 
                      (
                      support %in% c(0.02, 0.05) & # check when comparing support
                        alpha_paths == 0.5 & 
                        alpha_scores == 0.75 
                      )
                    )
                , ])

xyplot(stability.tt. ~ excl.cov.tt. | datasetname * rset_supp
      , alpha = 0.1
      , col = "black"
      , data = analysis
      , auto.key = TRUE)

histogram( ~f1.tt. | datasetname * rset_supp
       , data = analysis
       , auto.key = TRUE)


# stability
stability.tt.stats <- data.frame(
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
)

datasets <- unique(analysis$datasetname)
for (d in seq_along(datasets)) {

  stability.tt.wide <- filter(analysis, datasetname == datasets[d]) %>%
    select(instance_id, rset_supp, stability.tt.) %>%
    spread(rset_supp, stability.tt.)
  
  wlxac2.test <- wilcox.test(stability.tt.wide[, "CHIPS_0.02"]
                             , stability.tt.wide[, "anchors"]
                             , data=, paired=TRUE)
  wlxac5.test <- wilcox.test(stability.tt.wide[, "CHIPS_0.05"]
                             , stability.tt.wide[, "anchors"]
                             , data=, paired=TRUE)
  
  stability.tt.stats <- within(stability.tt.stats, {
    dataset[d] <- as.character(datasets[d])
    ma[d] <- mean(stability.tt.wide[["anchors"]])
    mc2[d] <- mean(stability.tt.wide[["CHIPS_0.02"]])
    mdiffac2[d] <- t.test(stability.tt.wide[, "CHIPS_0.02"]
                         , stability.tt.wide[, "anchors"]
                         , data=, paired=TRUE)$estimate
    wlxac2[d] <- wlxac2.test$statistic
    wlxac2p[d] <- wlxac2.test$p.value

    mc5[d] <- mean(stability.tt.wide[["CHIPS_0.05"]])
    mdiffac5[d] <- t.test(stability.tt.wide[, "CHIPS_0.05"]
                         , stability.tt.wide[, "anchors"]
                         , data=, paired=TRUE)$estimate
    wlxac5[d] <- wlxac5.test$statistic
    wlxac5p[d] <- wlxac5.test$p.value
  })
}

n_samples <- table(analysis$datasetname)/3
stability.tt.stats <- cbind(stability.tt.stats
    , N = as.vector(n_samples[s0.02.stats$dataset]))

s0.02.stats <- stability.tt.stats[,c(1,2,3,4,11,6)]
names(s0.02.stats) <- c("dataset"
                        , "anchors"
                        , "CHIRPS"
                        , "diff"
                        , "N"
                        , "p.val")

stargazer(s0.02.stats
         , summary = FALSE
         , rownames = FALSE)

s0.05.stats <- stability.tt.stats[,c(1,2,7,8,11,10)]
names(s0.05.stats) <- c("dataset"
                        , "anchors"
                        , "CHIRPS"
                        , "diff"
                        , "N"
                        , "p.val")

stargazer(s0.05.stats
          , summary = FALSE
          , rownames = FALSE)

# exclusive coverage
excl.cov.tt.stats <- data.frame(
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
)

datasets <- unique(analysis$datasetname)
for (d in seq_along(datasets)) {
  
  
  excl.cov.tt.wide <- filter(analysis, datasetname == datasets[d]) %>%
    select(instance_id, rset_supp, excl.cov.tt.) %>%
    spread(rset_supp, excl.cov.tt.)
  
  wlxac2.test <- wilcox.test(excl.cov.tt.wide[, "CHIPS_0.02"]
                             , excl.cov.tt.wide[, "anchors"]
                             , data=, paired=TRUE)
  wlxac5.test <- wilcox.test(excl.cov.tt.wide[, "CHIPS_0.05"]
                             , excl.cov.tt.wide[, "anchors"]
                             , data=, paired=TRUE)
  
  excl.cov.tt.stats <- within(excl.cov.tt.stats, {
    dataset[d] <- as.character(datasets[d])
    ma[d] <- mean(excl.cov.tt.wide[["anchors"]])
    mc2[d] <- mean(excl.cov.tt.wide[["CHIPS_0.02"]])
    mdiffac2[d] <- t.test(excl.cov.tt.wide[, "CHIPS_0.02"]
                          , excl.cov.tt.wide[, "anchors"]
                          , data=, paired=TRUE)$estimate
    wlxac2[d] <- wlxac2.test$statistic
    wlxac2p[d] <- wlxac2.test$p.value
    
    
    
    mc5[d] <- mean(excl.cov.tt.wide[["CHIPS_0.05"]])
    mdiffac5[d] <- t.test(excl.cov.tt.wide[, "CHIPS_0.05"]
                          , excl.cov.tt.wide[, "anchors"]
                          , data=, paired=TRUE)$estimate
    wlxac5[d] <- wlxac5.test$statistic
    wlxac5p[d] <- wlxac5.test$p.value
    
  })
  
}

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

ggplot(aes(y = time_per_exp
           , x = support)
       , data = analysis_supp) +
  geom_point(aes(colour = datasetname)) +
  geom_smooth(method = "gam"
              , formula = y~x
              , method.args = list(family = "Gamma")
              ) +
  theme_bw()

supp.null <- lmer(log(time_per_exp) ~ (1|datasetname)
      , data = analysis_supp
      )

supp.model <- lmer(log(time_per_exp) ~ lg(support) + 
                      (1|datasetname)
                   , data = analysis_supp)

anova(supp.null, supp.model)
summary(supp.model)



# timing results
# pull out the model performance stats
CHIRPS_acc <- time_results$acc..tt.
CHIRPS_ck <- time_results$coka..tt.
anch_acc <- time_results$acc..anch.
anch_ck <- time_results$coka..anch.

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
  mutate(rset_supp = factor(rset_supp)
         , model_accuracy = c(anch_acc, CHIRPS_acc, CHIRPS_acc)
         , model_kappa = c(anch_ck, CHIRPS_ck, CHIRPS_ck)
         , log_time_per_exp = log(time_per_exp)) %>%
  filter(precis.thresh == 0.95) %>%
  select(-precis.thresh)

densityplot(~time_per_exp, data = analysis_time)
densityplot(~time_per_exp, data = analysis_time)

time_aov_null <- aov(log_time_per_exp ~ datasetname, data = analysis_time)
time_aov_1 <- aov(log_time_per_exp ~ rset_supp * datasetname, data = analysis_time)
time_aov2 <- aov(log_time_per_exp ~ rset_supp * datasetname * model_kappa, data = analysis_time)

time.null <- lmer(log_time_per_exp ~ 
                     (1 + 1|datasetname)
  , data = analysis_time
  , REML = FALSE)
plot(time.null)                  

time.model1 <- lmer(log_time_per_exp ~ 
                  rset_supp + 
                  (1 + 1|datasetname)

  , data = analysis_time
  , REML = FALSE)
plot(time.model1)
summary(time.model1)

anova(time.null, time.model1)

time.model2 <- lmer(log_time_per_exp ~ 
                      rset_supp + 
                      (1 + rset_supp|datasetname)
                    
                    , data = analysis_time
                    , REML = FALSE)

plot(time.model2)
summary(time.model2)
anova(time.null, time.model2)
anova(time.model1, time.model2)
anova(time.null, time.model1, time.model2)

time.model3 <- lmer(log_time_per_exp ~ 
                      rset_supp * model_kappa +
                      (1 + rset_supp|datasetname)
                    
                    , data = analysis_time
                    , REML = FALSE)
plot(time.model3)
summary(time.model3)
anova(time.null, time.model3)
anova(time.model2, time.model3)

time.model4 <- lmer(log_time_per_exp ~ 
                      rset_supp + 
                      (1 + rset_supp|datasetname) +
                      (1|randst)
                    
                    , data = analysis_time
                    , REML = FALSE)

plot(time.model4)
summary(time.model4)
anova(time.null, time.model4)
anova(time.model2, time.model4)

time.null <- glmer(time_per_exp ~ 
                    (1|datasetname)
                  , data = analysis_time
                  , family = "Gamma")
plot(time.null)                  
summary(time.null)

time.model1 <- glmer(time_per_exp ~ 
                      rset_supp + 
                      (1|datasetname)
                    
                    , data = analysis_time
                    , family = "Gamma")
plot(time.model1)
summary(time.model1)

anova(time.null, time.model1)

time.model2 <- glmer(time_per_exp ~ 
                      rset_supp + 
                      (1 + rset_supp|datasetname)
                    
                    , data = analysis_time
                    , family = "Gamma")

plot(time.model2)
summary(time.model2)
anova(time.null, time.model2)
anova(time.model1, time.model2)
anova(time.null, time.model1, time.model2)

time.model3 <- glmer(time_per_exp ~ 
                      rset_supp * model_kappa +
                      (1 + rset_supp|datasetname)
                    
                    , data = analysis_time
                    , family = "Gamma")
plot(time.model3)
summary(time.model3)
anova(time.null, time.model3)
anova(time.model2, time.model3)

time.model4 <- glmer(time_per_exp ~ 
                      rset_supp + 
                      (1 + rset_supp|datasetname) +
                      (1|randst)
                    
                    , data = analysis_time
                    , family = "Gamma")

plot(time.model4)
summary(time.model4)
anova(time.null, time.model4)
anova(time.model2, time.model4)


# f1 score bc transform
analysis_f1 <- with(main_results
                 , main_results[
                   target_prec == 0.95 & 
                   random_state %in% 123:126 &             
                   (
                     result_set == "anchors" | 
                       (
                         support %in% c(0.02, 0.05) & # check when comparing support
                           alpha_paths == 0.5 & 
                           alpha_scores == 0.75 
                       )
                     )
                   , ])

for (d in seq_along(datasets)) {
  
  if(as.character(datasets[d]) != "rcdv") {
    f1_by_db <- filter(analysis_f1, datasetname == datasets[d]) %>%
      select(instance_id, datasetname, rset_supp, f1.tt., pred.class.label
             , model_acc, model_ck, random_state, support, target_prec)


    f1_ac <- asin(f1_by_db$f1.tt.) + 0.75
    lambda <- coef(powerTransform(1/f1_ac))
    
    f1_by_db$f1_trsf <- bcPower(f1_ac, lambda)
    
    if (d == 1) {
      out_f1 <- f1_by_db
    } else {
      out_f1 <- rbind(out_f1, f1_by_db)
    }
  }
}

histogram( ~f1_trsf | datasetname * rset_supp
           , data = out_f1
           , scales = list(relation = "free")
           , auto.key = TRUE)


precdiff.null = lmer(f1_trsf ~ 
                       (1 + rset_supp|instance_id) +
                       (1|datasetname) +
                       (1|random_state)
                     ,
                     data=out_f1,
                     REML=FALSE)

precdiff.model = lmer(f1_trsf ~ rset_supp +
                      (1 + rset_supp|instance_id) +
                      (1|datasetname) +
                      (1|random_state)
                      ,
                        data=out_f1,
                        REML=FALSE)

precdiff.model
summary(precdiff.model)
plot(precdiff.model)
anova(precdiff.null, precdiff.model)


precdiff.model2 = lmer(f1_trsf ~ rset_supp * model_ck +
                        (1 + rset_supp|instance_id) +
                        (1|datasetname) +
                        (1|random_state)
                      ,
                      data=out_f1,
                      REML=FALSE)

precdiff.model2
summary(precdiff.model2)
plot(precdiff.model2)
anova(precdiff.model, precdiff.model2)
anova(precdiff.null, precdiff.model2)
anova(precdiff.null, precdiff.model, precdiff.model2)
