source("collect_results_from_files_addt.R")
library(ggplot2)
library(lattice)
library(dplyr)
library(tidyr)
library(car)
library(lme4)
library(stargazer)
library(tikzDevice)
library(lmomco)
library(bootstrap)

unpretty <- function(r) {
  strsplit(r, " AND ")
}

with(main_results, main_results[datasetname == "adult" & addt == 0, "instance_id"])
with(main_results, main_results[datasetname == "adult" & instance_id == 606, ])





# one random state set
analysis <- with(main_results
                 , main_results[
                     datasetname == "adult"
                     , ])

datasets <- unique(analysis$datasetname)

xyplot(stability.tt. ~ excl.cov.tt. | datasetname
       , groups = pred.class.label
      , alpha = 0.1

      , data = analysis
      , auto.key = list(columns = 6))

histogram( ~f1.tt. | datasetname
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

n_samples <- table(analysis$datasetname)/3
excl.cov.tt.stats <- cbind(excl.cov.tt.stats
                            , N = as.vector(n_samples[s0.02.stats$dataset]))

s0.02.stats <- excl.cov.tt.stats[,c(1,2,3,4,11,6)]
names(s0.02.stats) <- c("dataset"
                        , "anchors"
                        , "CHIRPS"
                        , "diff"
                        , "N"
                        , "p.val")

stargazer(s0.02.stats
          , summary = FALSE
          , rownames = FALSE)

s0.05.stats <- excl.cov.tt.stats[,c(1,2,7,8,11,10)]
names(s0.05.stats) <- c("dataset"
                        , "anchors"
                        , "CHIRPS"
                        , "diff"
                        , "N"
                        , "p.val")

stargazer(s0.05.stats
          , summary = FALSE
          , rownames = FALSE)

# aggregated
stability_means <- with(analysis_allrand, tapply(stability.tt.
                              , list(instance_id, random_state, rset_supp, datasetname)
                              , FUN=mean))
stability.tt.means <- data.frame(randst = 1:5)
stability.tt.means <- cbind(stability.tt.means, matrix(NA, nrow = 5, ncol = 9))
names(stability.tt.means) = c("fold", as.character(datasets))

# aggregate for supp 0.02
for (iid in 1:5) {
  for (d in seq_along(datasets)) {
    stability.tt.wide <- filter(analysis_allrand, datasetname == datasets[d] & random_state == iid + 122) %>%
      select(instance_id, rset_supp, stability.tt.) %>%
      spread(rset_supp, stability.tt.)
    
    stability.tt.means[[iid, as.character(datasets[d])]] <- t.test(
      stability.tt.wide[, "CHIPS_0.02"]
     , stability.tt.wide[, "anchors"]
     , paired=TRUE)$estimate
  }
}

res <- sapply(select(stability.tt.means, -1), t.test)

diff2 <- as.numeric(unlist(res["estimate", ]))
p.value2 <- as.numeric(unlist(res["p.value", ]))

# aggregate for supp 0.02
for (iid in 1:5) {
  for (d in seq_along(datasets)) {
    stability.tt.wide <- filter(analysis_allrand, datasetname == datasets[d] & random_state == iid + 122) %>%
      select(instance_id, rset_supp, stability.tt.) %>%
      spread(rset_supp, stability.tt.)
    
    stability.tt.means[[iid, as.character(datasets[d])]] <- t.test(
      stability.tt.wide[, "CHIPS_0.05"]
      , stability.tt.wide[, "anchors"]
      , paired=TRUE)$estimate
  }
}

res <- sapply(select(stability.tt.means, -1), t.test)

diff5 <- as.numeric(unlist(res["estimate", ]))
p.value5 <- as.numeric(unlist(res["p.value", ]))

df = rep(4, 9)
stargazer(data.frame(datasets, diff2, p.value2, diff5, p.value5, df)
          , summary = FALSE
          , rownames = FALSE)


excl.cov_means <- with(analysis_allrand, tapply(excl.cov.tt.
                                                 , list(instance_id, random_state, rset_supp, datasetname)
                                                 , FUN=mean))
excl.cov.tt.means <- data.frame(randst = 1:5)
excl.cov.tt.means <- cbind(excl.cov.tt.means, matrix(NA, nrow = 5, ncol = 9))
names(excl.cov.tt.means) = c("fold", as.character(datasets))

# aggregate for supp 0.02
for (iid in 1:5) {
  for (d in seq_along(datasets)) {
    excl.cov.tt.wide <- filter(analysis_allrand, datasetname == datasets[d] & random_state == iid + 122) %>%
      select(instance_id, rset_supp, excl.cov.tt.) %>%
      spread(rset_supp, excl.cov.tt.)
    
    excl.cov.tt.means[[iid, as.character(datasets[d])]] <- t.test(
      excl.cov.tt.wide[, "CHIPS_0.02"]
      , excl.cov.tt.wide[, "anchors"]
      , paired=TRUE)$estimate
  }
}

res <- sapply(select(excl.cov.tt.means, -1), t.test)

diff2 <- as.numeric(unlist(res["estimate", ]))
p.value2 <- as.numeric(unlist(res["p.value", ]))

# aggregate for supp 0.02
for (iid in 1:5) {
  for (d in seq_along(datasets)) {
    excl.cov.tt.wide <- filter(analysis_allrand, datasetname == datasets[d] & random_state == iid + 122) %>%
      select(instance_id, rset_supp, excl.cov.tt.) %>%
      spread(rset_supp, excl.cov.tt.)
    
    excl.cov.tt.means[[iid, as.character(datasets[d])]] <- t.test(
      excl.cov.tt.wide[, "CHIPS_0.05"]
      , excl.cov.tt.wide[, "anchors"]
      , paired=TRUE)$estimate
  }
}

res <- sapply(select(excl.cov.tt.means, -1), t.test)

diff5 <- as.numeric(unlist(res["estimate", ]))
p.value5 <- as.numeric(unlist(res["p.value", ]))

df = rep(4, 9)
stargazer(data.frame(datasets, diff2, p.value2, diff5, p.value5, df)
          , summary = FALSE
          , rownames = FALSE)


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
  geom_smooth(method = "gam"
              , formula = y~x
              , method.args = list(family = "Gamma")
              ) +
  labs(y = "Time per explanation\n(seconds)") +
  theme(axis.text.x = element_text(size = rel(0.52))
        , axis.text.y = element_text(size = rel(0.25))) +
  theme_bw()
print(plot)
dev.off()



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
  mutate(datasetname = factor(tolower(datasetname))
         , rset_supp = factor(rset_supp)
         , model_accuracy = c(anch_acc, CHIRPS_acc, CHIRPS_acc)
         , model_kappa = c(anch_ck, CHIRPS_ck, CHIRPS_ck)
         , log_time_per_exp = log(time_per_exp)) %>%
  filter(precis.thresh == 0.95) %>%
  select(-precis.thresh)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

rf_hparams <- sapply(analysis_time[, c("max_depth", "min_samples", "n_trees")]
       , function(x) {tapply(x, analysis_time$datasetname, Mode)})

stargazer(rf_hparams)

densityplot(~time_per_exp, data = analysis_time)
densityplot(~log_time_per_exp, data = analysis_time)


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

ggplot(data = timing_means
       , aes(y = time_per_exp
             , ymin = time_per_exp_lwr
             , ymax = time_per_exp_upr
             , x = dataset
             , colour = rset_supp
             , linetype = rset_supp)) +
  geom_pointrange(position = position_dodge(width = 0.5)
                , stat = "identity"
                , size = 1) +
  #geom_point(aes(y = time_per_exp), position = position_dodge(0.9)) +
  labs(y = "Time per explanation\n(seconds)") +
  theme(axis.text.x = element_text(size = rel(0.52), angle=30, vjust = 0.5)
        , axis.text.y = element_text(size = rel(0.25))) +
  scale_color_grey() +
  theme_bw()
  


timing_harm <- t(with(timings
    , sapply(select(timings, -(1:2))
             , function(x) {tapply(x, list(rset_supp), FUN = harmonic.mean)})))
timing_harm <- as.data.frame(matrix(sapply(timing_harm, function(x) {x[[1]]}), ncol = 3)
                             , row.names = as.character(datasets))
names(timing_harm) <- c("anchors", "CHIRPS_0.02", "CHIRPS_0.05")

timing_harmsd <- as.data.frame(matrix(NA, nrow = 9, ncol = 3)
                               , row.names = as.character(datasets))
names(timing_harmsd) <- c("anchors", "CHIRPS_0.02", "CHIRPS_0.05")

for (d in seq_along(datasets)) {
  for(rs in unique(analysis_time$rset_supp)) {
    
    
    timing_harmsd[d, rs] <- with(analysis_time
                            , jackknife(time_per_exp[datasetname == as.character(datasets[d]) &
                                          rset_supp == rs]
                                       , theta = function(x) {
                                                 harmonic.mean(x)[[1]]
                                         }))[[1]]
  }
}

time_aov_null <- aov(log_time_per_exp ~ datasetname, data = analysis_time)
time_aov_1 <- aov(log_time_per_exp ~ rset_supp * datasetname, data = analysis_time)



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
                    , family = "Gamma"
                    , control = glmerControl(optimizer="bobyqa")
                    , start = list(fixef = c(0.37467
                                   , -0.02748
                                   , 0.54133))
                    )

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


#rule lengths
dataset <- analysis[analysis$datasetname == "cardiotography", ]

longest_rule <- max(dataset$rule.length)

tab <- table(dataset$rule.length)
dimnames(tab)[[1]] <- as.character(0:(longest_rule - 1))

plot(goodfit(tab, type = "poisson"), shade = TRUE)
plot(goodfit(tab, type = "nbinomial"), shade = TRUE)

fit <- as.data.frame(tab, stringsAsFactors = FALSE)
colnames(fit) <- c("rule.length", "freq")
fit$rule.length <- as.integer(fit$rule.length)
all_length <- data.frame(rule.length = 0:(longest_rule-1))
fit <- left_join(all_length, fit, by = "rule.length")
fit$freq <- ifelse(is.na(fit$freq), 0, fit$freq)

(lambda <- weighted.mean(as.numeric(fit$rule.length), w = fit$freq))

phat <- dpois(0 : (longest_rule-1), lambda = 2)
exp <- sum(fit$freq) * phat
chisq <- (fit$freq - exp)^2 / exp

GOF <- data.frame(fit, phat, exp, chisq)
GOF

sum(chisq)  # chi-square value
pchisq(sum(chisq), df = nrow(tab) - 2, lower.tail = FALSE)
summary(goodfit(GOF[, 2:1], type="nbinomial"))
summary(goodfit(GOF[, 2:1], type="poisson"))
rootogram(goodfit(GOF[, 2:1], type="nbinomial"), shade = TRUE)



dataset <- analysis[analysis$datasetname == "cardiotography", ]
longest_rule <- max(dataset$rule.length)


tabanch <- with(dataset
     , table(rule.length[rset_supp == "anchors"]))
tabch2 <- with(dataset
             , table(rule.length[rset_supp == "CHIRPS_0.02"]))
tabch5 <- with(dataset
               , table(rule.length[rset_supp == "CHIRPS_0.05"]))

anch <- as.data.frame(tabanch, stringsAsFactors = FALSE)
names(anch) <- c("rule.length", "Freq")
anch$rule.length <- as.integer(anch$rule.length)

ch2 <- as.data.frame(tabch2, stringsAsFactors = FALSE)
names(ch2) <- c("rule.length", "Freq")
ch2$rule.length <- as.integer(ch2$rule.length)
ch5 <- as.data.frame(tabch5, stringsAsFactors = FALSE)
names(ch5) <- c("rule.length", "Freq")
ch5$rule.length <- as.integer(ch5$rule.length)


all_length <- data.frame(rule.length = 0:(longest_rule))


fit <- left_join(all_length, anch, by = "rule.length") %>%
  left_join(ch2, by = "rule.length", suffix = c("anchors", "")) %>%
  left_join(ch5, by = "rule.length", suffix = c("CHIRPS_0.02", "CHIRPS_0.05"))
names(fit) <- gsub("Freq", "", names(fit))
fit <- as.data.frame(sapply(fit, function(x) {
  ifelse(is.na(x), 0, x)
}))

fit_gather <- gather(fit, rset_supp, Freq, -1)

rl <- matrix(rep(NA, 6), nrow = 2)
dimnames(rl) <- list(c("mean", "sd")
        , c("anchors", "CHIRPS_0.02", "CHIRPS_0.05"))
for (rs in c("anchors", "CHIRPS_0.02", "CHIRPS_0.05")) {
  freqs <- expand.dft(fit_gather[fit_gather$rset_supp == rs, c(1, 3)])
  
  rl[,rs] <- apply(freqs, 2
         , FUN = function(z) c(mean = mean(z), var = var(z)))
}
rl

sapply(fit, function(x) {
  x <- data.frame(Freq = x)
  
})



expand.dft(fit_gather)

sapply(fit[, 2:4], weighted.mean, w = fit$rule.length)
densityplot(~anchors+CHIRPS_0.02+CHIRPS_0.05, data = fit)

summary(glmer(rule.length~rset_supp +
                (1|datasetname)
              , data = analysis_allrand
              , family = negative.binomial(2)))

