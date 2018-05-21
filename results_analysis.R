source("collect_results_from_files.R")
library(ggplot2)
library(lattice)
library(car)
library(lme4)


analysis <- with(main_results
                , main_results[
                  target_prec >= 0.95 &
                    #random_state == 123 & # :127 # == 123
                    (
                    result_set == "anchors" | 
                      (
                      support %in% c(0.02, 0.05) & # check when comparing support
                        alpha_paths == 0.5 & 
                        alpha_scores == 0.75 
                      )
                    )
                , ])

xyplot(test_cons ~ test_covx | datasetname * rset_supp
      , alpha = 0.1
      , col = "black"
      , data = analysis
      , auto.key = TRUE)

histogram( ~f1.tt. | datasetname * rset_supp
       , data = analysis
       , auto.key = TRUE)


# consistency
test_cons_stats <- data.frame(
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
  
  
  test_cons_wide <- filter(analysis, datasetname == datasets[d]) %>%
    select(instance_id, rset_supp, test_cons) %>%
    spread(rset_supp, test_cons)
  
  wlxac2.test <- wilcox.test(test_cons_wide[, "CHIPS_0.02"]
                             , test_cons_wide[, "anchors"]
                             , data=, paired=TRUE)
  wlxac5.test <- wilcox.test(test_cons_wide[, "CHIPS_0.05"]
                             , test_cons_wide[, "anchors"]
                             , data=, paired=TRUE)
  
  test_cons_stats <- within(test_cons_stats, {
    dataset[d] <- as.character(datasets[d])
    ma[d] <- mean(test_cons_wide[["anchors"]])
    mc2[d] <- mean(test_cons_wide[["CHIPS_0.02"]])
    mdiffac2[d] <- t.test(test_cons_wide[, "CHIPS_0.02"]
                         , test_cons_wide[, "anchors"]
                         , data=, paired=TRUE)$estimate
    wlxac2[d] <- wlxac2.test$statistic
    wlxac2p[d] <- wlxac2.test$p.value
    
    
    
    mc5[d] <- mean(test_cons_wide[["CHIPS_0.05"]])
    mdiffac5[d] <- t.test(test_cons_wide[, "CHIPS_0.05"]
                         , test_cons_wide[, "anchors"]
                         , data=, paired=TRUE)$estimate
    wlxac5[d] <- wlxac5.test$statistic
    wlxac5p[d] <- wlxac5.test$p.value
    
  })
  
}

# exclusive coverage
test_covx_stats <- data.frame(
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
  
  
  test_covx_wide <- filter(analysis, datasetname == datasets[d]) %>%
    select(instance_id, rset_supp, test_covx) %>%
    spread(rset_supp, test_covx)
  
  wlxac2.test <- wilcox.test(test_covx_wide[, "CHIPS_0.02"]
                             , test_covx_wide[, "anchors"]
                             , data=, paired=TRUE)
  wlxac5.test <- wilcox.test(test_covx_wide[, "CHIPS_0.05"]
                             , test_covx_wide[, "anchors"]
                             , data=, paired=TRUE)
  
  test_covx_stats <- within(test_covx_stats, {
    dataset[d] <- as.character(datasets[d])
    ma[d] <- mean(test_covx_wide[["anchors"]])
    mc2[d] <- mean(test_covx_wide[["CHIPS_0.02"]])
    mdiffac2[d] <- t.test(test_covx_wide[, "CHIPS_0.02"]
                          , test_covx_wide[, "anchors"]
                          , data=, paired=TRUE)$estimate
    wlxac2[d] <- wlxac2.test$statistic
    wlxac2p[d] <- wlxac2.test$p.value
    
    
    
    mc5[d] <- mean(test_covx_wide[["CHIPS_0.05"]])
    mdiffac5[d] <- t.test(test_covx_wide[, "CHIPS_0.05"]
                          , test_covx_wide[, "anchors"]
                          , data=, paired=TRUE)$estimate
    wlxac5[d] <- wlxac5.test$statistic
    wlxac5p[d] <- wlxac5.test$p.value
    
  })
  
}


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
