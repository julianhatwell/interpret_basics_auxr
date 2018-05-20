source("collect_results_from_files.R")
library(ggplot2)
library(lattice)

analysis <- with(main_results
                , main_results[result_set == "anchors" | 
                              (
                                support %in% c(0.02, 0.05) & # check when comparing support
                                alpha_paths == 0.5 & 
                                alpha_scores == 0.75 
                              ) &
                               target_prec == 0.95 &
                               random_state %in% 123 # :127 # == 123
                              , ])

xyplot(test_cons ~ test_covx | datasetname * rset_supp
      #, groups = result_set
      , alpha = 0.05
      , col = "black"
      , data = analysis
      , auto.key = TRUE)

histogram( ~f1.tt. | datasetname * rset_supp
       , data = analysis
       , auto.key = TRUE)


test_cons_wide <- select(analysis, instance_id, rset_supp, test_cons) %>%
  spread(rset_supp, test_cons)

head(spread(analysis, key=rset_supp, value=test_cons))

datasets <- unique(analysis$datasetname)

test_cons_stats <- data.frame(
  ma = rep(NA, length(datasets))
  , mc2 = rep(NA, length(datasets))
  , mc5 = rep(NA, length(datasets))
  , mdiff2 = rep(NA, length(datasets))
  , mdiff5 = rep(NA, length(datasets))
  , w2 = rep(NA, length(datasets))
  , w5 = rep(NA, length(datasets))
)
rownames(analysis_stats) <- datasets

for (d in datasets) {
  test_cons_stats[d, "ma"] <- mean(analysis$test_cons[analysis$result_set=="anchors"])
  test_cons_stats[d, "mc2"] <- mean(analysis$test_cons[analysis$result_set=="CHIPS_0.02"])
  test_cons_stats[d, "mc5"] <- mean(analysis$test_cons[analysis$result_set=="CHIPS_0.05"])
  
    cons <- analysis[analysis$datasetname == d
                   , c("result_set", "instance_id", "test_cons")]
  print(c("anchors", mean(prec$precision.tt.[prec$result_set=="anchors"])))
  print(c("CHIPS", mean(prec$precision.tt.[prec$result_set=="CHIPS"])))
  print(wilcox.test(precision.tt. ~ result_set, paired = TRUE, data = f1[order(prec$instance_id), ]))
  print(t.test(precision.tt. ~ result_set, paired = TRUE, data = f1[order(prec$instance_id), ]))
}

for (d in datasets) {
  print(d)
  f1 <- analysis[analysis$datasetname == d
                   , c("result_set", "instance_id", "f1.tt.")]
  print(c("anchors", mean(f1$f1.tt.[f1$result_set=="anchors"])))
  print(c("CHIPS", mean(f1$f1.tt.[f1$result_set=="CHIPS"])))
  print(wilcox.test(f1.tt. ~ result_set, paired = TRUE, data = f1[order(f1$instance_id), ]))
  print(t.test(f1.tt. ~ result_set, paired = TRUE, data = f1[order(f1$instance_id), ]))
}

analysis <- with(main_results
                 , main_results[support == 0.02 &
                                  alpha_paths == 0.5 &
                                  alpha_scores == 0.75 &
                                  target_prec == 0.95 &
                                  random_state == 123, ])

precdiff.model = lmer(precision.tt. ~ result_set +
                          (1 + result_set|instance_id) +
                          (1 + result_set|datasetname),
                        data=analysis,
                        REML=FALSE)
precdiff.model2 = lmer(precision.tt. ~ result_set + model_ck +
                        (1 + result_set|instance_id) +
                        (1 + result_set|datasetname),
                      data=analysis,
                      REML=FALSE)
precdiff.null = lmer(precision.tt. ~ 
                       (1 + result_set|instance_id) +
                       (1 + result_set|datasetname),
                      data=analysis,
                      REML=FALSE)

precdiff.model
summary(precdiff.model)
anova(precdiff.null, precdiff.model)
precdiff.model2
summary(precdiff.model2)
anova(precdiff.model, precdiff.model2)
