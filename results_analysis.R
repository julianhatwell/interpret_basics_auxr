source("collect_results_from_files.R")
library(ggplot2)
library(lattice)

analysis <- with(main_results
                , main_results[support == 0.02 &
                                 alpha_paths == 0.5 &
                                 alpha_scores == 0.75 &
                                 target_prec == 0.95 &
                                 random_state == 123, ])

xyplot(precision.tt. ~ recall.tt. | datasetname * result_set
      #, groups = result_set
      , alpha = 0.5
      , data = analysis
      , auto.key = TRUE)

histogram( ~f1.tt. | datasetname * result_set
       , data = analysis
       , auto.key = TRUE)
bwplot(f1.tt. ~ result_set | datasetname
       , data = analysis
       , auto.key = TRUE)

datasets <- unique(analysis$datasetname)
for (d in datasets) {
  print(d)
  prec <- analysis[analysis$datasetname == d
                   , c("result_set", "instance_id", "precision.tt.")]
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
