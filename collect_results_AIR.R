# comparative analysis
rm(list=ls())

library(dplyr)
library(tidyr)
library(PMCMRplus)
library(cowplot)
library(rlang)
options(max.print=20*72)
algo_variant <- "gbt-HIPS" # there is a hardcoded version that must be changed. Find and replace.
algo_input_csv <- c("Anchors", "BRL", "CHIRPS", "defragTrees", "inTrees", "lore")
patt_input_csv <- paste0("(", paste(algo_input_csv, collapse = ")|("), ")")
algorithms <- c("Anchors", "BRL", algo_variant, "defragTrees", "inTrees", "LORE")
algorithms <- sort(algorithms)

source("data_files_mgmt_gbt.R")

resfilesdirs <- paste0(project_dir, datasetnames, pathsep)

first_comp <- TRUE
first_comp_summ <- TRUE

for (i in seq_along(resfilesdirs)) {
  datasetname <- datasetnames[i]
  filepath <- normalizePath(file.path(resfilesdirs[i]))
  filenames <- dir(filepath)
  for (fn in filenames) {
    if (fn == "y_test.csv") {
      test_set_size[datasetname] <- nrow(read.csv(normalizePath(file.path(filepath, fn))
                                                  , stringsAsFactors = FALSE))
    }
    if (fn == "y_train.csv") {
      train_set_size[datasetname] <- nrow(read.csv(normalizePath(file.path(filepath, fn))
                                                   , stringsAsFactors = FALSE))
    }
    if (grepl(patt_input_csv, fn)) {
      # load a sheet
      results <- read.csv(normalizePath(file.path(filepath, fn))
                          , stringsAsFactors = FALSE)
      
      if (!grepl("summary", fn)) {
        if (first_comp == TRUE) {
          comp_results <- results
          first_comp <- FALSE
        } else {
          names(results) <- names(comp_results)
          comp_results <- rbind(comp_results, results)
        }
      } else { # summary
        if (!any(grepl("model_kappa", names(results)))) { # no kappa model column
          results <- results %>%
            mutate(model_kappa = NA)
        }
        if (!any(grepl("proxy_kappa", names(results)))) { # no kappa model column
          results <- results %>%
            mutate(proxy_kappa = NA)
        }
        if (!any(grepl("median_rule_cascade", names(results)))) { # no med rc column
          results <- results %>%
            mutate(median_rule_cascade = NA)
        }
        if (!any(grepl("^fidelity", names(results)))) { # no kappa model column
          results <- results %>%
            mutate(fidelity = NA)
        }
        if (!any(grepl("^proxy_performance", names(results)))) { # no med rc column
          results <- results %>%
            mutate(proxy_performance = NA)
        }
        results <- results %>%
          dplyr::select(X, dataset_name, algorithm, n_instances, n_rules, n_rules_used
                        , median_rule_cascade, mean_rule_cascade, sd_rule_cascade
                        , mean_rulelen, sd_rulelen, begin_time, completion_time
                        , forest_performance, sd_forest_performance
                        , model_kappa, proxy_performance, sd_proxy_performance
                        , proxy_kappa, fidelity, sd_fidelity)
        if (first_comp_summ == TRUE) {
          comp_summ_results <- results
          first_comp_summ <- FALSE
        } else {
          comp_summ_results <- rbind(comp_summ_results, results)
        }
      }
    }
  }
}
names(comp_results) <- sub("_name", "", names(comp_results))
comp_results$algorithm <- ifelse(comp_results$algorithm == "greedy_stab", algo_variant, comp_results$algorithm)
comp_results$algorithm <- ifelse(comp_results$algorithm == "lore", "LORE", comp_results$algorithm)
# comp_results$algorithm <- ifelse(comp_results$algorithm == "GTB_HIPS", "gbt-HIPS", comp_results$algorithm)
comp_results$pretty.rule[is.na(comp_results$pretty.rule)] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == ""] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == "[]"] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == "X[,1]==X[,1]"] <- "{default}"
comp_results$dataset[comp_results$dataset == "lending_tiny_samp"] <- "lending"
comp_results$dataset[comp_results$dataset == "bankmark"] <- "bank"

# adjust 0-length rules to be valued as 1, {default}
comp_results$rule.length <- ifelse(comp_results$pretty.rule == "{default}"
                                   , 1
                                   , comp_results$rule.length)

# calculate weighted excl.cov
comp_results <- within(comp_results, {
  wxcoverage.tr. <- xcoverage.tr. * (nci.tr./(ci.tr. + nci.tr.))
  wxcoverage.tt. <- xcoverage.tt. * (nci.tt./(ci.tt. + nci.tt.))
})

# sort out the lengths of defragTrees rules
relength_pretty_rule <- function(pretty_rule) {
  length(unique(sapply(strsplit(pretty_rule, " AND ")[[1]]
                       , function(x) {
                         string <- substr(x, 1, regexpr(" ", x) - 1)
                         underscore <- regexpr("_", x)
                         if (underscore < 0) underscore <- nchar(string) + 1
                         string <- substr(x, 1, underscore - 1)
                         string
                       })))
}

comp_results[comp_results$algorithm == "defragTrees", "rule.length"] <- 
  sapply(comp_results[comp_results$algorithm == "defragTrees", "pretty.rule"]
         , relength_pretty_rule)

get_func_of <- function(func, metric, ...) {
  tapply(comp_results[[metric]]
         , list(comp_results$dataset
                , comp_results$algorithm)
         , func, ...)
}

get_mean_of <- function(metric) {
  get_func_of(mean, metric, na.rm = TRUE)
}
get_sd_of <- function(metric) {
  get_func_of(sd, metric, na.rm = TRUE)
}
get_median_of <- function(metric) {
  get_func_of(median, metric, na.rm = TRUE)
}

lwrq <- function(x) {
  quantile(x, probs = 0.25, na.rm = TRUE)
}
uprq <- function(x) {
  quantile(x, probs = 0.75, na.rm = TRUE)
}

get_lwrq_of <- function(metric) {
  get_func_of(lwrq, metric)
}
get_uprq_of <- function(metric) {
  get_func_of(uprq, metric)
}

above.75 <- function(x) {
  mean(x > 0.75) # shouldn't be any NA to rm.
}

get_above.75_of <- function(metric) {
  get_func_of(above.75, metric)
}

equal.zero <- function(x) {
  mean(x == 0) # shouldn't be any NA to rm.
}

get_equal.zero_of <- function(metric) {
  get_func_of(equal.zero, metric)
}

test_set_size_sqrt <- sapply(test_set_size, function(x) {
  sqrt(min(1000, x))
})

get_mean_ranks_of <- function(meas) {
  qmeas <- quo(meas)
  cres <- dplyr::select(comp_results, instance_id, dataset, !! qmeas, algorithm) %>%
    group_by(dataset) %>%
    spread(algorithm, !! qmeas)
  cres[, algorithms] <- t(apply(-cres[, algorithms], 1, rank))
  summarise(cres
            , Anchors = mean(Anchors)
            , BRL = mean(BRL)
            , `gbt-HIPS` = mean(`gbt-HIPS`)
            , defragTrees = mean(defragTrees)
            , inTrees = mean(inTrees)
            , LORE = mean(LORE))
}

post_hoc_ztest <- function(meas) {
  mr <- select(get_mean_ranks_of(meas), -dataset)
  
  top_rank <- as.vector(as.matrix(mr)[matrix(c(1:nrow(mr), apply(mr, 1, which.min)), ncol = 2)])
  scnd_rank <- c()
  for (i in 1:nrow(mr)) {
    scnd_rank[i] <- mr[i, ][, as.vector(mr[i, ] != top_rank[i])][which.min(mr[i, ][, as.vector(mr[i, ] != top_rank[i])])]
  }
  scnd_rank <- unlist(scnd_rank)
  ranks <- as.data.frame(t(matrix(c(top_rank, scnd_rank), ncol = 2)))
  names(ranks) <- datasetnames
  print(ranks)
  md <- scnd_rank - top_rank
  k <- c(6, 6, 6, 5, 6, 6, 4, 5, 6)
  z <- (md) / sqrt((k * (k + 1)) / (6 * test_set_size))
  print(z)
  ztest <- pnorm(z, lower.tail = FALSE)
  print(ztest)
  print("reject null bonferroni")
  print(ztest < 0.025/k)
}

get_fidelity <- function() {
  tapply(comp_results$predicted.class == comp_results$target.class
         , list(comp_results$dataset
                , comp_results$algorithm)
         , mean)
}

get_fidelity()

meas <- "stability.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)
get_above.75_of(meas)

meas <- "precision.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)

meas <- "coverage.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)

meas <- "wxcoverage.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)

meas <- "rule.length"
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_equal.zero_of("rule.length")

meas <- "elapsed_time"
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)

get_means_for_plotting <- function(meas, div = 1, select_algos = algorithms) {
  
  means <- as.data.frame(get_mean_of(meas) / div)[, select_algos]
  means$dataset <- rownames(means)
  means <- gather(means, algorithm, m, -dataset)
  
  st_errs <- as.data.frame((get_sd_of(meas) / test_set_size_sqrt) / div)[, select_algos]
  st_errs$dataset <- rownames(st_errs)
  st_errs <- gather(st_errs, algorithm, st_err, -dataset)
  
  means$lwr <- means$m - st_errs$st_err
  means$upr <- means$m + st_errs$st_err
  return(means)
}

get_points_for_plotting <- function(func_table, values_to_name, rows_to_name = "dataset") {
  vtn <- enquo(values_to_name)
  out <- as.data.frame(func_table)
  out[[rows_to_name]] <- dimnames(get_above.75_of(meas))[[1]]
  out <- gather(out, algorithm, !!vtn, -dataset)
  return(out)
}

# plot themes
source("KTheme.R")
reds <- k.grad.red.rev(2)
ongs <- k.grad.orange.rev(2)
blus <- k.grad.blue.rev(2)
grns <- k.grad.green.rev(2)
prps <- k.grad.purple.rev(2)
grys <- k.grad.grey.rev(2)

myPal1 <- c(reds[1], ongs[1], grns[1], blus[2], prps[1], grys[2])
myPal2 <- myPal1
myPal2 <- paste(myPal2, c(rep("55", 3), "", rep("55", 2)), sep = "")
myAlph1 <- c(rep(0.5, 3), 1, rep(0.5, 2))
myAlphFixed <- c(rep(0.25, 27), rep(1.0, 9), rep(0.25, 18))
myShap <- c(11, 2, 5, 15, 6, 1)
names(myPal1) <- algorithms
names(myPal2) <- algorithms
names(myAlph1) <- algorithms
names(myShap) <- algorithms
names(myAlphFixed) <- rep(algorithms, each = 9)
colos1 <- scale_colour_manual(
  values = myPal1)
colos2 <- scale_colour_manual(
  values = myPal2) # with alpha built in
alpas <- scale_alpha_manual(values = myAlph1)
shaps <- scale_shape_manual(values = myShap)
sens_colos_silent <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue)
                                         , guide = FALSE)
sens_colos <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue))


get_ylabel <- function(y) {
  ylabel <- gsub("cc", "coverage of target class"
                 , gsub("wx", "exclusive "
                   , gsub("stability", "reliability"
                          , gsub("_", " "
                              , gsub("rule.length", "Rule antecedent cardinality"
                                 , gsub(".tt.", "", y))))))
}

get_main_plot <- function(meas
                          , sc_y = scale_y_continuous()
                          , div = 1
                          , select_algos = algorithms) {
  myGeomPoint <- geom_point(size = 5)
  myGeomErrorBar <- geom_errorbar(width = 0.25)
  myGeomLine <- geom_line(size = 0.25, alpha = myAlphFixed[names(myAlphFixed) %in% select_algos])
  ggplot(data = get_means_for_plotting(meas, div, select_algos)
         , aes(y = m
               , ymin = lwr
               , ymax = upr
               , x = dataset
               , colour = algorithm
               , group = algorithm
               , alpha = algorithm
               , shape = algorithm)) +
    myGgTheme +
    myGeomLine +
    myGeomPoint +
    myGeomErrorBar +
    colos1 + 
    alpas +
    shaps +
    sc_y +
    ylab(get_ylabel(meas)) +
    xlab("data set") +
    theme(
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 26),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 26)
    )
}
get_main_plot("precision.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
get_main_plot("stability.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
get_main_plot("coverage.tt.", scale_y_continuous(limits = c(0.0, 1.00)))
get_main_plot("wxcoverage.tt.", scale_y_continuous(limits = c(0.0, 0.6)))
get_main_plot("rule.length")
get_main_plot("elapsed_time", scale_y_log10())

get_followup_plot <- function(
  summary_table
  , y_label
  #                         
  , sc_y = scale_y_continuous()
                          , select_algos = algorithms) {
  myGeomPoint <- geom_point(size = 5)
  myGeomLine <- geom_line(size = 0.25, alpha = myAlphFixed[names(myAlphFixed) %in% select_algos])
  ggplot(data = get_points_for_plotting(summary_table
                                        , values_to_name = "place_holder") # get_points_for_plotting(summary_table, values_to_name = values_to_name)
         , aes(y = place_holder
               , x = dataset
               , colour = algorithm
               , group = algorithm
               , alpha = algorithm
               , shape = algorithm)) +
    myGgTheme +
    myGeomLine +
    myGeomPoint +
    colos1 + 
    alpas +
    shaps +
    sc_y +
    ylab(y_label) +
    xlab("data set") +
    theme(
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 26),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 26)
    )
}
get_followup_plot(get_above.75_of("stability.tt.")
                  , y_label = "Proportion of reliability scores > 0.75")

algo_alphas <- myAlphFixed[
  as.vector(t(matrix(colnames(get_mean_of(meas))
       , byrow = TRUE
       , nrow = nrow(get_mean_of(meas))
       , ncol = ncol(get_mean_of(meas)))))[!is.na(as.vector(t(get_mean_of(meas))))]
  ]

get_main_boxplot <- function(meas
                             , sc_y = scale_y_continuous()) {
  ggplot(data = dplyr::select(comp_results, instance_id, dataset, !! enquo(meas), algorithm)
         , aes(y = !!enquo(meas)
               , colour = algorithm)) +
    geom_boxplot(outlier.alpha = 0.1
                 , alpha = algo_alphas
                 , position = position_dodge2(width = 2
                                              , padding = 0.2)
                 ) +
    colos2 +
    facet_wrap(dataset~.) +
    scale_x_continuous(labels = NULL
                       , breaks = 0) +
    sc_y +
    myGgTheme +
    theme(legend.position = "bottom"
          , legend.title = element_blank()
          , strip.text = element_text(size = 8
                                      , margin = ggplot2::margin(1,0,1,0, "pt"))
          ) +
    ylab(get_ylabel(as_label(enquo(meas)))) +
    theme(
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 26),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 26),
      strip.text = element_text(size = 30)
    )
}

get_main_boxplot(precision.tt.)
get_main_boxplot(stability.tt.)
get_main_boxplot(coverage.tt.)
get_main_boxplot(wxcoverage.tt.
                 , sc_y = scale_y_continuous(limits = c(0.0, 0.6)))
get_main_boxplot(rule.length)
get_main_boxplot(elapsed_time, scale_y_log10())


tikz(file = "cc.tikz", width = 6.85, height = 2)
get_main_plot("cc.tt."
              , sc_y = scale_y_continuous(limits = c(0.0, 1.0))
              , div = test_set_size
              , select_algos = c("Anchors", algo_variant))
dev.off()

tikz(file = "prec.tikz", width = 6.85, height = 2)
get_main_plot("precision.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

tikz(file = "stab.tikz", width = 6.85, height = 2)
get_main_plot("stability.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

tikz(file = "cov.tikz", width = 6.85, height = 2)
get_main_plot("coverage.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

tikz(file = "xcov.tikz", width = 6.85, height = 2)
get_main_plot("wxcoverage.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

tikz(file = "etime.tikz", width = 6.85, height = 2)
get_main_plot("elapsed_time", scale_y_log10())
dev.off()

tikz(file = "stabbox.tikz", width = 6.85, height = 3)
get_main_boxplot(stability.tt.)
dev.off()
tikz(file = "xcovbox.tikz", width = 6.85, height = 3)
get_main_boxplot(wxcoverage.tt.
                 , sc_y = scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

#########
get_lwrq_of(meas) / test_set_size
get_median_of(meas) / test_set_size
get_uprq_of(meas) / test_set_size

get_lwrq_of("stability.tt.")
get_median_of("stability.tt.")
get_uprq_of("stability.tt.")

get_meds_for_plotting <- function(meas, div = 1, select_algos = algorithms) {
  
  meds <- as.data.frame(get_median_of(meas)[, select_algos] / div)
  meds$dataset <- rownames(meds)
  meds <- gather(meds, algorithm, m, -dataset)
  lwrqs <- as.data.frame(get_lwrq_of(meas) / div)[, select_algos]
  lwrqs$dataset <- rownames(lwrqs)
  lwrqs <- gather(lwrqs, algorithm, lwr, -dataset)
  uprqs <- as.data.frame(get_uprq_of(meas) / div)[, select_algos]
  uprqs$dataset <- rownames(uprqs)
  uprqs <- gather(uprqs, algorithm, upr, -dataset)
  meds$lwr <- lwrqs$lwr
  meds$upr <- uprqs$upr
  return(meds)
  
}

get_med_plot <- function(meas
                         , sc_y = scale_y_continuous()
                         , div = 1
                         , select_algos = algorithms) {
  ylabel <- gsub("cc", "coverage of target class"
                 , gsub("wx", "exclusive "
                        , gsub("_", " "
                               , gsub(".tt.", "", meas))))
  myGeomPoint <- geom_point(size = 1.5)
  myGeomErrorBar <- geom_errorbar(width = 0.25)
  myGeomLine <- geom_line(size = 0.25, alpha = myAlphFixed[names(myAlphFixed) %in% select_algos])
  # y_label <- get_ylabel <- function(meas)
  ggplot(data = get_meds_for_plotting(meas, div, select_algos)
         , aes(y = m
               , ymin = lwr
               , ymax = upr
               , x = dataset
               , colour = algorithm
               , group = algorithm
               , alpha = algorithm
               , shape = algorithm)) +
    myGgTheme +
    myGeomLine +
    myGeomPoint +
    myGeomErrorBar +
    colos1 + 
    alpas +
    shaps +
    sc_y +
    # scale_x_discrete(labels = get_datasetname_stems(datasetnames)) +
    ylab(ylabel) +
    xlab("data set")
}
get_med_plot("stability.tt.", scale_y_continuous(limits = c(0.0, 1.0)))


get_med_plot("cc.tt."
           , sc_y = scale_y_continuous(limits = c(0.0, 1.0))
           , div = test_set_size
           , select_algos = c("Anchors", algo_variant))


meas <- "stability.tt."
for (ds in datasetnames) {
  qmeas <- quo(meas)
  cres <- comp_results %>% filter(dataset == ds) %>%
    dplyr::select(instance_id, dataset, !! qmeas, algorithm)
  cres <- spread(cres, algorithm, !! qmeas)
  mrs <- apply(apply(-cres[, algorithms], 1, rank), 1, mean)
  print(mrs)

  cres_frd <- friedman.test(as.matrix(cres))
  print(ds)
  print(cres_frd)
  
  chisqstat <- cres_frd$statistic
  df1 <- cres_frd$parameter
  N <- nrow(cres)
  fstat <- (chisqstat * (N - 1)) / (N * df1 - chisqstat)
  names(fstat) <- "Friedman F"
  df2 <- df1 * (N - 1)
  fpvalue <- pf(fstat, df1=df1, df2=df2, lower.tail = FALSE)
  print(fstat)
  print(fpvalue)
  
  cres_nem <- kwAllPairsNemenyiTest(cres[, algorithms])
  print(cres_nem)
  print(algorithms)
}

cres <- comp_results %>%
  dplyr::select(instance_id, dataset, !! qmeas, algorithm)
ggplot(data = cres # get_meds_for_plotting("stability.tt.")
       , aes(y = wxcoverage.tt.#m
             #, ymin = lwr
             #, ymax = upr
             , x = dataset
             , colour = algorithm
             #, group = algorithm
             , alpha = algorithm
             , shape = algorithm)) +
  myGgTheme +
  geom_boxplot(position = position_dodge()) +
  #geom_point() +
  #geom_errorbar(width = 0.2) +
  colos1 + 
  alpas

ggplot(data = comp_results, 
       aes(y = stability.tt., x = wxcoverage.tt.
           , colour = algorithm)) +
  geom_point(alpha = 0.1) +
  colos1 +
  myGgTheme_facets +
  labs(y = get_ylabel("stability.tt."), x = get_ylabel("wxcoverage.tt.")) +
  facet_grid(algorithm ~ dataset)

ggplot(data = comp_results, 
       aes(y = stability.tt., x = wxcoverage.tt.
           , colour = algorithm)) +
  geom_point(alpha = 0.1, size = 2) +
  colos1 +
  myGgTheme_facets +
  labs(y = get_ylabel("stability.tt."), x = get_ylabel("wxcoverage.tt.")) +
  scale_x_continuous(breaks = c(0.25, 0.75)) +
  scale_y_continuous(breaks = c(0.25, 0.75)) + 
  facet_grid(algorithm ~ dataset) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 30)
  )

ggplot(data = comp_results, 
       aes(y = rule.length, x = stability.tt.
           , colour = algorithm)) +
  geom_point(alpha = 0.1, size = 2) +
  colos1 +
  myGgTheme_facets +
  labs(x = get_ylabel("stability.tt."), y = get_ylabel("rule.length")) +
  scale_x_continuous(breaks = c(0.25, 0.75)) + 
  facet_grid(algorithm ~ dataset) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 30)
  )

ggplot(data = comp_results, 
       aes(y = rule.length, x = coverage.tt.
           , colour = algorithm)) +
  geom_point(alpha = 0.1, size = 2) +
  colos1 +
  myGgTheme_facets +
  labs(x = get_ylabel("wxcoverage.tt."), y = get_ylabel("rule.length")) +
  scale_x_continuous(breaks = c(0.25, 0.75)) + 
  facet_grid(algorithm ~ dataset) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 30)
  )
