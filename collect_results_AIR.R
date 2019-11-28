library(dplyr)
library(tidyr)
library(PMCMRplus)
library(cowplot)
options(max.print=20*72)
algorithms <- c("Anchors", "BRL", "CHIRPS", "defragTrees", "inTrees")

# setup
if (grep("linux", Sys.getenv()[['R_LIBS_USER']])) {
  pathsep <- "/"
} else {
  "\\"
}

# data management
source("data_files_mgmt.R")
datasets_master$difficulty <- c("large"
                              , "small"
                              , "small"
                              , "small"
                              , "small"
                              , "small"
                              , "small"
                              , "small"
                              ,"large")

sensdirs <- paste0(resfilesdirs, "rf_sensitivity", pathsep)

train_set_size <- integer(nrow(datasets_master))
test_set_size <- integer(nrow(datasets_master))
names(train_set_size) <- rownames(datasets_master)
names(test_set_size) <- rownames(datasets_master)

# sensitivity analysis
first_pass <- TRUE
for (i in seq_along(sensdirs)) {
  
  filepath <- normalizePath(file.path(sensdirs[i]))
  filenames <- dir(filepath)
  for (filename in filenames) {
    
    if (!grepl("summary", filename)) {
      supp <- gregexpr("sp\\_0.[0-9]{1,2}", filename)
      supp <- regmatches(filename, supp)[[1]]
      supp <- as.numeric(gsub("sp\\_", "", supp))
      alpp <- gregexpr("ap\\_0.[1-9]", filename)
      alpp <- regmatches(filename, alpp)[[1]]
      alpp <- as.numeric(gsub("ap\\_", "", alpp))
      dpb <- gregexpr("dpb\\_[1-9]{1,2}", filename)
      dpb <- regmatches(filename, dpb)[[1]]
      dpb <- as.numeric(gsub("dpb\\_", "", dpb))
      sf <- gregexpr("sf\\_[1-9]", filename)
      sf <- regmatches(filename, sf)[[1]]
      sf <- as.numeric(gsub("sf\\_", "", sf))
      wht <- gregexpr("w\\_(kldiv|nothing|lodds)", filename)
      wht <- regmatches(filename, wht)[[1]]
      wht <- gsub("w\\_", "", wht)
      
      results <- read.csv(paste0(sensdirs[i], filename), stringsAsFactors = FALSE)
      
      # datasetname <- rep(dsname, nrow(results))
      support <- rep(supp, nrow(results))
      alpha_paths <- rep(alpp, nrow(results))
      disc_path_bins <- rep(dpb, nrow(results))
      score_func <- rep(sf, nrow(results))
      weighting <- rep(wht, nrow(results))
      
      results <- cbind(results, as.data.frame(cbind(support, alpha_paths, disc_path_bins, score_func, weighting), stringsAsFactors = FALSE))
      if (first_pass) {
        main_results <- results
        first_pass <- FALSE
      } else {
        main_results <- rbind(main_results, results)
      }
    } else { # summary

    }
  }
}

sens_results <- main_results %>% mutate(
  dataset_name = factor(dataset_name)
  , true.class = factor(true.class)
  , true.class.label = factor(true.class.label)
  , predicted.class = factor(predicted.class)
  , predicted.class.label = factor(predicted.class.label)
  , target.class = factor(target.class)
  , target.class.label = factor(target.class.label)
  , support = factor(ifelse(ifelse(datasets_master[dataset_name, "difficulty"] == "large", as.numeric(support) - 0.1, as.numeric(support)) == 0.1, "A", "B"))
  , weighting = factor(weighting)
  , score_func = factor(score_func)
  , wxcoverage.tr. = xcoverage.tr. * (nci.tr./(ci.tr. + nci.tr.))
  , wxcoverage.tt. = xcoverage.tt. * (nci.tt./(ci.tt. + nci.tt.))
)

# get the listing of all the sensitivity analyses by grid
sens_groups <- with(sens_results, expand.grid(
  support = unique(support)
  , alpha = unique(alpha_paths)
  , bins = unique(disc_path_bins)
  , func = unique(score_func)
  , weights = unique(weighting)))
# provide id numbers
sens_groups$id <- as.numeric(rownames(sens_groups))

get_sensitivity <- function(measure) {
  sensitivity_analysis <- list()
  for (ds in rownames(datasets_master)) {
    sv <- filter(sens_results, dataset_name == ds)
    # order must be same as sens_groups above
    sens_values <- tapply(sv[[measure]], list(sv$support
                                              , sv$alpha_paths
                                              , sv$disc_path_bins
                                              , sv$score_func
                                              , sv$weighting), identity)
    
    sens_values <- matrix(unlist(sens_values), ncol = nrow(sens_groups))
    test_set_size[ds] <- nrow(sens_values)
    # get a first taste of results
    sens_groups$sd <- apply(sens_values, 2, sd)
    sens_groups$st.err <- sens_groups$sd / sqrt(test_set_size[ds])
    mn <- apply(sens_values, 2, mean)
    sens_groups$lwr_mn_ci <- mn - sens_groups$st.err
    sens_groups$mn <- mn
    sens_groups$upr_mn_ci <- mn + sens_groups$st.err
    qntl <- apply(sens_values, 2, quantile)
    sens_groups$lwrq <- qntl[2, ]
    sens_groups$med <- apply(sens_values, 2, median)
    sens_groups$uprq <- qntl[4, ]
    sens_groups$rank_mean <- colMeans(t(apply(sens_values, 1, rank)))
    sens_groups$rank_sum <- colSums(t(apply(sens_values, 1, rank)))
    
    sensitivity_analysis[[ds]] <- sens_groups
  }
  return(sensitivity_analysis)
}

get_best_sens <- function(sens) {
  for (ds in rownames(datasets_master)) {
    print(ds)
    print(sens[[ds]][which.max(sens[[ds]][["rank_mean"]]), ])
  }
}

stability_tr <- get_sensitivity("stability.tr.")
get_best_sens(stability_tr)

stability_sens <- get_sensitivity("stability.tt.")
wxcoverage_sens <- get_sensitivity("wxcoverage.tt.")
rule.length_sens <- get_sensitivity("rule.length")

# comparative analysis
patt <- paste0("(", paste(algorithms, collapse = ")|("), ")")

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
    if (grepl(patt, fn)) {
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
        # if (!any(grepl("sd_proxy_performance", names(results)))) { # no med rc column
        #   results <- results %>%
        #     mutate(median_rule_cascade = NA)
        # }
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
comp_results$algorithm <- ifelse(comp_results$algorithm == "greedy_stab", "CHIRPS", comp_results$algorithm)

comp_results$pretty.rule[is.na(comp_results$pretty.rule)] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == ""] <- "{default}"
comp_results$pretty.rule[comp_results$pretty.rule == "[]"] <- "{default}"

# adjust 0-length rules to be valued as 1, {default}
comp_results$rule.length <- ifelse(comp_results$pretty.rule == "{default}"
                                   , 0
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
            , CHIRPS = mean(CHIRPS)
            , defragTrees = mean(defragTrees)
            , inTrees = mean(inTrees))  
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
  k <- rep(5, 9)
  k[7] <- 3
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

meas <- "cc.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)

meas <- "elapsed_time"
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)

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

# plot themes
source("KTheme.R")
reds <- k.grad.red.rev(2)
ongs <- k.grad.orange.rev(2)
blus <- k.grad.blue.rev(2)
grns <- k.grad.green.rev(2)
prps <- k.grad.purple.rev(2)

myPal <- c(reds[2], ongs[2], blus[2], grns[2], k.pink)
myAlph <- c(0.3, 0.3, 1, 0.3, 0.3)
myShap <- c(1, 2, 15, 5, 6)
names(myPal) <- algorithms
names(myAlph) <- algorithms
names(myShap) <- algorithms
colos <- scale_colour_manual(
  values = myPal)
alpas <- scale_alpha_manual(values = myAlph)
shaps <- scale_shape_manual(values = myShap)
sens_colos_silent <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue)
                                         , guide = FALSE)
sens_colos <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue))

get_main_plot <- function(meas
                          , sc_y = scale_y_continuous()
                          , div = 1
                          , select_algos = algorithms) {
  ylabel <- gsub("cc", "coverage of target class"
                 , gsub("wx", "exclusive "
                        , gsub("_", " "
                               , gsub(".tt.", "", meas))))
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
    geom_line() +
    geom_point(size = 0.8) +
    geom_errorbar(width = 0.2) +
    colos + 
    alpas +
    shaps +
    sc_y +
    scale_x_discrete(labels = get_datasetname_stems(datasetnames)) +
    ylab(ylabel) +
    xlab("data set")
}

get_med_plot <- function(meas
                         , sc_y = scale_y_continuous()
                         , div = 1
                         , select_algos = algorithms) {
  ylabel <- gsub("cc", "coverage of target class"
                 , gsub("wx", "exclusive "
                        , gsub("_", " "
                               , gsub(".tt.", "", meas))))
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
    geom_line() +
    geom_point(size = 0.8) +
    geom_errorbar(width = 0.2) +
    colos + 
    alpas +
    shaps +
    sc_y +
    scale_x_discrete(labels = get_datasetname_stems(datasetnames)) +
    ylab(ylabel) +
    xlab("data set")
}

get_med_plot("cc.tt."
             , sc_y = scale_y_continuous(limits = c(0.0, 1.0))
             , div = test_set_size
             , select_algos = c("Anchors", "CHIRPS"))

tikz(file = "cc.tikz", width = 6.85, height = 2)
get_main_plot("cc.tt."
              , sc_y = scale_y_continuous(limits = c(0.0, 1.0))
              , div = test_set_size
              , select_algos = c("Anchors", "CHIRPS"))
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

# sensitivity plots
best_sens <- lapply(rownames(datasets_master), function(ds) {
  sens <- get_sensitivity("stability.tr.")
  sens[[ds]][which.max(sens[[ds]][["rank_mean"]]), ][["id"]]
})
names(best_sens) <- rownames(datasets_master)
idx <- c(1:6, 13:18, 25:30, 7:12, 19:24, 31:36)
idx <- c(idx, idx + 36)

for (ds in rownames(datasets_master)) {
  sens_plotting <- stability_sens[[ds]]
  sens_plotting$best <- sens_plotting$id == best_sens[[ds]]
  sens_plotting <- sens_plotting[idx, ] # re-order
  sens_plotting$id_plotting <- c(1:18, 20:37, 39:56, 58:75)
  tr_top_mean <- sens_plotting$id_plotting[which(sens_plotting$best)]
  arrowbase <- min(sens_plotting$lwr_mn_ci)
  arrowtip <- arrowbase + (max(sens_plotting$upr_mn_ci) - arrowbase) / 8
  g <- ggplot(data = sens_plotting
         , aes(y = mn, ymin = lwr_mn_ci, ymax = upr_mn_ci
               , x = id_plotting
               , shape = alpha
               , colour = func
               , size = support)) +
    geom_point() +
    geom_errorbar(size=0.5, width=1) +
    geom_segment(x = tr_top_mean
                 , xend = tr_top_mean
                 , y = arrowbase
                 , yend = arrowtip
                 , arrow = arrow(length = unit(0.025, "npc"))
                 , colour = myPalNeut[7]
                 , size = 0.5) +
    labs(title = get_datasetname_stems(ds)
         , x = NULL, y = NULL) +
    scale_size_discrete(range = c(0.8, 1.2), guide = FALSE) +
    scale_shape(guide = FALSE) +
    sens_colos_silent +
    geom_vline(xintercept = c(19, 57)
               , linetype = 2
               , colour = myPalNeut[4]) +
    geom_vline(xintercept = 38
               , linetype = 1
               , colour = myPalNeut[3]) +
    myGgTheme +
    theme(axis.text.x = element_blank()
          , axis.ticks.x = element_blank())
  
  tikz(file = paste0(get_datasetname_stems(ds), "_sens.tikz"), width = 2.2, height = 1.5)
  print(g)
  dev.off()
}
# plot to get only the legend
legend_plot <- ggplot(data = sens_plotting
             , aes(y = mn, ymin = lwr_mn_ci, ymax = upr_mn_ci
                   , x = id_plotting
                   , shape = alpha
                   , colour = func
                   , size = support)) +
  geom_point() +
  geom_errorbar(size=0.5, width=1) +
  scale_size_discrete(range = c(0.8, 1.2), labels = c("0.1", "0.2")) +
  sens_colos +
  theme(legend.box = "horizontal"
        , legend.key = element_blank())

lgnd <- get_legend(legend_plot)
grid.newpage()
tikz(file = "legend_sens.tikz", width = 2.2, height = 1.5)
grid.draw(lgnd)
dev.off()


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
  colos + 
  alpas
