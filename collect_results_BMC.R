library(dplyr)
library(tidyr)
library(ggplot2)
library(tikzDevice)
library(PMCMRplus)
library(rlang)
source("data_files_mgmt_bmc.R")

BMCdir <- "BMC" # "BMC2"
BMCsensdir <- "ada1_sensitivity" # "ada2_sensitivity"

resfilesdirs <- paste0(project_dir, BMCdir, pathsep, datasetnames, pathsep)
sensdirs <- paste0(resfilesdirs, BMCsensdir, pathsep)

# selector <- 5:6

# sensitivity analysis
# sensdirs <- sensdirs[selector] # temporary
first_pass <- TRUE
for (i in seq_along(sensdirs)) {
  
  filepath <- normalizePath(file.path(sensdirs[i]))
  filenames <- dir(filepath)
  for (filename in filenames) {
    
    if (!grepl("summary", filename)) {
      wcts <- gregexpr("wcts\\_(conf\\_weighted|majority)", filename)
      wcts <- regmatches(filename, wcts)[[1]]
      wcts <- substr(gsub("wcts\\_", "", wcts), 1, 3)
      supp <- gregexpr("sp\\_0.[0-9]{1,3}", filename)
      supp <- regmatches(filename, supp)[[1]]
      supp <- as.numeric(gsub("sp\\_", "", supp))
      alpp <- gregexpr("ap\\_0.[0-9]", filename)
      alpp <- regmatches(filename, alpp)[[1]]
      alpp <- as.numeric(gsub("ap\\_", "", alpp))
      dpb <- gregexpr("dpb\\_[1-9]{1,2}", filename)
      dpb <- regmatches(filename, dpb)[[1]]
      dpb <- as.numeric(gsub("dpb\\_", "", dpb))
      dpeq <- gregexpr("dpeq\\_(True|False)", filename)
      dpeq <- regmatches(filename, dpeq)[[1]]
      dpeq <- as.logical(toupper(gsub("dpeq\\_", "", dpeq)))
      sf <- gregexpr("sf\\_[1-9]", filename)
      sf <- regmatches(filename, sf)[[1]]
      sf <- as.numeric(gsub("sf\\_", "", sf))
      wht <- gregexpr("w\\_(kldiv|nothing|chisq)", filename)
      wht <- regmatches(filename, wht)[[1]]
      wht <- gsub("w\\_", "", wht)
      
      results <- read.csv(paste0(sensdirs[i], filename), stringsAsFactors = FALSE)
      
      # datasetname <- rep(dsname, nrow(results))
      wcts <- rep(wcts, nrow(results))
      support <- rep(supp, nrow(results))
      alpha_paths <- rep(alpp, nrow(results))
      disc_path_bins <- rep(dpb, nrow(results))
      disc_path_eqcounts <- rep(dpeq, nrow(results))
      score_func <- rep(sf, nrow(results))
      weighting <- rep(wht, nrow(results))
      
      results <- cbind(results, as.data.frame(cbind(wcts, support, alpha_paths, disc_path_bins, disc_path_eqcounts, score_func, weighting), stringsAsFactors = FALSE))
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

# datasets_master <- datasets_master[selector, ]

sens_results <- main_results %>% mutate(
  dataset_name = factor(dataset_name)
  , true.class = factor(true.class)
  , true.class.label = factor(true.class.label)
  , predicted.class = factor(predicted.class)
  , predicted.class.label = factor(predicted.class.label)
  , target.class = factor(target.class)
  , target.class.label = factor(target.class.label)
  , wcts = factor(wcts)
  , support = factor(support)
  , weighting = factor(weighting)
  , score_func = factor(score_func)
  , wxcoverage.tr. = xcoverage.tr. * (nci.tr./(ci.tr. + nci.tr.))
  , wxcoverage.tt. = xcoverage.tt. * (nci.tt./(ci.tt. + nci.tt.))
)

# get the listing of all the sensitivity analyses by grid
sens_groups <- with(sens_results, expand.grid(
  wcts = unique(wcts)
  , support = unique(support)
  , alpha = unique(alpha_paths)
  , bins = unique(disc_path_bins)
  , eqcounts = unique(disc_path_eqcounts)
  , func = unique(score_func)
  , weights = unique(weighting)))
# provide id numbers
sens_groups$id <- as.numeric(rownames(sens_groups))

get_sensitivity <- function(measure) {
  sensitivity_analysis <- list()
  for (ds in rownames(datasets_master)) {
    sv <- filter(sens_results, dataset_name == ds)
    # order must be same as sens_groups above
    sens_values <- tapply(sv[[measure]], list(sv$wcts
                                              , sv$support
                                              , sv$disc_path_bins
                                              , sv$disc_path_eqcounts
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

# comparative analysis
source("data_files_mgmt_bmc.R")
BMCdir <- "BMC" # "BMC2"

resfilesdirs <- paste0(project_dir, BMCdir, pathsep, datasetnames, pathsep)

algorithms <- c("Anchors", "CHIRPS", "lore")
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
            , CHIRPS = mean(CHIRPS)
            , lore = mean(lore))  
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

get_zeros_count <- function(meas) {
  tapply(comp_results[[meas]] == 0.0
         , list(comp_results$dataset
                , comp_results$algorithm)
         , mean)
}

get_fidelity()
get_zeros_count("precision.tt.")

get_wilcox_tests <- function(meas) {
  qmeas <- quo(meas)
  cres <- dplyr::select(comp_results, instance_id, dataset, !! qmeas, algorithm) %>%
    group_by(dataset) %>%
    spread(algorithm, !! qmeas)

  for (ds in datasetnames[c(5:9)]) {
    wcox <- cres %>% filter(dataset == ds) %>% ungroup() %>%
      dplyr::select(-instance_id, -dataset)
    wxtest <- wilcox.test(wcox[[algorithms[1]]], wcox[[algorithms[2]]], paired = TRUE)
    print(ds)
    print(wxtest)
  }
}

get_sign_tests <- function(meas) {
  qmeas <- quo(meas)
  cres <- dplyr::select(comp_results, instance_id, dataset, !! qmeas, algorithm) %>%
    group_by(dataset) %>%
    spread(algorithm, !! qmeas)
  
  for (ds in datasetnames[c(5:9)]) {
    wcox <- cres %>% filter(dataset == ds) %>% ungroup() %>%
      dplyr::select(-instance_id, -dataset)
    sgn <- sign(wcox[[algorithms[1]]] - wcox[[algorithms[2]]])
    s <- c(sum(sgn == 1, na.rm = TRUE), sum(sgn == -1, na.rm = TRUE))
    t <- sum(sgn != 0, na.rm = TRUE)
    wxtest <- binom.test(s, t)
    print(ds)
    print(wxtest)
  }
}

get_t_tests <- function(meas) {
  qmeas <- quo(meas)
  cres <- dplyr::select(comp_results, instance_id, dataset, !! qmeas, algorithm) %>%
    group_by(dataset) %>%
    spread(algorithm, !! qmeas)
  
  for (ds in datasetnames[c(5, 7:9)]) {
    wcox <- cres %>% filter(dataset == ds) %>% ungroup() %>%
      dplyr::select(-instance_id, -dataset)
    wxtest <- t.test(wcox[[algorithms[1]]], wcox[[algorithms[2]]], paired = TRUE)
    print(ds)
    print(wxtest)
  }
}
get_friedman_tests <- function(meas) {
  qmeas <- quo(meas)
  cres <- dplyr::select(comp_results, instance_id, dataset, !! qmeas, algorithm) %>%
    group_by(dataset) %>%
    spread(algorithm, !! qmeas)
  for (ds in datasetnames[c(1:4, 6)]) {
    fried <- cres %>% filter(dataset == ds) %>% ungroup() %>%
      dplyr::select(-instance_id, -dataset)
    fried_test <- friedman.test(as.matrix(fried))
    print(ds)
    print(fried_test)
    chisqstat <- fried_test$statistic
    df1 <- fried_test$parameter
    N <- nrow(fried)
    fstat <- (chisqstat * (N - 1)) / (N * df1 - chisqstat)
    names(fstat) <- "Friedman F"
    df2 <- df1 * (N - 1)
    fpvalue <- pf(fstat, df1=df1, df2=df2, lower.tail = FALSE)
    print(fstat)
    print(fpvalue)
    fried_nem <- kwAllPairsNemenyiTest(fried[, algorithms])
    print(fried_nem)
    print(0.025/(df1+1))
  }
}

get_means_for_plotting <- function(meas, div = 1, select_algos = algorithms) {
  
  means <- as.data.frame(get_mean_of(meas) / div)[, select_algos]
  names(means) <- ifelse(names(means) == "lore", "lore", names(means))
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

algo2 <- c("Anchors", "Ada-WHIPS", "LORE")
myPal1 <- c(reds[1], blus[2], ongs[2])
myPal2 <- myPal1
myPal2 <- paste(myPal2, c("99", "", "99"), sep = "")
myAlph1 <- c(0.6, 1, 0.6)
myAlphFixed <- c(rep(0.25, 9), rep(0.75, 9), rep(0.25, 9))
myShap <- c(1, 15, 2)
names(myPal1) <- algorithms
names(myPal2) <- algorithms
names(myAlph1) <- algorithms
names(myShap) <- algorithms
names(myAlphFixed) <- rep(algorithms, each = 9)
colos1 <- scale_colour_manual(
  values = myPal1, labels = algo2)
colos2 <- scale_colour_manual(
  values = myPal2, labels = algo2) # with alpha built in
alpas <- scale_alpha_manual(values = myAlph1, labels = algo2)
shaps <- scale_shape_manual(values = myShap, labels = algo2)
sens_colos_silent <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue)
                                         , guide = FALSE)
sens_colos <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue))

get_ylabel <- function(y) {
  ylabel <- gsub("cc", "coverage of target class"
                 , gsub("wx", "exclusive "
                        , gsub("_", " "
                               , gsub("\\.tt\\.", ""
                                      , gsub("e\\.l", "e l", y)))))
}

datasetlabels <- c("breast cancer"
                   , "cardiotocography"
                   , "diabetic retinopathy"
                   , "cleveland heart"
                   , "mental health survey 14"
                   , "mental health survey 16"
                   , "hospital readmission"
                   , "thyroid"
                   , "understanding society")
names(datasetlabels) <- datasetnames                      

dataset_labeller <- function(variable,value){
  return(datasetlabels[value])
}

get_main_plot <- function(meas
                          , sc_y = scale_y_continuous()
                          , div = 1
                          , select_algos = algorithms) {
  myGeomPoint <- geom_point(size = 1.5)
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
    scale_x_discrete(labels = datasetlabels) +
    ylab(get_ylabel(meas)) +
    xlab("data set") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

get_main_boxplot <- function(meas
                             , sc_y = scale_y_continuous()
                             , oa = 0.1) {
  
  ggplot(data = dplyr::select(comp_results, instance_id, dataset, !! enquo(meas), algorithm)
         , aes(y = !!enquo(meas)
               , colour = algorithm)) +
    geom_boxplot(outlier.alpha = oa
                 , position = position_dodge2(width = 2
                                              , padding = 0.2)
    ) +
    colos2 +
    facet_wrap(dataset~., labeller = dataset_labeller) +
    scale_x_continuous(labels = NULL
                       , breaks = NULL) +
    sc_y  +
    ylab(get_ylabel(as_label(enquo(meas)))) +
    myGgTheme +
    theme(legend.position = "bottom"
          , legend.title = element_blank()
          , strip.text = element_text(size = 8
                                      , margin = margin(1,0,1,0, "pt"))
    )
  
}

meas <- "coverage.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)[, c(2, 1, 3)]
round(st_err, 4)[, c(2, 1, 3)]
get_mean_ranks_of(meas)[, c(1, 3, 2, 4)]
post_hoc_ztest(meas)
get_wilcox_tests(meas)
get_sign_tests(meas)
# tikz(file = paste0(BMCdir, pathsep, "meancov.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "meancov.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
png(file = paste0(BMCdir, pathsep, "meancov.png")
    , width = 6.85, height = 3, units = "in", res = 1800)
get_main_plot("coverage.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()
# tikz(file = paste0(BMCdir, pathsep, "covbox.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "covbox.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
get_main_boxplot(coverage.tt.
                 , oa = 0.3
                 , sc_y = scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

meas <- "precision.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)[, c(2, 1, 3)]
round(st_err, 4)[, c(2, 1, 3)]
get_mean_ranks_of(meas)[, c(1, 3, 2, 4)]
get_t_tests(meas)
get_wilcox_tests(meas)
get_sign_tests(meas)
get_friedman_tests(meas)
post_hoc_ztest(meas)
get_zeros_count(meas)[, c(2,1,3)]
# tikz(file = paste0(BMCdir, pathsep, "meanprec.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "meanprec.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
get_main_plot("precision.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()
# tikz(file = paste0(BMCdir, pathsep, "precbox.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "precbox.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
get_main_boxplot(precision.tt., oa = 0.3)
dev.off()

meas <- "stability.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)[, c(2, 1, 3)]
round(st_err, 4)[, c(2, 1, 3)]
get_mean_ranks_of(meas)[, c(1, 3, 2, 4)]
post_hoc_ztest(meas)
get_wilcox_tests(meas)
get_t_tests(meas)
get_friedman_tests(meas)
# tikz(file = paste0(BMCdir, pathsep, "meanstab.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "meanstab.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
get_main_plot("stability.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()
# tikz(file = paste0(BMCdir, pathsep, "stabbox.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "stabbox.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
get_main_boxplot(stability.tt.)
dev.off()

meas <- "wxcoverage.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)
get_main_plot("wxcoverage.tt.", scale_y_continuous(limits = c(0.0, 1.0)))


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
# tikz(file = paste0(BMCdir, pathsep, "meantime.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "meantime.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
get_main_plot("elapsed_time", scale_y_continuous(trans = "log10"))
dev.off()
# tikz(file = paste0(BMCdir, pathsep, "timebox.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "timebox.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
get_main_boxplot(elapsed_time, oa = 0.5
                 , sc_y = scale_y_continuous(trans = "log10"))

dev.off()

meas <- "rule.length"
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)
post_hoc_ztest(meas)
# tikz(file = paste0(BMCdir, pathsep, "rllnbox.tikz"), width = 6.85, height = 3)
pdf(file = paste0(BMCdir, pathsep, "rllnbox.pdf")
    , width = 6.85, height = 3
    , compress = FALSE
    , title = NULL
    , colormodel = "cmyk")
get_main_boxplot(elapsed_time, oa = 0.5
                 , sc_y = scale_y_log10())
dev.off()


# more examples
get_main_boxplot(rule.length, oa = 1)
get_main_boxplot(coverage.tt.
                 , sc_y = scale_y_continuous(limits = c(0.0, 1.0)))



