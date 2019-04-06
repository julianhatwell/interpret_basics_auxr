library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(tikzDevice)
library(cowplot) # get legend out of a ggplot

# plot themes
source("C:\\Users\\id126493\\OneDrive\\Documents\\PhD\\KTheme.R")
reds <- k.grad.red.rev(2)
ongs <- k.grad.orange.rev(2)
blus <- k.grad.blue.rev(2)
grns <- k.grad.green.rev(2)
prps <- k.grad.purple.rev(2)
nuts <- myPalNeut
algs <- c(algorithms, "CHIRPS")
myPal <- c(reds[2], ongs[2], blus[2], grns[2], k.pink)
myAlph <- c(rep(0.4, length(algs) - 1), 1)
names(myPal) <- algs
names(myAlph) <- algs
colos <- scale_colour_manual(
  values = myPal)
alpas <- scale_alpha_manual(values = myAlph)
sens_colos_silent <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue)
                                          , guide = FALSE)
sens_colos <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue))
# fudge until added absolute coverage, train and test set sizes to results
train_set_size <- integer(length(datasetnames))
test_set_size <- integer(length(datasetnames))

for (i in seq_along(datasetnames)) {
  train_set_size[i] <- get_train_test_sizes(i)[1]
  test_set_size[i] <- get_train_test_sizes(i)[2]
}

names(train_set_size) <- datasetnames
names(test_set_size) <- datasetnames

calc_wxcov <- function(results_set, datasetname, sensitivity = TRUE) {
  # filter all results to one dataset
  analysis_out <- results_set %>% 
    filter(dataset == datasetname)
  
  if (sensitivity == TRUE) {
    analysis_out <- analysis_out %>% mutate(support = factor(support)
                            , alpha_paths = factor(alpha_paths)
                            , disc_path_bins = factor(disc_path_bins)
                            , score_func = factor(score_func)
                            , weighting = factor(weighting)
                            )
  }
  
  # analysis_out <- analysis_out %>% mutate(
  #          # covered is train set size * coverage
  #          covered.tr. = (train_set_size[datasetname] * coverage.tr.)
  #          # covered correct is (covered + 1 (explanandum)) * stability
  #          # , cc.tr. = (covered.tr. + 1) * stability.tr.
  #          , cc.tr. = covered.tr. * precision.tr.
  #          # covered incorrect is covered - cc
  #          , ci.tr. = covered.tr. - cc.tr.
  #          , stability.tr.lapl. = (cc.tr. + 1) / (covered.tr. + n_classes[datasetname] + 1) # laplace and stability
  #          , xcoverage.tr.lapl. = (covered.tr. + 1) / (train_set_size[datasetname] + n_classes[datasetname] + 1) # laplace and xcoverage
  #          , odds.tr.lapl. = (cc.tr. + 1 + 1/n_classes[datasetname])/(ci.tr. + (n_classes[datasetname] - 1) + (n_classes[datasetname] - 1)/n_classes[datasetname])
  #          , lodds.tr.lapl. = log(odds.tr.lapl.)
  #          # covered is test set size - 1  * coverage (because of LOO)
  #          , covered.tt. = ((test_set_size[datasetname] - 1) * coverage.tt.)
  #          # covered correct is (covered + 1 (explanandum)) * stability
  #          # , cc.tt. = (covered.tt. + 1) * stability.tt.
  #          , cc.tt. = covered.tt. * precision.tt.
  #          # covered incorrect is covered - cc
  #          , ci.tt. = covered.tt. - cc.tt.
  #          , stability.tt.lapl. = (cc.tt. + 1) / (covered.tt. + n_classes[datasetname] + 1) # laplace and stability
  #          , xcoverage.tt.lapl. = ((covered.tt. + 1) / (test_set_size[datasetname] + n_classes[datasetname] + 1)) * # laplace and xcoverage
  #            ( ci.tt. / (ci.tt. + (test_set_size[datasetname] - covered.tt.)))
  #          , odds.tt.lapl. = (cc.tt. + 1 + 1/n_classes[datasetname])/(ci.tt. + (n_classes[datasetname] - 1) + (n_classes[datasetname] - 1)/n_classes[datasetname])
  #          , lodds.tt.lapl. = log(odds.tt.lapl.))
  analysis_out <- analysis_out %>% mutate(tnr.tr. = nci.tr./(ci.tr. + nci.tr.)
    , tnr.tt. = nci.tt./(ci.tt. + nci.tt.)
    , wxcoverage.tr. = tnr.tr. * xcoverage.tr.
    , wxcoverage.tt. = tnr.tt. * xcoverage.tt.
    , npv.tr. = nci.tr./(ncc.tr. + nci.tr.)
    , npv.tt. = nci.tt./(ncc.tt. + nci.tt.)
    , wxcoverage2.tr. = npv.tr. * xcoverage.tr.
    , wxcoverage2.tt. = npv.tt. * xcoverage.tt.)
  return(analysis_out)
}

get_CHIRPS_analysis <- function(measure, results_set
                                , top_mean_block = NA
                                , top_ranksum_block = NA) {
  # sensitivity analysis
  analysis_out <- list()

  for (ds in datasetnames) {
  
    # results collection
    analysis_out[[ds]] <- list()
    
    # filter all results to one dataset
    analysis <- calc_wxcov(results_set, ds)
    
    # get the listing of all the sensitivity analyses by grid
    analysis_groups <- with(analysis, expand.grid(
      support = unique(support)
      , alpha = unique(alpha_paths)
      , bins = unique(disc_path_bins)
      , func = unique(score_func)
      , weights = unique(weighting)))
    
    # provide id numbers
    analysis_groups$id <- as.numeric(rownames(analysis_groups))
    
    # extract the values of "measure"
    analysis_values <- tapply(analysis[, measure]
                              , list(analysis$support
                                     , analysis$alpha_paths
                                     , analysis$disc_path_bins
                                     , analysis$score_func
                                     , analysis$weighting)
                              , identity)
  
    analysis_values <- matrix(unlist(analysis_values), ncol = nrow(analysis_groups))
    analysis_out[[ds]][["analysis_values"]] <- analysis_values
    
    # get a first taste of results
    analysis_groups$mean <- apply(analysis_values, 2, mean)
    analysis_groups$sd <- apply(analysis_values, 2, sd)
    analysis_groups$st.err <- analysis_groups$sd / sqrt(test_set_size[ds])
    analysis_groups$rank_mean <- colMeans(t(apply(analysis_values, 1, rank)))
    analysis_groups$rank_sum <- colSums(t(apply(analysis_values, 1, rank)))
    
    # collect groups results
    analysis_out[[ds]][["analysis_groups"]] <- analysis_groups
    
    # which is the best mean measure discover or feed in from best stability
    if (is.na(top_mean_block)) {
      ds_top_mean_block <- which.max(analysis_groups$mean)
    } else {
      # extract the values of "top_mean_block"
      tr_values <- tapply(analysis[, top_mean_block]
                                , list(analysis$support
                                       , analysis$alpha_paths
                                       , analysis$disc_path_bins
                                       , analysis$score_func
                                       , analysis$weighting)
                                , identity)
      tr_values <- matrix(unlist(tr_values), ncol = nrow(analysis_groups))
      
      # get a first taste of results
      tr_rank_mean <- colMeans(t(apply(tr_values, 1, rank)))
      ds_top_mean_block <- which.max(tr_rank_mean)
    }
    analysis_out[[ds]][["top_mean_block"]] <- ds_top_mean_block
    
    if (is.na(top_ranksum_block)) {
      ds_top_ranksum_block <- which.max(analysis_groups$rank_sum)
    } else {
      # extract the values of "top_mean_block"
      tr_values <- tapply(analysis[, top_ranksum_block]
                          , list(analysis$support
                                 , analysis$alpha_paths
                                 , analysis$disc_path_bins
                                 , analysis$score_func
                                 , analysis$weighting)
                          , identity)
      tr_values <- matrix(unlist(tr_values), ncol = nrow(analysis_groups))
      
      # get a first taste of results
      tr_rank_sum <- colSums(t(apply(tr_values, 1, rank)))
      ds_top_ranksum_block <- which.max(tr_rank_sum)
    }
    analysis_out[[ds]][["top_ranksum_block"]] <- ds_top_ranksum_block
    
    top_group_stats <- filter(analysis_groups, id == ds_top_ranksum_block) %>%
      dplyr::select(support, alpha, bins, func, weights)
    
    analysis_out[[ds]][["ds_bestmean_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out[[ds]][["top_mean_block"]]]
    analysis_out[[ds]][["ds_bestranksum_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out[[ds]][["top_ranksum_block"]]]
    
    sens_raw <- analysis %>% filter(support == top_group_stats$support
                                    , alpha_paths == top_group_stats$alpha
                                    , disc_path_bins == top_group_stats$bins
                                    , score_func == top_group_stats$func
                                    , weighting == top_group_stats$weights)
    
    analysis_out[[ds]][["sens_raw"]] <- sens_raw
  }
  mean_all_ds <- matrix(NA
                        , ncol = length(datasetnames)
                        , nrow = ncol(analysis_out[[1]]$analysis_values))
  ranksum_all_ds <- matrix(NA
                           , ncol = length(datasetnames)
                           , nrow = ncol(analysis_out[[1]]$analysis_values))
  
  for (j in seq_along(datasetnames)) {
    mean_all_ds[, j] <- analysis_out[[datasetnames[j]]]$analysis_groups$mean
    ranksum_all_ds[, j] <- analysis_out[[datasetnames[j]]]$analysis_groups$rank_sum
  }
  
  analysis_out$overall_bestmean_block <- which.max(apply(mean_all_ds, 1, mean)) # 5
  analysis_out$overall_bestranksum_block <- which.max(apply(ranksum_all_ds, 1, mean)) # 5
  
  for (ds in datasetnames) {
    analysis_out[[ds]][["overall_bestmean_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out$overall_bestmean_block]
    analysis_out[[ds]][["overall_bestranksum_values"]] <- analysis_out[[ds]]$analysis_values[, analysis_out$overall_bestranksum_block]
  }
  return(analysis_out)
}

# comparative analysis
get_comparative_analysis <- function(measure, results_set, CHIRPS_analysis) {
  
  for (ds in datasetnames) {
    
    # isolate and transform for one dataset
    analysis <- calc_wxcov(results_set, ds, sensitivity = FALSE)

    # transpose the raw values for analysis
    analysis_values <- tapply(analysis[, measure]
                              , analysis$algorithm
                              , identity)
    
    algos <- names(analysis_values) # no brl for cardio or nursery, no lending for brl or intrees

    analysis_values <- matrix(unlist(analysis_values), ncol = length(algos))
    analysis_values <- cbind(analysis_values, CHIRPS_analysis[[ds]]$ds_bestranksum_values)
    
    # collect
    algos <- c(algos, "CHIRPS")
    colnames(analysis_values) <- algos
    CHIRPS_analysis[[ds]][["comp_raw"]] <- analysis
    CHIRPS_analysis[[ds]][["comp_values"]] <- analysis_values
    
    # friedman test
    CHIRPS_analysis[[ds]][["comp_frd.tt"]] <- friedman.test(analysis_values)
    chisqstat <- CHIRPS_analysis[[ds]][["comp_frd.tt"]]$statistic
    df1 <- CHIRPS_analysis[[ds]][["comp_frd.tt"]]$parameter
    N <- nrow(analysis_values)
    fstat <- (chisqstat * (N - 1)) / (N * df1 - chisqstat)
    names(fstat) <- "Friedman F"
    df2 <- df1 * (N - 1)
    fpvalue <- pf(fstat, df1=df1, df2=df2, lower.tail = FALSE)
    
    # collect
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.statistic"]] <- fstat
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.p.value"]] <- fpvalue
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.N"]] <- N
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.df1"]] <- df1
    CHIRPS_analysis[[ds]][["comp_frd.tt"]][["F.df2"]] <- df2
  }
  return(CHIRPS_analysis)
}

results_in_detail <- function(analysis_in, rounding = 2, sgn = -1) {
  for (ds in datasetnames) {
    print(ds)
    print("mean qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, mean), rounding))
    print("sd qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, sd), rounding))
    print("min qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, min), rounding))
    print("max qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, max), rounding))
    print("med qm")
    print(round(apply(analysis_in[[ds]]$comp_values, 2, median), rounding))
    mr <- apply(apply(sgn * analysis_in[[ds]]$comp_values, 1, rank), 1, mean)
    print("mean ranks")
    print(round(mr, rounding))
    print("Fried test")
    fstat <- analysis_in[[ds]][["comp_frd.tt"]][["F.statistic"]]
    df1 <- analysis_in[[ds]][["comp_frd.tt"]][["F.df1"]]
    df2 <- analysis_in[[ds]][["comp_frd.tt"]][["F.df2"]]
    N <- analysis_in[[ds]][["comp_frd.tt"]][["F.N"]]
    f.p.value <- analysis_in[[ds]][["comp_frd.tt"]][["F.p.value"]]
    print(c(fstat, df1, df2, N, f.p.value))
    print("post hoc tests")
    algs <- names(mr)[names(mr) != "CHIRPS"]
    k <- analysis_in[[ds]]$comp_frd.tt$parameter + 1
    print("top and 2nd rank")
    top_rank <- which.min(mr)
    scnd_rank <- which.min(mr[-top_rank])
    print(c(top_rank, scnd_rank))
    print("rank diff")
    md <- mr[names(scnd_rank)] - mr[names(top_rank)]
    print(round(md, rounding))
    print("z stat")
    z <- (md) / sqrt((k * (k + 1)) / (6 * N))
    ztest <- pnorm(z, lower.tail = FALSE)
    print(z)
    print("post-hoc z-test")
    print(ztest)
    print("reject null bonferroni")
    print(ztest < 0.025/df1)
    print(0.025/df1)
    print("presence of ones")
    print(apply(analysis_in[[ds]]$comp_values
          , 2
          , function(x) {sum(ifelse(x == 1, TRUE, FALSE))}))
    print("presence of zeroes")
    print(apply(analysis_in[[ds]]$comp_values
          , 2
          , function(x) {sum(ifelse(x == 0, TRUE, FALSE))}))
  }
}

results_in_plots <- function(measure, analysis_in, rounding = 2, sgn = -1, y_lims = NA) {
  nds <- length(datasetnames)
  dmnms <- list(datasetnames, algs)
  crt_mat <- function() {
    matrix(NA
          , ncol = length(algs)
          , nrow = nds
          , dimnames = dmnms)
  }
  mn_meas <- crt_mat()
  sd_meas <- crt_mat()
  se_meas <- crt_mat()
  mr_meas <- crt_mat()
  N <- numeric(nds)
  for (i in seq_along(datasetnames)) {
    mns <- apply(analysis_in[[datasetnames[i]]]$comp_values, 2, mean)
    mn_meas[i, names(mns)] <- mns
    sds <- apply(analysis_in[[datasetnames[i]]]$comp_values, 2, sd)
    sd_meas[i, names(sds)] <- sds
    mrs <- apply(apply(sgn * analysis_in[[datasetnames[i]]]$comp_values, 1, rank), 1, mean)
    mr_meas[i, names(mrs)] <- mrs
    N[i] <- analysis_in[[datasetnames[i]]][["comp_frd.tt"]][["F.N"]]
  }
  
  dataset <- get_datasetname_stems(datasetnames)
  mn_meas <- as.data.frame(mn_meas)
  mn_meas$dataset <- dataset
  mn_meas <- mn_meas %>% gather(algorithm, mean, -dataset)
  se_meas <- sd_meas/sqrt(N)
  se_meas <- as.data.frame(se_meas)
  se_meas$dataset <- dataset
  se_meas <- se_meas %>% gather(algorithm, st.err, -dataset)
  mn_meas <- full_join(mn_meas, se_meas)
  # mn_meas$alpha <- ifelse(mn_meas$algorithm == "CHIRPS"
  #                         , "weight", "light")
  
  if (is.na(y_lims)) {
    sc_y <- scale_y_continuous()
  } else {
    sc_y <- scale_y_continuous(limits = y_lims)
  }
  g1 <- ggplot(data = mn_meas
              , aes(x = dataset, y = mean
                    , ymin = I(mean-st.err)
                    , ymax = I(mean+st.err)
                    , colour = algorithm
                    , group = algorithm
                    , alpha = algorithm)) +
    geom_line() +
    geom_point() +
    geom_errorbar(width = 0.2) +
    colos + 
    alpas +
    labs(y = get_measure_stem(measure)) +
    sc_y +
    myGgTheme
  
  return(g1)
}

sens_plot <- function(analysis_in) {
  
  # quirk of the way they are organised
  idx <- c(1:6, 13:18, 25:30, 7:12, 19:24, 31:36)
  idx <- c(idx, idx + 36)
  
  analysis_plotting <- analysis_in$analysis_groups
  analysis_plotting <- analysis_plotting[idx, ] # re-order
  analysis_plotting$id_plotting <- c(1:18, 20:37, 39:56, 58:75)
  
  tr_top_mean <- analysis_plotting[analysis_plotting$id == analysis_in$top_mean_block, "id_plotting"]
  arrowbase <- range(analysis_plotting$mean)[1]
  arrowtip <- arrowbase + (range(analysis_plotting$mean)[2]-range(analysis_plotting$mean)[1]) / 8
  
  g <- ggplot(
    data = analysis_plotting
    , aes(y = mean
          , ymin = I(mean-st.err)
          , ymax = I(mean+st.err)
          , x = id_plotting
          , shape = alpha
          , colour = func
          , size = support)) +
    geom_point() +
    geom_errorbar(size=0.5, width=1) +
    geom_vline(xintercept = c(19, 57)
               , linetype = 2
               , colour = myPalNeut[8]) +
    geom_vline(xintercept = 38
               , linetype = 1
               , colour = myPalNeut[7]) +
    geom_segment(x = tr_top_mean
                 , xend = tr_top_mean
                 , y = arrowbase
                 , yend = arrowtip
                 , arrow = arrow(length = unit(0.025, "npc"))
                 , colour = myPalNeut[7]
                 , size = 0.5) +
    labs(title = get_datasetname_stems(ds)
         , x = NULL, y = NULL) +
    scale_size_discrete(range = c(0.5, 1), guide = FALSE) +
    scale_shape(guide = FALSE) +
    sens_colos_silent +
    myGgTheme +
    theme(axis.text.x = element_blank()
          , panel.grid.major = element_blank()
          , panel.grid.minor = element_blank())
  return(g)
}

measure <- "stability.tr."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results)

for (ds in datasetnames) {
  print(c(CHIRPS_analysis[[ds]][["top_mean_block"]]
          , CHIRPS_analysis[[ds]][["top_ranksum_block"]])
  )
}

# before each measure, get the top rank mean from stability.tr.
measure <- "stability.tt."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = "stability.tr."
                                       , top_ranksum_block = "stability.tr.")
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)
stab_plot <- results_in_plots(measure, tt_analysis, rounding = 4, y_lims = c(0,1))

measure <- "precision.tt."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = "stability.tr."
                                       , top_ranksum_block = "stability.tr.")
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)
prec_plot <- results_in_plots(measure, tt_analysis, rounding = 4, y_lims = c(0,1))

tikz(file = "prec_stab.tikz", width = 6.85, height = 5)
grid.arrange(prec_plot, stab_plot, nrow = 2)
dev.off()

measure <- "coverage.tt."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = "stability.tr."
                                       , top_ranksum_block = "stability.tr.")
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)
cov_plot <- results_in_plots(measure, tt_analysis, rounding = 4, y_lims = c(0,1))

measure <- "wxcoverage.tt."
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = "stability.tr."
                                       , top_ranksum_block = "stability.tr.")
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4)
wxcov_plot <- results_in_plots(measure, tt_analysis, rounding = 4, y_lims = c(0,1))


tikz(file = "cov_wxcov.tikz", width = 6.85, height = 5)
grid.arrange(cov_plot, wxcov_plot, nrow = 2)
dev.off()

measure <- "rule.length"
CHIRPS_analysis <- get_CHIRPS_analysis(measure, sens_results
                                       , top_mean_block = "stability.tr."
                                       , top_ranksum_block = "stability.tr.")
tt_analysis <- get_comparative_analysis(measure, comp_results, CHIRPS_analysis)
results_in_detail(tt_analysis, rounding = 4, sgn = 1) # positive ranking. shortest wins
results_in_plots(measure, tt_analysis, rounding = 4)

mean_rulelens <- matrix(NA, ncol= 5, nrow = 9)
dimnames(mean_rulelens) <- list(datasetnames
                                , c(algorithms, "CHIRPS"))
for (ds in datasetnames) {
  mn_rl <- apply(tt_analysis[[ds]]$comp_values, 2, mean)
  mean_rulelens[ds, ] <- mn_rl[c(algorithms, "CHIRPS")]
}

measure <- "mean_rule_cascade"
analysis <- comp_summ_results
# transpose the raw values for analysis
analysis_values <- tapply(analysis[, measure]
                          , list(analysis$dataset_name
                                 , analysis$algorithm)
                          , identity)

mean_rulelens
analysis_values <- cbind(analysis_values, CHIRPS = 1)
analysis_values
mean_rulelens * analysis_values

measure <- "fidelity"
analysis <- comp_summ_results
# transpose the raw values for analysis
analysis_values <- tapply(analysis[, measure]
                          , list(analysis$dataset_name
                                 , analysis$algorithm)
                          , identity)
analysis_values


for (ds in datasetnames) {
  A_rules <- with(tt_analysis[[ds]]$comp_raw
                  , tt_analysis[[ds]]$comp_raw[algorithm == "Anchors"
                                              , c("pretty.rule", "rule.length")])

  B_rules <- with(tt_analysis[[ds]]$comp_raw
                  , tt_analysis[[ds]]$comp_raw[algorithm == "BRL"
                                              , c("pretty.rule", "rule.length")])

  D_rules <- with(tt_analysis[[ds]]$comp_raw
                  , tt_analysis[[ds]]$comp_raw[algorithm == "defragTrees"
                                              , c("pretty.rule", "rule.length")])
  
  I_rules <- with(tt_analysis[[ds]]$comp_raw
                  , tt_analysis[[ds]]$comp_raw[algorithm == "inTrees"
                                              , c("pretty.rule", "rule.length")])
  
  C_rules <- with(tt_analysis[[ds]]$sens_raw
                  , tt_analysis[[ds]]$sens_raw[, c("pretty.rule", "rule.length")])
  
  uniques <- c(length(unique(A_rules$pretty.rule))
               , length(unique(B_rules$pretty.rule))
               , length(unique(D_rules$pretty.rule))
               , length(unique(I_rules$pretty.rule))
               , length(unique(C_rules$pretty.rule)))
  
  names(uniques) <- c("A", "B", "D", "I", "C")
  print(ds)
  print(uniques)
}

res <- (gather(as.data.frame(tt_analysis$adult_small_samp$comp_values)
               , "algorithm", measure))
densityplot(~measure
            , data = res
            , groups = algorithm
            , col=myPal[c(1,2,5,3,4)]
            , alpha = c(0.4, 0.4, 1, 0.4, 0.4)
            , auto.key = list(columns=5
                              , col=myPal[c(1,2,5,3,4)]
                              , alpha = c(0.4, 0.4, 1, 0.4, 0.4)
                              , lines=FALSE)
            , plot.points=FALSE
            , par.settings = MyLatticeTheme)

res <- calc_wxcov(sens_results, "adult_small_samp", sensitivity = TRUE)

# sensitivity plots
for(ds in datasetnames) {
  tikz(file = paste0(ds, "_sens.tikz"), width = 2.2, height = 1.5)
  print(sens_plot(CHIRPS_analysis[[ds]]))
  dev.off()
}  

# plot to get only the legend
legend_plot <- ggplot(
  data = analysis_in
  , aes(y = mean
        , ymin = I(mean-st.err)
        , ymax = I(mean+st.err)
        , x = id_plotting
        , shape = alpha
        , colour = func
        , size = support)) +
  geom_point() +
  geom_errorbar(size=0.5, width=1) +
  geom_vline(xintercept = c(19, 57)
             , linetype = 2
             , colour = myPalNeut[8]) +
  geom_vline(xintercept = 38
             , linetype = 1
             , colour = myPalNeut[7]) +
  labs(title = get_datasetname_stems(ds)
       , x = NULL, y = NULL) +
  scale_size_discrete(range = c(0.5, 1)) +
  sens_colos +
  myGgTheme +
  theme(axis.text.x = element_blank()
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , legend.box = "horizontal")

lgnd <- cowplot::get_legend(legend_plot)
grid.newpage()
tikz(file = "legend_sens.tikz", width = 2.2, height = 1.5)
grid.draw(lgnd)
dev.off()

analysis_in <- CHIRPS_analysis[[ds]]


