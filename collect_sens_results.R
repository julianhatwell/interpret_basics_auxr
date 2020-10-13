rm(list=ls())

library(dplyr)
library(tidyr)
library(PMCMRplus)
library(cowplot)
library(rlang)
options(max.print=20*72)
algorithms <- c("Anchors", "BRL", "CHIRPS", "defragTrees", "inTrees", "lore")

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

maindir <- "sens"
sensdir <- "ada2_sensitivity"

resfilesdirs <- paste0(project_dir, maindir, pathsep, datasetnames, pathsep)
sensdirs <- paste0(resfilesdirs, sensdir, pathsep)
# sensdirs <- paste0(resfilesdirs, "rf_sensitivity", pathsep)

# sensitivity analysis
first_pass <- TRUE
for (i in seq_along(sensdirs)) {
  
  filepath <- normalizePath(file.path(sensdirs[i]))
  filenames <- dir(filepath)
  for (filename in filenames) {
    
    if (!grepl("summary", filename)) {
      wcts <- gregexpr("wcts\\_(conf\\_weighted|majority|targetclass|signdelta)", filename)
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
    geom_errorbar(size=0.255, width=1) +
    geom_segment(x = tr_top_mean
                 , xend = tr_top_mean
                 , y = arrowbase
                 , yend = arrowtip
                 , arrow = arrow(length = unit(0.025, "npc"))
                 , colour = myPalNeut[7]
                 , size = 0.5) +
    labs(title = ds # get_datasetname_stems(ds)
         , x = NULL, y = NULL) +
    scale_size_discrete(range = c(0.7, 1.3), guide = FALSE) +
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
  
  tikz(file = paste0(
    ds # get_datasetname_stems(ds)
    , "_sens.tikz"), width = 2.2, height = 1.5)
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
  scale_size_discrete(range = c(0.7, 1.3), labels = c("0.1", "0.2")) +
  sens_colos +
  theme(legend.box = "horizontal"
        , legend.key = element_blank())

lgnd <- get_legend(legend_plot)
grid.newpage()
tikz(file = "legend_sens.tikz", width = 2.2, height = 1.5)
grid.draw(lgnd)
dev.off()