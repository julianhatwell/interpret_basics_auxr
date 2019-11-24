library(dplyr)
library(tidyr)
library(PMCMRplus)
options(max.print=20*72)
algorithms <- c("Anchors", "BRL", "defragTrees", "inTrees", "CHIRPS") # 

# plot themes
source("KTheme.R")
# source("C:\\Users\\id126493\\OneDrive\\Documents\\PhD\\KTheme.R")
# source("C:\\Users\\Crutt\\OneDrive\\Documents\\PhD\\KTheme.R")
reds <- k.grad.red.rev(2)
ongs <- k.grad.orange.rev(2)
blus <- k.grad.blue.rev(2)
grns <- k.grad.green.rev(2)
prps <- k.grad.purple.rev(2)
nuts <- myPalNeut

myPal <- c(reds[2], ongs[2], blus[2], grns[2], k.pink)
myAlph <- c(rep(0.4, length(algorithms) - 1), 1)
names(myPal) <- algorithms
names(myAlph) <- algorithms
colos <- scale_colour_manual(
  values = myPal)
alpas <- scale_alpha_manual(values = myAlph)


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
    sens_groups$lwr_mean_ci <- mn - sens_groups$st.err
    sens_groups$mean <- mn
    sens_groups$upr_mean_ci <- mn + sens_groups$st.err
    qntl <- apply(sens_values, 2, quantile)
    sens_groups$lwrq <- qntl[2, ]
    sens_groups$median <- apply(sens_values, 2, median)
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

stability_sens <- get_sensitivity("stability.tt.")
wxcoverage_sens <- get_sensitivity("wxcoverage.tt.")
rule.length_sens <- get_sensitivity("rule.length")

get_best_sens(stability_sens)

# ds <- "cardio"
# stability_sens[[ds]]
# wxcoverage_sens[[ds]]
# rule.length_sens[[ds]]
# 
# which.min(sensitivity_analysis[["german"]][["rank_mean"]])
# ph[[3]][which.min(ph[[3]])]
# 
# cn <- apply(select(sens_groups, support, alpha, bins, func, weights), 1, paste, collapse = "_")
# dimnames(sens_values) <- list(NULL, cn)
# rowMeans(apply(-cbind(sens_values[, cn[1:(length(cn) -1)][which.min(ph[[3]]) %/% (length(cn) -1)]],
# sens_values[, cn[2:length(cn)][which.min(ph[[3]]) %% (length(cn) -1)]]), 1, rank))
# 
# apply(sens_values, 2, mean)
# ds <- "car"
# sensitivity_analysis[[ds]][sensitivity_analysis[[ds]]$rank_mean == min(sensitivity_analysis[[ds]]$rank_mean), ] 


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

meas <- "stability.tt."
st_err <- get_sd_of(meas) / test_set_size_sqrt
round(get_mean_of(meas), 4)
round(st_err, 4)
get_mean_ranks_of(meas)

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

get_lwrq_of("precision.tt.")
get_median_of("precision.tt.")
get_uprq_of("precision.tt.")

get_lwrq_of("stability.tt.")
get_median_of("stability.tt.")
get_uprq_of("stability.tt.")



cres <- dplyr::select(comp_results, instance_id, dataset, !! qmeas, algorithm) %>%
  mutate(dataset = factor(dataset), algorithm = factor(algorithm))
ggplot(data = cres, aes(y = stability.tt.
                        , colour = algorithm)) +
  geom_boxplot() +
  facet_wrap(~dataset) +
  myGgTheme

get_meds_for_plotting <- function(meas) {
  
  meds <- as.data.frame(get_median_of(meas))
  meds$dataset <- rownames(meds)
  meds <- gather(meds, algorithm, m, -dataset)
  lwrqs <- as.data.frame(get_lwrq_of(meas))
  lwrqs$dataset <- rownames(lwrqs)
  lwrqs <- gather(lwrqs, algorithm, lwr, -dataset)
  uprqs <- as.data.frame(get_uprq_of(meas))
  uprqs$dataset <- rownames(uprqs)
  uprqs <- gather(uprqs, algorithm, upr, -dataset)
  meds$lwr <- lwrqs$lwr
  meds$upr <- uprqs$upr
  return(meds)
  
}

get_means_for_plotting <- function(meas) {
  
  means <- as.data.frame(get_mean_of(meas))
  means$dataset <- rownames(means)
  means <- gather(means, algorithm, m, -dataset)
  
  st_errs <- as.data.frame(get_sd_of(meas) / test_set_size_sqrt)
  st_errs$dataset <- rownames(st_errs)
  st_errs <- gather(st_errs, algorithm, st_err, -dataset)
  
  means$lwr <- means$m - st_errs$st_err
  means$upr <- means$m + st_errs$st_err
  return(means)
}

lwrmedupr <- get_meds_for_plotting("stability.tt.")
lwrmnupr <- get_means_for_plotting("stability.tt.")

ggplot(data = get_means_for_plotting("stability.tt.")
       , aes(y = m
             , ymin = lwr
             , ymax = upr
             , x = dataset
             , colour = algorithm
             , group = algorithm
             , alpha = algorithm)) +
  myGgTheme +
  geom_line() +
  geom_point() +
  geom_errorbar(width = 0.2) +
  colos + 
  alpas

meas <- "stability.tt."
for (ds in datasetnames) {
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
}



# to do list
# tabulate for each dataset - mean & mean rank; prec, stab, cov, xcov
# plot for each of prec, stab, cov, xcov, line through mean or median with boxplot superimposed
# fidelity for three global methods (n - number of zeros) / n
# plot elapsed time on log scale
# settings used for each run
# sensitivity analysis plots (boxplots)
# note that BRL requires pre-discretised data