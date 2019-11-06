library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(tikzDevice)
library(cowplot) # get legend out of a ggplot

algorithms <- c("Anchors", "lore", "CHIRPS")

# plot themes
source("C:\\Users\\id126493\\OneDrive\\Documents\\PhD\\KTheme.R")
# source("C:\\Users\\Crutt\\OneDrive\\Documents\\PhD\\KTheme.R")
reds <- k.grad.red.rev(2)
ongs <- k.grad.orange.rev(2)
blus <- k.grad.blue.rev(2)
grns <- k.grad.green.rev(2)
prps <- k.grad.purple.rev(2)
nuts <- myPalNeut

myPal <- c(blus[2], grns[2], k.pink)
myAlph <- c(rep(0.4, length(algorithms) - 1), 1)
names(myPal) <- algorithms
names(myAlph) <- algorithms
colos <- scale_colour_manual(
  values = myPal)
alpas <- scale_alpha_manual(values = myAlph)

# data management
project_dir <- "V:\\whiteboxing\\BMC2\\"
# project_dir <- "C:\\Users\\Crutt\\Documents\\whiteboxing\\BMC\\"
# project_dir <- "/home/julian/whiteboxing/"

n_classes <- c(breast = 2
               , cardio = 3
               , diaretino = 2
               , heart = 3
               , readmit = 2
               , thyroid = 2
               , mhtech14 = 2
               , mh2tech16 = 2
               , usoc2 = 2
)

datasetnames <- names(n_classes)

class_cols <- c(
  "mb"
  , "NSP"
  , "dr"
  , "HDisease"
  , "readmitted"
  , "diagnosis"
  , "treatment"
  , "mh2"
  , "mh"
)

names(class_cols) <- names(n_classes)

resfilesdirs <- paste0(project_dir, datasetnames, "\\")

# comparative analysis
patt <- paste0("(", paste(algorithms, collapse = ")|("), ")")

train_set_size <- integer(length(datasetnames))
test_set_size <- integer(length(datasetnames))
names(train_set_size) <- datasetnames
names(test_set_size) <- datasetnames

first_comp <- TRUE
first_comp_summ <- TRUE

for (i in seq_along(resfilesdirs)) {
  datasetname <- datasetnames[i]
  filepath <- normalizePath(file.path(resfilesdirs[i]))
  filenames <- dir(filepath)
  for (fn in filenames) {
    if (fn == "y_test.csv") {
      test_set_size[datasetname] <- nrow(read.csv(paste0(filepath, fn)))
    }
    if (fn == "y_train.csv") {
      train_set_size[datasetname] <- nrow(read.csv(paste0(filepath, fn)))
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

precs <- with(comp_results, tapply(precision.tt.
                          , list(comp_results$dataset
                                 , comp_results$algorithm)
                          , mean))

if (!("lore" %in% dimnames(precs)[[2]])) {
  precs <- cbind(precs, lore = NA)
}

stabs <- with(comp_results, tapply(stability.tt.
                          , list(comp_results$dataset
                                 , comp_results$algorithm)
                          , mean))

if (!("lore" %in% dimnames(stabs)[[2]])) {
  stabs <- cbind(stabs, lore = NA)
}


covs <- with(comp_results, tapply(coverage.tt.
                          , list(comp_results$dataset
                                 , comp_results$algorithm)
                          , mean))

if (!("lore" %in% dimnames(covs)[[2]])) {
  covs <- cbind(covs, lore = NA)
}


wxcovs <- with(comp_results, tapply(wxcoverage.tt.
                          , list(comp_results$dataset
                                 , comp_results$algorithm)
                          , mean))

if (!("lore" %in% dimnames(wxcovs)[[2]])) {
  wxcovs <- cbind(wxcovs, lore = NA)
}


for (ds in datasetnames) {
  print(ds)
  cres <- comp_results %>% filter(dataset == ds) %>%
    dplyr::select(instance_id, dataset, stability.tr., algorithm)
  cres <- spread(cres, algorithm, stability.tr.)
  
  if (is.na(stabs[ds, "lore"])) {
    
    print("mean ranks")
    mrs <- apply(apply(-cres[, algorithms[c(1, 3)]], 1, rank), 1, mean)
    print(mrs)
    cres_ttest <- t.test(cres[, algorithms[1]], cres[, algorithms[3]], paired = TRUE)
    cres_wxtest <- wilcox.test(cres[, algorithms[1]], cres[, algorithms[3]], paired = TRUE)
    print(cres_ttest)
    print(cres_wxtest)
  
    } else {
    
    print("mean ranks")
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
    
    cres_nem <- posthoc.friedman.nemenyi.test(as.matrix(cres[, algorithms]))
    print(cres_nem)  
    
  }
  
  
}
