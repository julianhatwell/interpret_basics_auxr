library(tidyr)
library(dplyr)
library(lattice)
library(ggplot2)
library(rjson)
library(modeest)

# main results collection
folders <- c("adult_small_samp"
             , "bankmark_samp"
             , "car"
             , "cardio"
             , "credit"
             , "german"
             , "lending_tiny_samp"
             , "nursery_samp"
             , "rcdv_samp")

meths <- c("Anchors", "BRL", "CHIRPS", "defragTrees", "inTrees")

metrics_names <- c("precision.tt.", "stability.tt.", "recall.tt."
                   , "f1.tt.", "accuracy.tt.", "lift.tt."
                   , "coverage.tt.", "xcoverage.tt.", "kl_div.tt.")
metrics <- sub(".tt.", "", metrics_names)


np <- function(fld) {
  return(function(fl) {
    normalizePath(file.path(resfilesdir, fld, fl))  
  })
}

is_which_method <- function(x) {
  for (meth in meths) {
    res <- grepl(meth, x)
    if (res) return(meth)
  }
  return(NA)
}

is_forestperf_file <- function(x) {
  grepl("forest_performance", fl)
}

is_bestparams_file <- function(x) {
  grepl("best_params", fl)
}

get_rnst_from_filename <- function(x) {
  pos <- regexpr("rnst_", x)[1]
  rnst <- substr(x
                 , start = pos + 5
                 , stop = pos + 7)
}

# temp - correcting for mismatch in file names
for(fld in folders) {
  normpath <- np(fld)
  
  fls <- dir(normpath(""))
  for (fl in fls) {
    if (grepl("results_", fl)) {
      file.rename(normpath(fl), normpath(sub("results_", "", fl)))
    }
  }
}

# temp - correcting for mismatch in file names
for(fld in folders) {
  normpath <- np(fld)
  
  fls <- dir(normpath(""))
  for (fl in fls) {
    if (grepl("rndst_", fl)) {
      file.rename(normpath(fl), normpath(sub("rndst_", "rnst_", fl)))
    }
  }
}

results_mat <- function(d, m) {
  tapply(d[[m]], list(d$dataset_name
                 , d$algorithm
                 , d$random_state)
         , mean)
}

rfperf_mat <- function(d, m, func) {
    tapply(d[[m]], list(d$dataset_name)
         , func)
}

mean_narm <- function(x) {
  mean(x, na.rm = TRUE)
}

sd_narm <- function(x) {
  sd(x, na.rm = TRUE)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

prep_res <- function(x, func) {
  res <- apply(x, c(1, 2), func)
  res <- as.data.frame(res)
  res$datasetname <- rownames(res)
  rownames(res) <- NULL
  return(res)
}

mean_df <- function(x) {
  res <- prep_res(x,  mean_narm)
  mean_res <- gather(res
                     , key="method", value = "mn", -datasetname)
  res <- prep_res(x,  sd_narm)
  sd_res <- gather(res
                   , key="method", value = "sd", -datasetname)
  
  res <- inner_join(mean_res, sd_res)
  res$upr <- res$mn + 1.96 * res$sd
  res$lwr <- res$mn - 1.96 * res$sd
  
  res$datasetname <- factor(res$datasetname)
  res$method <- factor(res$method)
  return(res)
}

my_dotplot <- function(dat) {
  pd <- position_dodge(0.5)
  g <- ggplot(data=dat, aes(x=method, y=mn
                            , ymin = lwr, ymax=upr
                            , colour=datasetname))
  g <- g + geom_point(position=pd)
  g <- g + geom_errorbar(position=pd, width=0.5)
  g <- g + facet_wrap(~qmeasure)
  g <- g + labs(y = "mean")
  g <- g + theme_bw()
  return(g)
}

my_histo <- function(dat) {
  g <- ggplot(data = dat
              , aes(x=ntree)) +
    stat_bin(bins=8) +
    theme(axis.title = element_blank()
                 , axis.text = element_blank()
                 , axis.ticks = element_blank()
                 , panel.grid.major = element_blank()
                 , panel.grid.minor = element_blank()
                 , panel.background = element_blank())
  return(g)
}
