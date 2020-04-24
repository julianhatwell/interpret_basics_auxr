library(tidyr)
library(dplyr)
library(lattice)
library(ggplot2)

# main results collection
folders <- c("adult"
             , "bankmark"
             , "car"
             , "cardio"
             , "credit"
             , "german"
             , "lending_tiny_samp"
             , "nursery"
             , "rcdv")

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

results_mat <- function(x) {
  tapply(x, list(main_results$dataset_name
                 , main_results$algorithm
                 , main_results$random_state)
         , mean)
}

mean_narm <- function(x) {
  mean(x, na.rm = TRUE)
}

sd_narm <- function(x) {
  sd(x, na.rm = TRUE)
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

my_dotplot <- function(dat, metric) {
  pd <- position_dodge(0.5)
  g <- ggplot(data=dat, aes(x=method, y=mn
                            , ymin = lwr, ymax=upr
                            , colour=datasetname))
  g <- g + geom_point(position=pd)
  g <- g + geom_errorbar(position=pd, width=0.5)
  g <- g + labs(y = metric)
  g <- g + theme_bw()
  return(g)
}
