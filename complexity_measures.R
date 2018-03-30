library(dplyr)
source("dataset_structures.R")

dsets <- data.frame(
  dnames = c("accident", "adult", "bankmark", "car", "cardiotography", "credit" , "german", "lending", "nursery", "rcdv")
  , class_cols = c("Accident_Severity", "income", "y", "acceptability", "NSP", "A16", "rating", "loan_status", "decision", "recid")
  , stringsAsFactors = FALSE
  )

# full collection of measures from subsampled data sets
source("general_setup.R")
for (i in 1:nrow(dsets)) {
  dname <- dsets$dnames[i]
  class_col <- dsets$class_cols[i]
  meas <- dataset_structure(dname, class_col)
  meas_names <- names(meas)
  
  if (i == 1) {
    # first time use
    complexity_measures <- list()
    for(name in meas_names) {
      complexity_measures[[name]] <- meas[name]
    }
    complexity_measures <- as.data.frame(complexity_measures)
    rownames(complexity_measures) <- dname
  } else {
    # subequent uses
    complexity_measures <- rbind(complexity_measures, meas)
    rownames(complexity_measures)[i] <- dname
  }
}
complexity_measures
save(complexity_measures, file="complexity_measures.RData")

# just shape measures from subsampled data sets
source("general_setup_full_size.R")
for (i in 1:nrow(dsets)) {
  dname <- dsets$dnames[i]
  class_col <- dsets$class_cols[i]
  meas <- dataset_structure(dname, class_col, what = "shp")
  meas_names <- names(meas)
  
  if (i == 1) {
    # first time use
    complexity_measures <- list()
    for(name in meas_names) {
      complexity_measures[[name]] <- meas[name]
    }
    complexity_measures <- as.data.frame(complexity_measures)
    rownames(complexity_measures) <- dname
  } else {
    # subequent uses
    complexity_measures <- rbind(complexity_measures, meas)
    rownames(complexity_measures)[i] <- dname
  }
}
shape_measures <- complexity_measures
save(shape_measures, file="shape_measures.RData")
