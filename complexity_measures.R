# source("general_setup.R")
# source("dataset_structures.R)

dsets <- data.frame(
  dnames = c("car", "cardiotography", "credit" , "german", "nursery", "adult", "rcdv", "lending")
  , class_cols = c("acceptability", "NSP", "A16", "rating", "decision", "income", "recid", "loan_status")
  , stringsAsFactors = FALSE
  )


# subequent uses
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
    complexity_measures <- rbind(complexity_measures, meas)
    rownames(complexity_measures)[i] <- dname
  }
}

complexity_measures
