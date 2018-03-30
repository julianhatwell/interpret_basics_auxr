# previously completed results
load("shape_measures_30032018.RData")
shape_measures <- shape_measures[c("n_row", "asso", "prnc")]
names(shp_measures) <- c("n_row_full", "asso_full", "prnc_full")

load("complexity_measures_30032018.RData")
complexity_measures <- cbind(shp_measures, complexity_measures)
complexity_measures <- complexity_measures[c(1, 4:8, 2, 9, 3, 10:31)]
