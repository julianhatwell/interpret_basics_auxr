# previously completed results
load("shape_measures_30032018.RData")
shape_measures <- shape_measures[c("n_row", "asso", "prnc")]
names(shape_measures) <- c("n_row_full", "asso_full", "prnc_full")

load("complexity_measures_30032018.RData")
complexity_measures <- cbind(shape_measures, complexity_measures)
complexity_measures <- complexity_measures[c(1, 4:8, 2, 9, 3, 10:31)]

save(complexity_measures, file="complexity_measures.RData")
