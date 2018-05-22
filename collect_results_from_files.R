source("dirs.R")
library(xlsx)

# main results collection
exps <- c("adult_small_samp_pickles"
          , "bankmark_samp_pickles"
          , "car_pickles"
          , "cardiotography_pickles"
          , "credit_pickles"
          , "german_pickles"
          , "lending_tiny_samp_pickles"
          , "nursery_samp_pickles"
          , "rcdv_samp_pickles")

count <- 1
for (exp in exps) {
  filenames <- grep("pr\\_0\\.9[59]?\\.csv", dir(normalizePath(file.path(resfilesdir, exp))), value = TRUE)
  for (filename in filenames) {
    datasetname <- sub("\\_pickles", "", exp)
    randst <- gregexpr("rnst\\_12[34567]", filename)
    randst <- regmatches(filename, randst)[[1]]
    randst <- as.numeric(gsub("rnst\\_", "", randst))
    supp <- gregexpr("sp\\_0.[0-9]{1,2}", filename)
    supp <- regmatches(filename, supp)[[1]]
    supp <- as.numeric(gsub("sp\\_", "", supp))
    alpp <- gregexpr("ap\\_0.[1-9]", filename)
    alpp <- regmatches(filename, alpp)[[1]]
    alpp <- as.numeric(gsub("ap\\_", "", alpp))
    alsc <- gregexpr("as\\_0.[0-9]{1,2}", filename)
    alsc <- regmatches(filename, alsc)[[1]]
    alsc <- as.numeric(gsub("as\\_", "", alsc))
    prpp <- gregexpr("pr\\_0\\.9[59]?", filename)
    prpp <- regmatches(filename, prpp)[[1]]
    prpp <- as.numeric(gsub("pr\\_", "", prpp))
    
    results <- read.csv(normalizePath(file.path(resfilesdir, exp, filename)), stringsAsFactors = FALSE)
    
    datasetname <- rep(datasetname, nrow(results))
    random_state <- rep(randst, nrow(results))
    support <- rep(supp, nrow(results))
    alpha_paths <- rep(alpp, nrow(results))
    alpha_scores <- rep(alsc, nrow(results))
    target_prec <- rep(prpp, nrow(results))
    
    results <- cbind(results, datasetname, as.data.frame(cbind(random_state, support, alpha_paths, alpha_scores, target_prec)))

    if (count == 1) {
      main_results <- results
    } else {
      main_results <- rbind(main_results, results)
    }
    count <- count + 1
  }
}

# tidy up
main_results <- within(main_results, {
  instance_id <- factor(instance_id)
  result_set <- ifelse(result_set == "greedy_prec", "CHIPS", result_set)
  rset_supp <- ifelse(result_set == "CHIPS", paste0("CHIPS_", support), result_set)
  
  result_set <- factor(result_set)
  rset_supp <- factor(rset_supp)
  
  # tidy names
  datasetname <- gsub("\\_samp", "", datasetname)
  datasetname <- gsub("\\_small", "", datasetname)
  datasetname <- gsub("\\_tiny", "", datasetname)
  datasetname <- factor(datasetname)
  pred.class <- factor(pred.class)
  pred.class.label <- factor(pred.class.label)
  target.class <- factor(target.class)
  target.class.label <- factor(target.class.label)
  random_state <- factor(random_state)
  
})

main_results[main_results$result_set == "anchors", c("majority.vote.share", "support", "alpha_paths", "alpha_scores")] <- NA


names(main_results) <- gsub("total.coverage", "excl.cov", gsub("precision", "stability", names(main_results)))

# results timings
time_results <- read.xlsx(normalizePath(file.path(resfilesdir, "base model stats.xlsx"))
          , sheetIndex = 1)

time_results <- within(time_results, {
  datasetname <- factor(datasetname)
  randst <- factor(randst)
})

# support timings
support_results <- read.xlsx(normalizePath(file.path(resfilesdir, "support sensitivity.xlsx"))
                          , sheetIndex = 1)

support_results <- within(support_results, {
  datasetname <- factor(datasetname)
})

    