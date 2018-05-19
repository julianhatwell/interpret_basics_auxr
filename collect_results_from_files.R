source("dirs.R")
dir(resultsfilesdir)

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

# set up the factor vars
main_results$instance_id <- factor(main_results$instance_id)
main_results$result_set <- ifelse(main_results$result_set == "greedy_prec", "CHIPS", main_results$result_set)
main_results$datasetname <- gsub("\\_samp", "", main_results$datasetname)
main_results$datasetname <- gsub("\\_small", "", main_results$datasetname)
main_results$datasetname <- gsub("\\_tiny", "", main_results$datasetname)
main_results$result_set <- factor(main_results$result_set)
main_results$pred.class <- factor(main_results$pred.class)
main_results$pred.class.label <- factor(main_results$pred.class.label)
main_results$target.class <- factor(main_results$target.class)
main_results$target.class.label <- factor(main_results$target.class.label)
main_results$random_state <- factor(main_results$random_state)