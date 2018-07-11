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
  filenames <- grep("_addt_.*\\.csv", dir(normalizePath(file.path(resfilesdir, exp))), value = TRUE)
  for (filename in filenames) {
    datasetname <- sub("\\_pickles", "", exp)
    randst <- gregexpr("rnst\\_1[2-5][0-9]", filename)
    randst <- regmatches(filename, randst)[[1]]
    randst <- as.numeric(gsub("rnst\\_", "", randst))
    
    addt <- gregexpr("addt\\_[0-9]{1,4}", filename)
    addt <- regmatches(filename, addt)[[1]]
    addt <- as.numeric(gsub("addt\\_", "", addt))
    
    results <- read.csv(normalizePath(file.path(resfilesdir, exp, filename)), stringsAsFactors = FALSE)
    
    datasetname <- rep(datasetname, nrow(results))
    random_state <- rep(randst, nrow(results))
    add_trees <- rep(addt, nrow(results))
    
    results <- cbind(results, datasetname, random_state, addt)

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
  
  result_set <- NULL
  
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

names(main_results) <- gsub("total.coverage", "excl.cov", gsub("precision", "stability", names(main_results)))