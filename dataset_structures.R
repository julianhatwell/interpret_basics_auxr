dataset_structure <- function(dname, class_col) {
  
  library(ECoL)
  library(nFactors)
  library(caret)
  library(vcd)
  
  
  dat <- get(dname)
  
  fmla <- as.formula(paste(class_col, "~ ."))
  
  cats <- sapply(dat, class)
  facs <- names(dat)[cats == "factor" & names(dat) != class_col]
  nums <- names(dat)[cats %in% c("numeric", "integer") & names(dat) != class_col]

  asso <- 0
  for (f in facs) {
    if (assocstats(table(dat[[f]], dat[[class_col]]))$chisq_test[2,3] <= 0.05) {
      asso <- asso + 1
    }
  }
  
  prnc <- 0
  if (length(nums) > 0)
    if (length(nums) == 1) {
      prnc <- 1
    } else {
      preproc <- preProcess(dat[nums], method = c("BoxCox", "center", "scale"))
      trans <- predict(preproc, dat[nums])
      prnc <- max(nScree(cor(trans))$components
                  , sum(eigen(cor(trans))$values >= 1.0))
    }

  shape <- c(n_row = nrow(dat)
             , n_col = ncol(dat)
             , n_facs = length(facs)
             , n_num = length(nums)
             , asso = asso
             , prnc = prnc)
  
  bal <- tryCatch(balance(fmla, dat)
                  , error = function(e) return(c(C1 = NA, C2 = NA)))
  print(paste("balance measures:", dname))
  print(bal)
  
  dmn <- tryCatch(dimensionality(fmla, dat)
                  , error = function(e) return(c(T2 = NA, T3 = NA, T4 = NA)))
  print(paste("dimensionality measures:", dname))
  print(dmn)
  
  # takes a long time to run
  lnr <- tryCatch(linearity(fmla, dat)
                  , error = function(e) return(c(L1 = NA, L2 = NA, L3 = NA)))
  print(paste("linearity measures:", dname))
  print(lnr)
  
  ovr <- tryCatch(overlapping(fmla, dat)
                  , error = function(e) return(F1 = NA, F1v = NA, F2 = NA, F3 = NA, F4 = NA))
  print(paste("overlapping measures:", dname))
  print(ovr)
  
  # requires very large vector
  nei <- tryCatch(neighborhood(fmla, dat)
                  , error = function(e) return(c(N1 = NA, N2 = NA, N3 = NA, N4 = NA, T1 = NA, LSCAvg = NA)))
  print(paste("neighbourhood measures:", dname))
  print(nei)
  
  # requires very large vector
  nwk <- tryCatch(network(fmla, dat)
                  , error = function(e) return(c(Density = NA, ClsCoef = NA, Hubs = NA)))
  print(paste("network measures:", dname))
  print(nwk)
  
  return(c(shape, bal, dmn, lnr, ovr, nei, nwk))
}