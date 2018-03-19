library(ECoL)
library(nFactors)
# source("general_setup.R")
# now check the complexity measures

# car
print(paste("Complexity measures of adult"))
bal <- balance(acceptability~., data=car)
print(bal)

dmn <- dimensionality(acceptability~., data=car)
print(dmn)

# takes a long time to run
lnr <- linearity(acceptability~., data=car)
print(lnr)

ovr <- overlapping(acceptability~., data=car)
print(ovr)

# requires very large vector
nei <- neighborhood(acceptability~., data=car)
print(nei)

# requires very large vector
nwk <- network(acceptability~., data=car)
print(nwk)

class_col <- "acceptability"
cats <- sapply(car, class)
facs <- names(car)[cats == "factor" & names(car) != class_col]
nums <- names(car)[cats %in% c("numeric", "integer") & names(car) != class_col]

asso <- 0
for (f in facs) {
  if (assocstats(table(car[[f]], car[[class_col]]))$chisq_test[2,3] <= 0.05) {
    asso <- asso + 1
  }
}

prnc <- 0
if (length(nums) > 0)
  if (length(nums) == 1) {
    prnc <- 1
  } else {
    # do princomp analysis
    # return optimum number
  }

shape <- c(n_row = nrow(car)
           , n_col = ncol(car)
           , n_facs = length(facs)
           , n_num = length(nums)
           , asso = asso
           , prnc = prnc)

meas <- c(shape, bal, dmn, lnr, ovr, nei, nwk)
meas_names <- names(meas)

complexity_measures <- list()
for(name in meas_names) {
  complexity_measures[[name]] <- meas[name]
}
complexity_measures <- as.data.frame(complexity_measures)
rownames(complexity_measures) <- "car"

# cardiotography
bal <- balance(NSP~., data = cardiotography) # C2 is best used for multiclass problems
print(bal)

dmn <- dimensionality(NSP~., data = cardiotography)
print(dmn)

lnr <- linearity(NSP~., data = cardiotography)
print(lnr)

ovr <- overlapping(NSP~., data = cardiotography)
print(ovr)

nei <- neighborhood(NSP~., data = cardiotography)
print(nei)

nwk <- network(NSP~., data = cardiotography)
print(nwk)

class_col <- "NSP"
cats <- sapply(cardiotography, class)
facs <- names(cardiotography)[cats == "factor" & names(cardiotography) != class_col]
nums <- names(cardiotography)[cats %in% c("numeric", "integer") & names(cardiotography) != class_col]

asso <- 0
if (length(facs) > 0) {
  for (f in facs) {
    if (assocstats(table(cardiotography[[f]], cardiotography[[class_col]]))$chisq_test[2,3] <= 0.05) {
      asso <- asso + 1
    }
  }
}

prnc <- 0
if (length(nums) > 0)
  if (length(nums) == 1) {
    prnc <- 1
  } else {
    prnc <- max(nScree(cor(scale(cardiotography[nums])))$components
                 , sum(eigen(cor(scale(cardiotography[nums])))$values >= 1.0))
  }

shape <- c(n_row = nrow(cardiotography)
           , n_col = ncol(cardiotography)
           , n_facs = length(facs)
           , n_num = length(nums)
           , asso = asso
           , prnc = prnc)

meas <- c(shape, bal, dmn, lnr, ovr, nei, nwk)

complexity_measures <- rbind(complexity_measures, cardiotography = meas)

# adult
# bal <- balance(income~., data = adult) # C2 is best used for multiclass problems
print(bal)

# dmn <- dimensionality(income~., data = adult)
print(dmn)

# lnr <- linearity(income~., data = adult)
print(lnr)

# ovr <- overlapping(income~., data = adult)
print(ovr)

# nei <- neighborhood(income~., data = adult)
print(nei)

# nwk <- network(income~., data = adult)
print(nwk)
