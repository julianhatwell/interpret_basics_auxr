library(dplyr)
library(xlsx)

# utility code
onehot_encode <- function(dframe) {
  
  input_names <- names(dframe)
  output_names <- list()
  numerics <- list()
  for (name in input_names) {
    if (class(dframe[[name]]) %in% c("numeric", "integer") ) {
      numerics[[name]] <- dframe[[name]]
      dframe[[name]] <- NULL
    } else {
      output_names[[name]] <- paste(name, unique(dframe[, name]), sep = "_")  
    }
  }
  dframe <- as.data.frame(
    lapply(dframe, function(x) {
      
      vals <- unique(x)
      one_hots <- matrix(nrow = nrow(dframe)
                         , ncol = length(vals))
      
      for (i in seq_along(vals)) {
        one_hots[, i] <- x == vals[i]
      }
      one_hots
    })
  )
  names(dframe) <- unlist(output_names)
  # add any untreated numerics back to the end of the new frame
  if(length(numerics) > 0) { # there are some numeric columns
    if (length(dframe) > 0) { # there are also some non-numeric columns
      dframe <- cbind(dframe, numerics)  
    } else {
      dframe <- as.data.frame(numerics) # there were no non-numeric columns
    }
  }
  dframe
}

homedir <- "c:\\Dev\\Study\\Python\\interpret_basics2\\"

print("Common setup for nursery data")
nursery <- read.csv(gzfile(paste0(homedir, "nursery_pickles\\nursery.csv.gz")))

print("nursery: basic dataset")

nursery_enc <- nursery %>% select(-decision)
nursery_enc <- onehot_encode(nursery_enc)
# add the unencoded class col back
nursery_enc$decision <- nursery$decision
print("nursery_enc: one hot encoded dataset")

# german
print("Common setup for german data")

german <- read.csv(gzfile(paste0(homedir, "german_pickles\\german.csv.gz")))
print("german: basic dataset")

german_enc <- german %>% select(-rating)
german_enc <- onehot_encode(german_enc)
german_enc$rating <- german$rating
print("german_enc: one hot encoded dataset")

# credit
print("Common setup for credit data")

credit <- read.csv(gzfile(paste0(homedir, "credit_pickles\\credit.csv.gz")))
print("credit: basic dataset")

credit_enc <- credit %>% select(-A16)
credit_enc <- onehot_encode(credit_enc)
credit_enc$A16 <- credit$A16
print("credit_enc: one hot encoded dataset")

# adult
print("Common setup for adult data")

adult <- read.csv(gzfile(paste0(homedir, "adult_pickles\\adult.csv.gz")))
print("adult: basic dataset")

adult_enc <- adult %>% select(-income)
adult_enc <- onehot_encode(adult_enc)
adult_enc$income <- adult$income
print("adult_enc: one hot encoded dataset")

# car
print("Common setup for car data")

car <- read.csv(gzfile(paste0(homedir, "car_pickles\\car.csv.gz")))
print("car: basic dataset")

car_enc <- car %>% select(-acceptability)
car_enc <- onehot_encode(car_enc)
car_enc$income <- car$income
print("car_enc: one hot encoded dataset")

# cardiotrography
print("Common setup for cardiotography data")

cardiotography <- read.csv(gzfile(paste0(homedir, "cardiotography_pickles\\cardiotography.csv.gz")))
print("cardiotography: basic dataset")

cardiotography_enc <- cardiotography %>% select(-NSP) # encode without the class column
cardiotography_enc <- onehot_encode(cardiotography_enc)
cardiotography_enc$NSP <- factor(cardiotography$NSP) # add back the class column
print("cardiotography_enc: one hot encoded dataset")

# cardiotrography
print("Common setup for lending data")
lending <- read.csv(gzfile(paste0(homedir, "lending_pickles\\lend_samp.csv.gz")))
print("lending: basic dataset")

lending_enc <- lending %>% select(-loan_status)
lending_enc <- onehot_encode(lending_enc)
lending_enc$loan_status <- lending$loan_status
print("lending_enc: one hot encoded dataset")

# rcdv
print("Common setup for rcdv data")

rcdv <- read.csv(gzfile(paste0(homedir, "rcdv_pickles\\rcdv.csv.gz")))
print("rcdv: basic dataset")

rcdv_enc <- rcdv %>% select(-recid)
rcdv_enc <- onehot_encode(rcdv_enc)
rcdv_enc$recid <- rcdv$recid
print("rcdv_enc: one hot encoded dataset")
