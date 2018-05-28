library(dplyr)
source("dirs.R")

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

# accident
print("Common setup for accident data")

accident <- read.csv(gzfile(paste0(datafilesdir, "accident_small_samp.csv.gz")))
print("accident: basic dataset")

accident_enc <- accident %>% dplyr::select(-Accident_Severity)
accident_enc <- onehot_encode(accident_enc)
accident_enc$Accident_Severity <- accident$Accident_Severity
print("accident_enc: one hot encoded dataset")

# adult
print("Common setup for adult data")

adult <- read.csv(gzfile(paste0(datafilesdir, "adult_small_samp.csv.gz")))
print("adult: basic dataset")

adult_enc <- adult %>% dplyr::select(-income)
adult_enc <- onehot_encode(adult_enc)
adult_enc$income <- adult$income
print("adult_enc: one hot encoded dataset")

# lending
print("Common setup for lending data")
lending <- read.csv(gzfile(paste0(datafilesdir, "lending_tiny_samp.csv.gz")))
print("lending: basic dataset")

lending_enc <- lending %>% dplyr::select(-loan_status)
lending_enc <- onehot_encode(lending_enc)
lending_enc$loan_status <- lending$loan_status
print("lending_enc: one hot encoded dataset")


# rcdv
print("Common setup for rcdv data")

rcdv <- read.csv(gzfile(paste0(datafilesdir, "rcdv_samp.csv.gz")))
print("rcdv: basic dataset")

rcdv_enc <- rcdv %>% dplyr::select(-recid)
rcdv_enc <- onehot_encode(rcdv_enc)
rcdv_enc$recid <- rcdv$recid
print("rcdv_enc: one hot encoded dataset")