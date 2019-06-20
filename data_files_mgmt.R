# data settings
# data_dir <- "C:\\Users\\id126493\\Documents\\GitHub\\explain_te\\CHIRPS\\datafiles\\"
data_dir <- "C:\\Users\\Crutt\\Documents\\GitHub\\explain_te\\CHIRPS\\datafiles\\"
# data_dir <- "~/github/explain_te/CHIRPS/datafiles/"
# datafilesdir <- "c:\\Dev\\Study\\Python\\explain_te\\forest_surveyor\\datafiles\\"

# project_dir <- "V:\\whiteboxing\\"
project_dir <- "C:\\Users\\Crutt\\Documents\\whiteboxing\\tests\\"
#project_dir <- "/home/julian/whiteboxing/"

n_classes <- c(adult_small_samp = 2
               , bankmark_samp = 2
               , car = 2
               , cardio = 3
               , credit = 2
               , german = 2
               , lending_tiny_samp = 2
               , nursery_samp = 4
               , rcdv_samp = 2
               )

class_cols <- c(
  "income"
  , "y"
  , "acceptability"
  , "NSP"
  , "A16"
  , "rating"
  , "loan_status"
  , "decision"
  , "recid"
)

data_files <- c(
  "adult_small_samp.csv.gz"
  , "bankmark_samp.csv.gz"
  , "car.csv.gz"
  , "cardio.csv.gz"
  , "credit.csv.gz"
  , "german.csv.gz"
  , "lending_tiny_samp.csv.gz"
  , "nursery_samp.csv.gz"
  , "rcdv_samp.csv.gz"
)

get_datasetnames <- function(x) {
  sub(".csv.gz", "", x)
}

get_datasetname_stems <- function(x) {
  sub("_samp", "", sub("_small|_tiny", "", x))
}

get_measure_stem <- function(x) {
  paste("mean", sub("wx", "exc. ", sub(".tt.|.tr.", "", measure)))
}

datasetnames <- sapply(data_files, get_datasetnames)
resfilesdirs <- paste0(project_dir, datasetnames, "\\")
resfilesdirs <- paste0(project_dir, datasetnames, "/")

