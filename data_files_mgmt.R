# setup
if (grepl("linux", Sys.getenv()[['R_LIBS_USER']])) {
  pathsep <- "/"
} else {
  pathsep <- "\\"
}

# data settings
# data_dir <- "C:\\Users\\id126493\\Documents\\GitHub\\explain_te\\CHIRPS\\datafiles\\"
# data_dir <- "C:\\Users\\Crutt\\Documents\\GitHub\\explain_te\\CHIRPS\\datafiles\\"
data_dir <- "~/Documents/github/explain_te/CHIRPS/datafiles/"
# datafilesdir <- "c:\\Dev\\Study\\Python\\explain_te\\forest_surveyor\\datafiles\\"

# project_dir <- "V:\\whiteboxing\\tests\\"
# project_dir <- "V:\\whiteboxing\\"
# project_dir <- "C:\\Users\\Crutt\\Documents\\whiteboxing\\tests\\"
project_dir <- "/datadisk/whiteboxing/2020Ada2_copy/"

datasetnames <- c("adult"
                  , "bankmark"
                  , "car"
                  , "cardio"
                  , "credit"
                  , "german"
                  , "lending_tiny_samp"
                  , "nursery"
                  , "rcdv")

n_classes <- c(2
               , 2
               , 2
               , 3
               , 2
               , 2
               , 2
               , 4
               , 2
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

difficulty <- c("large"
                , "small"
                , "small"
                , "small"
                , "small"
                , "small"
                , "small"
                , "small"
                ,"large")

positive_classes <- c(">50K", "yes", "acc", NA, "plus", "good", "Fully Paid", NA, "Y")
negative_classes <- c("<=50K", "no", "unacc", NA, "minus", "bad", "Charged Off", NA, "N")

datasets_master <- data.frame(class_cols, n_classes, positive_classes, negative_classes, difficulty, stringsAsFactors = FALSE)



rownames(datasets_master) <- datasetnames

get_data_files <- function(x = NA) {
  if (is.na(x)) x <- 1:nrow(datasets_master)
  paste0(rownames(datasets_master[x, ]), ".csv.gz")
}

data_files <- get_data_files()

get_measure_stem <- function(x) {
  paste("mean", sub("wx", "exc. ", sub(".tt.|.tr.", "", measure)))
}

train_set_size <- integer(nrow(datasets_master))
test_set_size <- integer(nrow(datasets_master))
names(train_set_size) <- rownames(datasets_master)
names(test_set_size) <- rownames(datasets_master)

get_best_params <- function() {
  fromJSON(readLines(file(paste0(
  project_dir
  , datasetnames[i]
  , pathsep
  , model
  , "_best_params_rnst_"
  , random_states[r]
  , ".json"))))
}
