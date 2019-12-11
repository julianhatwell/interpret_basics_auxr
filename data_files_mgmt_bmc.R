# data settings
data_dir <- "C:\\Users\\id126493\\Documents\\GitHub\\explain_te\\CHIRPS\\datafiles\\"

project_dir <- "V:\\whiteboxing\\"

# setup
if (grepl("linux", Sys.getenv()[['R_LIBS_USER']])) {
  pathsep <- "/"
} else {
  pathsep <- "\\"
}

datasetnames <- c("breast"
                  , "cardio"
                  , "diaretino"
                  , "heart"
                  , "mhtech14"
                  , "mh2tech16"
                  , "readmit"
                  , "thyroid"
                  , "usoc2")

n_classes <- c(2
               , 3
               , 2
               , 2
               , 2
               , 2
               , 2
               , 2
               , 3
               )

class_cols <- c(
  "mb"
  , "NSP"
  , "dr"
  , "HDisease"
  , "treatment"
  , "mh2"
  , "readmit"
  , "diagnosis"
  , "mh"
)

datasets_master <- data.frame(class_cols, n_classes, stringsAsFactors = FALSE)

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
