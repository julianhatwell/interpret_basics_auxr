library(ECoL)
nursery <- read.csv("c:\\Dev\\Study\\Python\\interpret_basics\\nursery.csv"
                    , col.names = c('parents'
                                    , 'has_nurs'
                                    , 'form'
                                    , 'children'
                                    , 'housing'
                                    , 'finance'
                                    , 'social'
                                    , 'health'
                                    , 'decision')
)

# C2 is best used for multiclass problems
# ratio, should be close to 1
balance(decision~., data = nursery)["C2"]
# we know nursery is very unbalanced
table(nursery$decision)

# this is the same for all variations
dimensionality(decision~., data = nursery)

class_names <- as.character(unique(nursery$decision))

for (cn in class_names) {
  nursery_temp <- nursery
  nursery_temp[[cn]] <- FALSE
  nursery_temp[[cn]][as.character(nursery$decision) == cn] <- TRUE
  nursery_temp$decision <- NULL
  assign(paste0("nursery_", cn)
         , nursery_temp) 
         
}

# C2 is best used for multiclass problems
# ratio, should be close to 1
balance(not_recom~., data = nursery_not_recom)["C2"]
balance(priority~., data = nursery_priority)["C2"]
balance(recommend~., data = nursery_recommend)["C2"]
balance(spec_prior~., data = nursery_spec_prior)["C2"]
balance(very_recom~., data = nursery_very_recom)["C2"]

linearity(not_recom~., data = nursery_not_recom)
linearity(priority~., data = nursery_priority)
linearity(recommend~., data = nursery_recommend)
linearity(spec_prior~., data = nursery_spec_prior)
linearity(very_recom~., data = nursery_very_recom)

overlapping(not_recom~., data = nursery_not_recom)
overlapping(priority~., data = nursery_priority)
overlapping(recommend~., data = nursery_recommend)
overlapping(spec_prior~., data = nursery_spec_prior)
overlapping(very_recom~., data = nursery_very_recom)

target_url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/car/car.data'

car <- read.csv(target_url
                , col.names = c('buying'
                               , 'maint'
                               , 'doors'
                               , 'persons'
                               , 'lug_boot'
                               , 'safety'
                               , 'acceptability')
                , stringsAsFactors = FALSE
                
                    
)

# C2 is best used for multiclass problems
# ratio, should be close to 1
balance(acceptability~., data = car)["C2"]
# we know nursery is very unbalanced
table(car$acceptability)

# this is the same for all variations
dimensionality(acceptability~., data = car)

class_names <- as.character(unique(car$acceptability))

for (cn in class_names) {
  car_temp <- car
  car_temp[[cn]] <- FALSE
  car_temp[[cn]][as.character(car$acceptability) == cn] <- TRUE
  car_temp$acceptability <- NULL
  assign(paste0("car_", cn)
         , car_temp) 
  
}

# C2 is best used for multiclass problems
# ratio, should be close to 1
balance(acc~., data = car_acc)["C2"]
balance(good~., data = car_good)["C2"]
balance(unacc~., data = car_unacc)["C2"]
balance(vgood~., data = car_vgood)["C2"]

linearity(acc~., data = car_acc)
linearity(good~., data = car_good)
linearity(unacc~., data = car_unacc)
linearity(vgood~., data = car_vgood)

overlapping(acc~., data = car_acc)
overlapping(good~., data = car_good)
overlapping(unacc~., data = car_unacc)
overlapping(vgood~., data = car_vgood)

neighborhood(acc~., data = car_acc)
neighborhood(good~., data = car_good)
neighborhood(unacc~., data = car_unacc)
neighborhood(vgood~., data = car_vgood)
