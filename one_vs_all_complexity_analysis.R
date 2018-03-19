library(ECoL)
source("general_setup.R")
# now check the complexity measures

# C2 is best used for multiclass problems
# ratio, should be close to 1
balance(decision~., data = nursery_enc)["C2"]
# we know nursery_enc is very unbalanced
table(nursery_enc$decision)

# this is the same for all variations
dimensionality(decision~., data = nursery_enc)

# one versus all analysis

class_names <- as.character(unique(nursery_enc$decision))
for (cn in class_names) {
  nursery_enc_temp <- nursery_enc
  nursery_enc_temp[[cn]] <- FALSE
  nursery_enc_temp[[cn]][as.character(nursery_enc$decision) == cn] <- TRUE
  nursery_enc_temp$decision <- NULL
  assign(paste0("nursery_enc_", cn)
         , nursery_enc_temp) 
         
}

# C2 is best used for multiclass problems
# ratio, should be close to 1
balance(not_recom~., data = nursery_enc_not_recom)["C2"]
balance(priority~., data = nursery_enc_priority)["C2"]
balance(spec_prior~., data = nursery_enc_spec_prior)["C2"]
balance(very_recom~., data = nursery_enc_very_recom)["C2"]

linearity(decision~., data = nursery_enc)
linearity(not_recom~., data = nursery_enc_not_recom)
linearity(priority~., data = nursery_enc_priority)
linearity(spec_prior~., data = nursery_enc_spec_prior)
linearity(very_recom~., data = nursery_enc_very_recom)

overlapping(not_recom~., data = nursery_enc_not_recom)
overlapping(priority~., data = nursery_enc_priority)
overlapping(spec_prior~., data = nursery_enc_spec_prior)
overlapping(very_recom~., data = nursery_enc_very_recom)


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
# we know nursery_enc is very unbalanced
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
