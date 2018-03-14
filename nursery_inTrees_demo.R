source("nursery_setup.R")
source("inTrees_wrapper.R")
library(randomForest)
library(caret)

set.seed(23)
training_ids <- createDataPartition(y = nursery$decision, p = 0.8, list = FALSE)

training_set <- nursery[training_ids,]
validation_set <- nursery[-training_ids,]
training_labels <- nursery$decision[training_ids]
validation_labels <- nursery$decision[-training_ids]
n_tree <- 1000

rf <- randomForest(training_labels~.
                   , data = training_set
                   , n_tree = n_tree)

learner <- inTrees_wrapper(training_set
                           , training_labels
                           , model = rf
                           , ntree = n_tree)

# present the rules with a more readable format
readableLearner <- presentRules(learner
                                , colnames(training_set))
readableLearner

# predict and and assess
pred_it <- applyLearner(learner, validation_set)
cm_it <- confusionMatrix(table(pred_it, true = validation_labels))

pred_rf <- predict(rf, newdata = validation_set)
cm_rf <- confusionMatrix(table(pred_rf, true = validation_labels))

cm_it$table
cm_rf$table

readableLearner[1,]
readableLearner[2,]
readableLearner[3,]