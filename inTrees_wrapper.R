library(inTrees)
inTrees_wrapper <- function(X_train
                           , y_train
                           , model
                           , ntree) {
  # extract rules
  ruleExec0 <- extractRules(RF2List(model)
                            , X_train
                            , ntree = ntree) # hard limit on 5000 
  ruleExec <- unique(ruleExec0)
  ruleMetric <- getRuleMetric(ruleExec
                              , X_train
                              , y_train)
  
  #build the simplified tree ensemble learner
  learner <- buildLearner(ruleMetric
                          , X_train
                          , y_train)
  
  return(learner)
}
print("Created inTrees_wrapper function")
