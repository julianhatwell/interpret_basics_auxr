source("nursery_setup.R")
library(vcd)
library(FactoMineR)
library(factoextra)
library(extracat)

# test instance 5
('parents_usual', False, 0.5),
('form_complete', False, 0.5),
('health_recommended', False, 0.5),
('housing_less_conv', False, 0.5),
('has_nurs_improper', False, 0.5),
('children_1', False, 0.5)

with(nursery_enc, 
     table(nursery_enc[parents_usual > 0.5 &
                         form_complete > 0.5 &
                         health_recommended > 0.5 &
                         housing_less_conv > 0.5 &
                         has_nurs_improper > 0.5 &
                         children_1 > 0.5, "decision"]
     ))

struc <- structable(decision~
        parents_usual+
        form_complete+
        health_recommended+
        housing_less_conv+
        has_nurs_improper+
        children_1
      , data = nursery_enc)

# how many samples remain out of total (training set, test set), when isolated by rule?
# feature depth and path length to use in scoring?
# more than one compound rule

# stopping when rule covers/discriminates
# shortest rule when..... how to select? depths, weightings, GA and optimization
# how do rules relate to depths?
# analysing how RF tuning parameters affect efficiency.


# work out the probabilities that feature appears in a tree
# and how gini index affects it (ranking of features)
# permutation effect








