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

balance(decision~., data = nursery)
# C2 is best used for multiclass problems
balance(decision~., data = nursery)["C2"]

# we know nursery is very unbalanced
table(nursery$decision)

# what's the effect of removing the minority classes?
nursery_simple <- nursery[!(nursery$decision %in% c("recommend", "very_recom")), ]
nursery_simple$decision <- factor(as.character(nursery_simple$decision))
table(nursery_simple$decision)


nursery_simple2 <- nursery[nursery$decision != "recommend", ]
nursery_simple2$decision <- factor(as.character(nursery_simple2$decision))
table(nursery_simple2$decision)


# the balance ratio is much closer to 1
balance(decision~., data = nursery_simple)["C2"]


dimensionality(decision~., data = nursery)
dimensionality(decision~., data = nursery_simple) # T2 and T3 are slightly less
# try dimensionality on one hot encoded data
nursery_mm <- model.matrix(decision~., data = nursery)
nursery_simple_mm <- model.matrix(decision~., data = nursery_simple)
dimensionality(nursery$decision~., data = as.data.frame(nursery_mm))
dimensionality(nursery_simple$decision~., data = as.data.frame(nursery_simple_mm)) # T2 and T3 are slightly less
# no effect? must already be converting frame to mm!


# takes time to run
linearity(decision~., data = nursery_simple2) # one class has only one member
linearity(decision~., data = nursery_simple)
# removing the remaining minority class returns higher results
# an increment of 0.01-0.02 on each measure

# these fail because the calculation is too large
# neighborhood(decision~., data = nursery_simple2, measures = "N1") # one class has only one member
# neighborhood(decision~., data = nursery_simple, measures = "N1") # T2 and T3 are slightly less

network(decision~., data = nursery_simple2) # one class has only one member
network(decision~., data = nursery_simple) # T2 and T3 are slightly less

overlapping(decision~., data = nursery_simple2) # one class has only one member
overlapping(decision~., data = nursery_simple) # T2 and T3 are slightly less

# F1v is much larger for simp2
# F2 is smaller
# F3 is larger
# F4 is larger

# doesn't work well for this data set 
# because of problems in the individual measurements
# complexity(decision~., data = nursery) 


sup <- round((rpois(10, 2) + 1) * 100 + rnorm(10, 0, 25))
len <- rpois(10, 2) + 1
scr <- function(s, l, a) {
  log(s) * (l - a)/l
}
a_vec <- seq(0, 1, 0.1)
scrs <- t(sapply(a_vec, function(x) {mapply(FUN = scr, sup, len, a = x)}))

matplot(scrs, pch = as.character(len)
        , type = "o"
        , col = rainbow(ncol(scrs))
        , xaxt = "n")
axis(1, at=seq_along(a_vec), labels = a_vec)
abline(v = which(a_vec == 0))

