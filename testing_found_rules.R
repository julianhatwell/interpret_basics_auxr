# investigate first frequent patterns

library(vcd)
library(FactoMineR)
library(factoextra)
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

# utility code
onehot_encode <- function(dframe) {
  
  input_names <- names(dframe)
  output_names <- list()
  for (name in input_names) {
    output_names[[name]] <- paste(name, unique(dframe[, name]), sep = "_")
  }
  
  dframe <- as.data.frame(
    sapply(dframe, function(x) {
      
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
  dframe  
}

nursery_onehot <- onehot_encode(nursery)
# add the unencoded class col back
nursery_onehot$decision <- as.character(nursery$decision)

# test instance 0 - frequent patterns
# apriori(transactions = paths_0, support = 0.2, max_itemset_size = 5)
# sorted_0 = sort_fp(freq_patt_0, by_score=True, alpha=0.4)
# freq patts scored with alpha = 0.4
# [(('health_not_recom',), 687),
#   (('health_priority', 'health_recommended'), 305),
#   (('has_nurs_very_crit',), 243),
#   (('has_nurs_proper',), 240),
#   (('parents_great_pret',), 206),
#   (('has_nurs_less_proper',), 200)]

# looks like health not recom is the most important by far
xtabs(~decision+health_not_recom+health_priority+health_recommended
      , data = nursery_onehot)

cotabplot(xtabs(~decision+health_not_recom+health_priority+health_recommended
                , data = nursery_onehot)
          , panel = cotab_loddsratio
)
cotabplot(xtabs(~decision+health_not_recom+health_priority+health_recommended
                , data = nursery_onehot)
          , panel = cotab_loddsratio
          , cond = c(4, 3)
)



structable(decision~health_not_recom+health_priority+health_recommended
           , data = nursery_onehot)

# go to full screen
ftable(nursery_onehot[c("health_not_recom", "health_priority", "health_recommended", "decision")])

margin.table(table(nursery_onehot[c("health_not_recom", "health_priority", "health_recommended", "decision")]), margin = c(2, 3, 4))
margin.table(table(nursery_onehot[c("health_not_recom", "health_priority", "health_recommended", "decision")]), margin = c(1, 4))
prop.table(table(nursery_onehot[c("health_not_recom", "decision")]), margin = 1)


mosaic(xtabs(~decision+health_not_recom+health_priority+health_recommended
      , data = nursery_onehot), shade = TRUE, rot_labels = c(0, 0, 0, 0))



with(nursery_onehot,
     mosaic(table(decision
                  , health_not_recom)
            , shade = TRUE
            , rot_labels = c(0, 0))
)


with(nursery_onehot,
     mosaic(table(decision
           , health_not_recom
           , health_priority
           , health_recommended)
           , shade = TRUE, rot_labels = c(0, 0, 0, 0))
     )

mosaic(structable(~decision+health_not_recom+health_priority+health_recommended
                 , data = nursery_onehot)
      , shade = TRUE, rot_labels = c(0, 0, 0, 0))


mosaic(table(nursery[c("decision", "health", "has_nurs")]), shade = TRUE)
cotabplot(table(nursery[c("decision", "health", "has_nurs")])
          , shade = TRUE
          , cond = 2)
assoc(table(nursery[c("decision", "health")]), shade = TRUE)
assoc(table(nursery[c("decision", "has_nurs")]), shade = TRUE)          


cotabplot(xtabs(~decision+health_not_recom+has_nurs_very_crit
                , data = nursery_onehot)
          , panel = cotab_loddsratio
)
cotabplot(xtabs(~decision+health_not_recom+has_nurs_very_crit
                , data = nursery_onehot)
          , panel = cotab_mosaic
          , cond = 1)

prop.table(table(nursery[c("health", "decision")]), margin = 1)
prop.table(table(nursery[c("has_nurs", "decision")]), margin = 1)
prop.table(table(nursery[c("parents", "decision")]), margin = 1)

prop.table(table(nursery[c("health", "has_nurs", "decision")]), c(2,3))


with(nursery_onehot, 
table(nursery_onehot[has_nurs_less_proper < 0.5 &
                       has_nurs_proper < 0.5 &
                       has_nurs_very_crit < 0.5 &
                       health_priority > 0.5, "decision"]
))

with(nursery_onehot, 
     table(nursery_onehot[has_nurs_less_proper < 0.5 &
                            has_nurs_proper < 0.5 &
                            has_nurs_very_crit < 0.5 &
                            parents_great_pret > 0.5, "decision"]
     ))

with(nursery_onehot, 
     round(prop.table(table(nursery_onehot[has_nurs_less_proper < 0.5 &
                            has_nurs_critical < 0.5 &
                            has_nurs_very_crit < 0.5
                          , "decision"])), 4)
     )

with(nursery_onehot, 
     round(prop.table(table(nursery_onehot[has_nurs_less_proper < 0.5 &
                                             has_nurs_proper < 0.5 &
                                             parents_usual > 0.5
                                           , "decision"])), 4)
)

with(nursery_onehot, 
     round(prop.table(table(nursery_onehot[has_nurs_less_proper < 0.5 &
                                             has_nurs_critical < 0.5 &
                                             has_nurs_very_crit < 0.5 &
                                             has_nurs_proper < 0.5 &
                                             parents_usual > 0.5
                                           , "decision"])), 4)
)

with(nursery_onehot, 
     round(prop.table(table(nursery_onehot[has_nurs_less_proper < 0.5 &
                                             has_nurs_critical < 0.5 &
                                             has_nurs_very_crit < 0.5 &
                                             has_nurs_proper < 0.5 &
                                             parents_usual > 0.5 &
                                             health_not_recom < 0.5
                                           , "decision"])), 4)
)

with(nursery_onehot, 
     round(prop.table(table(nursery_onehot[has_nurs_improper > 0.5 &
                                             parents_usual > 0.5 &
                                             health_not_recom < 0.5
                                           , "decision"])), 4)
)
with(nursery_onehot, 
     round(prop.table(table(nursery_onehot[has_nurs_improper > 0.5 &
                                             parents_usual > 0.5 &
                                             health_priority > 0.5
                                           , "decision"])), 4)
)


with(nursery_onehot, 
     round(prop.table(table(nursery_onehot[children_1 > 0.5 &
                                             has_nurs_very_crit < 0.5 &
                                             parents_great_pret < 0.5 &
                                             has_nurs_less_proper < 0.5 &
                                             form_foster < 0.5
                                           , "decision"])), 4)
)

# binary complement
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








