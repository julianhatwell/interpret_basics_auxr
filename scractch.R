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


# scoring function
precis <- c(0.5507246376811594,
            0.5518018018018018,
            0.9587628865979382,
            0.9587628865979382,
            0.9608938547486033,
            0.9849624060150376,
            0.9846153846153847,
            1.0)
cover <- c(1,
           0.9192546583850931,
           0.40165631469979296,
           0.40165631469979296,
           0.37060041407867494,
           0.2753623188405797,
           0.2691511387163561,
           0.0)
card <- 0:7


acc2 <- c(0.4389233954451346,
          0.2546583850931677,
          0.6376811594202898,
          0.6231884057971014,
          0.6231884057971014,
          0.6107660455486542,
          0.5817805383022774,
          0.5817805383022774,
          0.0)
cover2 <- c(1.0,
            0.4858490566037736,
            0.4009433962264151,
            0.19811320754716982,
            0.19811320754716982,
            0.13679245283018868,
            0.05188679245283019,
            0.05188679245283019,
            0.0)
precis2 <- c(0.4389233954451346,
             0.2909604519774011,
             0.6390977443609023,
             0.7777777777777778,
             0.7777777777777778,
             0.8529411764705882,
             0.9166666666666666,
             0.9166666666666666,
             1.0)

card2 <- 0:8

alpha <- 0.5

lf <- function(x) {
  log2(x + 1)
}

sf <- function(prec, cove, crd, alph) {
  ((alph * lf(prec) + (1 - alph) * lf(cove)) * crd / (1 + crd^2)) 
}

# score <- (alpha * precis + (1 - alpha) * cover) 
score <- mapply(sf, precis2, cover2, card2, alpha) # * card2 / 2
score2 <- ((lf(precis2) + lf(cover2)) ) * card2 / (1 + card2^2)
score3 <- lf(acc2 * card2 / (1 + card2^2))
matplot(matrix(c(precis2, cover2, acc2, score, score2, score3), ncol = 6), type = "l")


# analysis_groups$mean <- apply(analysis_values, 2, mean)
# analysis_groups$rank_mean <- colMeans(t(apply(-analysis_values, 1, rank)))
# analysis_groups$rank_sum <- colSums(t(apply(-analysis_values, 1, rank)))
# 
# bins_plot <- ggplot(
#   data = analysis_groups
#   , aes(y = mean
#         , x = id
#         , shape = alpha
#         , colour = func
#         , size = support)) +
#   geom_point() + 
#   # scale_color_grey() +
#   scale_size_discrete(range = c(2, 4)) +
#   theme_bw() +
#   facet_grid(weights~bins)

analysis_groups_8 <- analysis_groups %>% filter(bins == 8)


analysis_bins_check <- inner_join(analysis_groups_4, analysis_groups_8
                                  , by = c("support", "alpha", "func", "weights")
                                  , suffix = c("_4", "_8"))

rank_total <- tail(cumsum(1:nrow(analysis_values)), 1)
wilcos <- list()
for (i in 1:nrow(analysis_bins_check)) {
  wilcos[[i]] <- wilcox.test(analysis_values[, analysis_bins_check$id_4[i]]
                             , analysis_values[, analysis_bins_check$id_8[i]]
                             , paired = TRUE, conf.int = TRUE)
}

bins_wilcos <- logical(length(wilcos))
for (i in 1:length(wilcos)) { bins_wilcos[i] <- wilcos[[i]]$p.value < 0.05 / length(wilcos) }
mean(bins_wilcos)

eff_sizes <- numeric(0)  
for(wb in which(bins_wilcos)) {eff_sizes <- c(eff_sizes, wilcos[[wb]]$statistic/rank_total)}

wilcos <- list()
for (i in 1:nrow(analysis_weights_check)) {
  wilcos[[i]] <- wilcox.test(analysis_values[, analysis_bins_check$id_4[i]]
                             , analysis_values[, analysis_bins_check$id_8[i]]
                             , paired = TRUE, conf.int = TRUE)
}

bins_wilcos <- logical(length(wilcos))
for (i in 1:length(wilcos)) { bins_wilcos[i] <- wilcos[[i]]$p.value < 0.05 / length(wilcos) }
mean(bins_wilcos)

# bins_wilcox <- wilcox.test(analysis_bins_check$mean_4
#                           , analysis_bins_check$mean_8
#                           , paired = TRUE)
bins_wilcox <- wilcox.test(analysis_values[, analysis_bins_check$id_4]
                           , analysis_values[, analysis_bins_check$id_8]
                           , paired = TRUE
                           , conf.int = TRUE
                           , conf.level = 0.99)

bins_effect_size <- mean(analysis_values[, analysis_bins_check$id_4] -
                           analysis_values[, analysis_bins_check$id_8]) /
  sd(analysis_values[, analysis_bins_check$id_4] -
       analysis_values[, analysis_bins_check$id_8])

# we know we only need to keep bins_4
analysis_groups_chisq <- analysis_groups_4 %>% filter(weights == "chisq")
analysis_groups_nothing <- analysis_groups_4 %>% filter(weights == "nothing")

analysis_weights_check <- inner_join(analysis_groups_chisq, analysis_groups_nothing
                                     , by = c("support", "alpha", "func")
                                     , suffix = c("_chisq", "_nothing")) %>%
  select(-weights_chisq, -weights_nothing)

# weights_wilcox <- wilcox.test(analysis_weights_check$mean_chisq
#                              , analysis_weights_check$mean_nothing
#                              , paired = TRUE)

weights_wilcox <- wilcox.test(analysis_values[, analysis_weights_check$id_chisq]
                              , analysis_values[, analysis_weights_check$id_nothing]
                              , paired = TRUE
                              , conf.int = TRUE
                              , conf.level = 0.99)

weights_effect_size <- mean(analysis_values[, analysis_weights_check$id_chisq] -
                              analysis_values[, analysis_weights_check$id_nothing]) /
  sd(analysis_values[, analysis_weights_check$id_chisq] -
       analysis_values[, analysis_weights_check$id_nothing])

analysis_groups_chisq$mean <- apply(analysis_included_values, 2, mean)
analysis_groups_chisq$rank_mean <- colMeans(t(apply(-analysis_included_values, 1, rank)))

best_blocks_plot <- ggplot(
  data = analysis_groups_chisq
  , aes(y = mean
        , x = id
        , shape = alpha
        , colour = func
        , size = support)) +
  geom_point() + 
  scale_size_discrete(range = c(2, 4)) +
  theme_bw()
