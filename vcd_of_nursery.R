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
head(nursery)

table(nursery$decision)
nursery <- nursery[!(nursery$decision %in% c("recommend", "very_recom")), ]
nursery$decision <- factor(as.character(nursery$decision))

nur_dec_heal <- xtabs(~decision+health, data = nursery)
nur_dec_heal
chisq.test(nur_dec_heal)
residuals(chisq.test(nur_dec_heal))
assocstats(nur_dec_heal)
assoc(nur_dec_heal, shade=TRUE)
sieve(nur_dec_heal, sievetype = "expected", labeling = labeling_values, value_type = "expected", shade = TRUE)
sieve(nur_dec_heal, labeling = labeling_values, shade = TRUE)
fluctile(nur_dec_heal)

mosaic(nur_dec_heal
       , labeling = labeling_values
       , rot_labels = c(top = -90)
       , shade = TRUE)

res.ca <- CA(nur_dec_heal, graph = FALSE)
# fviz_ca_row(res.ca, repel = TRUE)
# fviz_ca_col(res.ca) + theme_minimal()
fviz_ca_biplot(res.ca, repel = TRUE)

nur_dec_heal_hasn <- xtabs(~decision+health+has_nurs, data = nursery)
mosaic(nur_dec_heal_hasn, shade = TRUE)

pairs(nur_dec_heal_hasn, gp=shading_Friendly2)

# correspondence analysis
# see also http://www.sthda.com/english/articles/22-principal-component-methods-videos/67-ca-in-r-using-factominer-quick-scripts-and-videos/