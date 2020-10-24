get_means_for_plotting <- function(meas, div = 1, select_algos = algorithms) {
  
  means <- as.data.frame(get_mean_of(meas) / div)[, select_algos]
  means$dataset <- rownames(means)
  means <- gather(means, algorithm, m, -dataset)
  
  st_errs <- as.data.frame((get_sd_of(meas) / test_set_size_sqrt) / div)[, select_algos]
  st_errs$dataset <- rownames(st_errs)
  st_errs <- gather(st_errs, algorithm, st_err, -dataset)
  
  means$lwr <- means$m - st_errs$st_err
  means$upr <- means$m + st_errs$st_err
  return(means)
}

get_points_for_plotting <- function(func_table, values_to_name, rows_to_name = "dataset") {
  vtn <- enquo(values_to_name)
  out <- as.data.frame(func_table)
  out[[rows_to_name]] <- dimnames(get_above.75_of(meas))[[1]]
  out <- gather(out, algorithm, !!vtn, -dataset)
  return(out)
}

# plot themes
source("KTheme.R")
reds <- k.grad.red.rev(2)
ongs <- k.grad.orange.rev(2)
blus <- k.grad.blue.rev(2)
grns <- k.grad.green.rev(2)
grys <- k.grad.grey.rev(2)

myPal1 <- c(blus[2], reds[2], ongs[2], grns[2], grys[2])
myPal2 <- myPal1
myPal2 <- paste(myPal2, c("", rep("55", 4)), sep = "")
myAlph1 <- c(1, rep(0.7, 4))
myAlphFixed <- c(rep(1.0, 9), rep(0.5, 36))
myShap <- c(15, 11, 2, 5, 1)
names(myPal1) <- algorithms
names(myPal2) <- algorithms
names(myAlph1) <- algorithms
names(myShap) <- algorithms
names(myAlphFixed) <- rep(algorithms, each = 9)
colos1 <- scale_colour_manual(
  values = myPal1)
colos2 <- scale_colour_manual(
  values = myPal2) # with alpha built in
alpas <- scale_alpha_manual(values = myAlph1)
shaps <- scale_shape_manual(values = myShap)
sens_colos_silent <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue)
                                         , guide = FALSE)
sens_colos <- scale_colour_manual(values = c(k.purple, k.pink, k.brightblue))

get_ylabel <- function(y) {
  ylabel <- gsub("cc", "coverage of target class"
                 , gsub("wx", "exclusive "
                        , gsub("stability", "reliability"
                               , gsub("time", "time (sec)"
                                      , gsub("_", " "
                                             , gsub("rule.length", "rule length (antecedent cardinality)"
                                                    , gsub(".tt.", "", y)))))))
}

get_main_plot <- function(meas
                          , sc_y = scale_y_continuous()
                          , div = 1
                          , select_algos = algorithms) {
  myGeomPoint <- geom_point(size = 5)
  myGeomErrorBar <- geom_errorbar(width = 0.25)
  myGeomLine <- geom_line(size = 0.5, alpha = myAlphFixed[names(myAlphFixed) %in% select_algos])
  ggplot(data = get_means_for_plotting(meas, div, select_algos)
         , aes(y = m
               , ymin = lwr
               , ymax = upr
               , x = dataset
               , colour = algorithm
               , group = algorithm
               , alpha = algorithm
               , shape = algorithm)) +
    myGgTheme +
    myGeomLine +
    myGeomPoint +
    myGeomErrorBar +
    colos1 + 
    alpas +
    shaps +
    sc_y +
    ylab(get_ylabel(meas)) +
    xlab("data set") +
    theme(
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 26),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 26)
    )
}
get_main_plot("precision.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
get_main_plot("stability.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
get_main_plot("coverage.tt.", scale_y_continuous(limits = c(0.0, 1.00)))
get_main_plot("wxcoverage.tt.", scale_y_continuous(limits = c(0.0, 0.6)))
get_main_plot("rule.length")
get_main_plot("elapsed_time", scale_y_log10())

get_followup_plot <- function(
  summary_table
  , y_label
  #                         
  , sc_y = scale_y_continuous()
  , select_algos = algorithms) {
  myGeomPoint <- geom_point(size = 5)
  myGeomLine <- geom_line(size = 0.5, alpha = myAlphFixed[names(myAlphFixed) %in% select_algos])
  ggplot(data = get_points_for_plotting(summary_table
                                        , values_to_name = "place_holder") # get_points_for_plotting(summary_table, values_to_name = values_to_name)
         , aes(y = place_holder
               , x = dataset
               , colour = algorithm
               , group = algorithm
               , alpha = algorithm
               , shape = algorithm)) +
    myGgTheme +
    myGeomLine +
    myGeomPoint +
    colos1 + 
    alpas +
    shaps +
    sc_y +
    ylab(y_label) +
    xlab("data set") +
    theme(
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 26),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 26)
    )
}
get_followup_plot(get_above.75_of("stability.tt.")
                  , y_label = "Reliability Floor (threshold 0.75)")

get_followup_plot(get_equal.zero_of("rule.length")
                  , sc_y = scale_y_continuous(limits = c(0.0, 1.0))
                  , y_label = "Rule Length Floor (threshold 1)")


algo_alphas <- myAlphFixed[
  as.vector(t(matrix(colnames(get_mean_of(meas))
                     , byrow = TRUE
                     , nrow = nrow(get_mean_of(meas))
                     , ncol = ncol(get_mean_of(meas)))))[!is.na(as.vector(t(get_mean_of(meas))))]
  ]

get_main_boxplot <- function(meas
                             , sc_y = scale_y_continuous()) {
  ggplot(data = dplyr::select(comp_results, instance_id, dataset, !! enquo(meas), algorithm)
         , aes(y = !!enquo(meas)
               , colour = algorithm)) +
    geom_boxplot(outlier.alpha = 0.1
                 , alpha = algo_alphas
                 , position = position_dodge2(width = 2
                                              , padding = 0.2)
    ) +
    colos2 +
    facet_wrap(dataset~.) +
    scale_x_continuous(labels = NULL
                       , breaks = 0) +
    sc_y +
    myGgTheme +
    theme(legend.position = "bottom"
          , legend.title = element_blank()
          , strip.text = element_text(size = 8
                                      , margin = ggplot2::margin(1,0,1,0, "pt"))
    ) +
    ylab(get_ylabel(as_label(enquo(meas)))) +
    theme(
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 26),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 26),
      strip.text = element_text(size = 30)
    )
}

get_main_boxplot(precision.tt.)
get_main_boxplot(stability.tt.)
get_main_boxplot(coverage.tt.)
get_main_boxplot(wxcoverage.tt.
                 , sc_y = scale_y_continuous(limits = c(0.0, 0.7)))
get_main_boxplot(rule.length)
get_main_boxplot(elapsed_time, scale_y_log10())


tikz(file = "cc.tikz", width = 6.85, height = 2)
get_main_plot("cc.tt."
              , sc_y = scale_y_continuous(limits = c(0.0, 1.0))
              , div = test_set_size
              , select_algos = c("Anchors", algo_variant))
dev.off()

tikz(file = "prec.tikz", width = 6.85, height = 2)
get_main_plot("precision.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

tikz(file = "stab.tikz", width = 6.85, height = 2)
get_main_plot("stability.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

tikz(file = "cov.tikz", width = 6.85, height = 2)
get_main_plot("coverage.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

tikz(file = "xcov.tikz", width = 6.85, height = 2)
get_main_plot("wxcoverage.tt.", scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

tikz(file = "etime.tikz", width = 6.85, height = 2)
get_main_plot("elapsed_time", scale_y_log10())
dev.off()

tikz(file = "stabbox.tikz", width = 6.85, height = 3)
get_main_boxplot(stability.tt.)
dev.off()
tikz(file = "xcovbox.tikz", width = 6.85, height = 3)
get_main_boxplot(wxcoverage.tt.
                 , sc_y = scale_y_continuous(limits = c(0.0, 1.0)))
dev.off()

#########
get_lwrq_of(meas) / test_set_size
get_median_of(meas) / test_set_size
get_uprq_of(meas) / test_set_size

get_lwrq_of("stability.tt.")
get_median_of("stability.tt.")
get_uprq_of("stability.tt.")

get_meds_for_plotting <- function(meas, div = 1, select_algos = algorithms) {
  
  meds <- as.data.frame(get_median_of(meas)[, select_algos] / div)
  meds$dataset <- rownames(meds)
  meds <- gather(meds, algorithm, m, -dataset)
  lwrqs <- as.data.frame(get_lwrq_of(meas) / div)[, select_algos]
  lwrqs$dataset <- rownames(lwrqs)
  lwrqs <- gather(lwrqs, algorithm, lwr, -dataset)
  uprqs <- as.data.frame(get_uprq_of(meas) / div)[, select_algos]
  uprqs$dataset <- rownames(uprqs)
  uprqs <- gather(uprqs, algorithm, upr, -dataset)
  meds$lwr <- lwrqs$lwr
  meds$upr <- uprqs$upr
  return(meds)
  
}

get_med_plot <- function(meas
                         , sc_y = scale_y_continuous()
                         , div = 1
                         , select_algos = algorithms) {
  ylabel <- gsub("cc", "coverage of target class"
                 , gsub("wx", "exclusive "
                        , gsub("_", " "
                               , gsub(".tt.", "", meas))))
  myGeomPoint <- geom_point(size = 1.5)
  myGeomErrorBar <- geom_errorbar(width = 0.25)
  myGeomLine <- geom_line(size = 0.25, alpha = myAlphFixed[names(myAlphFixed) %in% select_algos])
  # y_label <- get_ylabel <- function(meas)
  ggplot(data = get_meds_for_plotting(meas, div, select_algos)
         , aes(y = m
               , ymin = lwr
               , ymax = upr
               , x = dataset
               , colour = algorithm
               , group = algorithm
               , alpha = algorithm
               , shape = algorithm)) +
    myGgTheme +
    myGeomLine +
    myGeomPoint +
    myGeomErrorBar +
    colos1 + 
    alpas +
    shaps +
    sc_y +
    # scale_x_discrete(labels = get_datasetname_stems(datasetnames)) +
    ylab(ylabel) +
    xlab("data set")
}
get_med_plot("stability.tt.", scale_y_continuous(limits = c(0.0, 1.0)))


get_med_plot("cc.tt."
             , sc_y = scale_y_continuous(limits = c(0.0, 1.0))
             , div = test_set_size
             , select_algos = c("Anchors", algo_variant))


meas <- "stability.tt."
for (ds in datasetnames) {
  qmeas <- quo(meas)
  cres <- comp_results %>% filter(dataset == ds) %>%
    dplyr::select(instance_id, dataset, !! qmeas, algorithm)
  cres <- spread(cres, algorithm, !! qmeas)
  mrs <- apply(apply(-cres[, algorithms], 1, rank), 1, mean)
  print(mrs)
  
  cres_frd <- friedman.test(as.matrix(cres))
  print(ds)
  print(cres_frd)
  
  chisqstat <- cres_frd$statistic
  df1 <- cres_frd$parameter
  N <- nrow(cres)
  fstat <- (chisqstat * (N - 1)) / (N * df1 - chisqstat)
  names(fstat) <- "Friedman F"
  df2 <- df1 * (N - 1)
  fpvalue <- pf(fstat, df1=df1, df2=df2, lower.tail = FALSE)
  print(fstat)
  print(fpvalue)
  
  cres_nem <- kwAllPairsNemenyiTest(cres[, algorithms])
  print(cres_nem)
  print(algorithms)
}

cres <- comp_results %>%
  dplyr::select(instance_id, dataset, !! qmeas, algorithm)
ggplot(data = cres # get_meds_for_plotting("stability.tt.")
       , aes(y = wxcoverage.tt.#m
             #, ymin = lwr
             #, ymax = upr
             , x = dataset
             , colour = algorithm
             #, group = algorithm
             , alpha = algorithm
             , shape = algorithm)) +
  myGgTheme +
  geom_boxplot(position = position_dodge()) +
  #geom_point() +
  #geom_errorbar(width = 0.2) +
  colos1 + 
  alpas

ggplot(data = comp_results, 
       aes(y = stability.tt., x = wxcoverage.tt.
           , colour = algorithm)) +
  geom_point(alpha = 0.1) +
  colos1 +
  myGgTheme_facets +
  labs(y = get_ylabel("stability.tt."), x = get_ylabel("wxcoverage.tt.")) +
  facet_grid(algorithm ~ dataset)

ggplot(data = comp_results, 
       aes(y = stability.tt., x = wxcoverage.tt.
           , colour = algorithm)) +
  geom_point(alpha = 0.1, size = 2) +
  colos1 +
  myGgTheme_facets +
  labs(y = get_ylabel("stability.tt."), x = get_ylabel("wxcoverage.tt.")) +
  scale_x_continuous(breaks = c(0.25, 0.75)) +
  scale_y_continuous(breaks = c(0.25, 0.75)) + 
  facet_grid(algorithm ~ dataset) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 30)
  )

ggplot(data = comp_results, 
       aes(y = rule.length, x = stability.tt.
           , colour = algorithm)) +
  geom_point(alpha = 0.1, size = 2) +
  colos1 +
  myGgTheme_facets +
  labs(x = get_ylabel("stability.tt."), y = get_ylabel("rule.length")) +
  scale_x_continuous(breaks = c(0.25, 0.75)) + 
  facet_grid(algorithm ~ dataset) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 30)
  )

ggplot(data = comp_results, 
       aes(y = rule.length, x = coverage.tt.
           , colour = algorithm)) +
  geom_point(alpha = 0.1, size = 2) +
  colos1 +
  myGgTheme_facets +
  labs(x = get_ylabel("wxcoverage.tt."), y = get_ylabel("rule.length")) +
  scale_x_continuous(breaks = c(0.25, 0.75)) + 
  facet_grid(algorithm ~ dataset) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 30)
  )
