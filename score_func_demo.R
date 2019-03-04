library(tikzDevice)

# score func demo
set.seed(123)
sup <- runif(10) # over-written
# sup <- seq(0.05, 0.95, length.out = 10)
len <- rpois(10, 2) + 1
scr <- function(s, l, a) {
  s * (l - a)/l
}
a_vec <- seq(0, 1, 0.1)
scrs <- t(sapply(a_vec, function(x) {mapply(FUN = scr, sup, len2, a = x)}))

tikz(file = "score_func_sim.tikz", height = 4, width = 7)
matplot(scrs, pch = as.character(len2)
        , type = "o"
        , col = rainbow(ncol(scrs))
        , xaxt = "n"
        , main = "Score function simulation"
        , xlab = expression(alpha)
        , ylab = "Score")
axis(1, at=seq_along(a_vec), labels = a_vec)
text(8, 0.9, labels = "N = rule length", col = 2)
dev.off()

# abline(v = which(a_vec == 0))
