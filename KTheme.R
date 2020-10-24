library(ggplot2)
library(grid)
library(gridExtra)
library(tikzDevice)
library(cowplot) # get legend out of a ggplot
library(lattice)

# Kaplan Theme

# brand
k.purple <- "#1F0477"
k.darkblue <- "#00539F"
k.brightblue <- "#0098CD"
k.darkred <- "#D6083B"
k.orange <- "#E55302"
k.deepgreen <- "#3D9B35"
k.lime <- "#BED600"
k.pink <- "#E62899"
k.brightorange <- "#FFC82E"

# neutral
k.warmgrey9 <- "#837870"
k.warmgrey2 <- "#D7D2CB"
k.coolgrey9 <- "#747679"
k.coolgrey4 <- "#BDBDBC"
k.wheat <- "#E8ED70"
k.stone <- "#DCCEAB"

myPalBrand <- c(k.purple, k.darkblue, k.brightblue
           , k.darkred, k.orange, k.deepgreen
           , k.lime, k.pink, k.brightorange)
myPalBrand2 <- c(k.purple, k.orange, k.deepgreen
                 , k.darkblue, k.darkred, k.lime
                 , k.brightblue, k.pink, k.brightorange)
myPalBrandContrasts <- c(k.purple, k.brightorange, k.lime
                    , k.darkblue, k.pink, k.deepgreen
                    , k.orange, k.brightblue, k.darkred)

myPalNeut <- c(k.warmgrey9, k.warmgrey2, k.coolgrey9, k.coolgrey4, k.wheat, k.stone,  "#414346")
myPalNeut2 <- c(k.coolgrey4, k.coolgrey9, k.warmgrey9, k.warmgrey2, k.stone, k.wheat)

grad.purp <- c("#6D066D", "#8C6DD6")
k.grad.purple <- colorRampPalette(grad.purp)
k.grad.purple.rev <- colorRampPalette(rev(grad.purp))

grad.red <- c("#A80212", "#F71043")
k.grad.red <- colorRampPalette(grad.red)
k.grad.red.rev <- colorRampPalette(rev(grad.red))

grad.grey <- c("#2D2E23", "#AAABAD")
k.grad.grey <- colorRampPalette(grad.grey)
k.grad.grey.rev <- colorRampPalette(rev(grad.grey))

grad.blue <- c("#0071B2", "#10C8F4")
k.grad.blue <- colorRampPalette(grad.blue)
k.grad.blue.rev <- colorRampPalette(rev(grad.blue))

grad.orange <- c("#DE4B00", "#FFA241")
k.grad.orange <- colorRampPalette(grad.orange)
k.grad.orange.rev <- colorRampPalette(rev(grad.orange))

grad.green <- c("#2A5E02", "#B3E879")
k.grad.green <- colorRampPalette(grad.green)
k.grad.green.rev <- colorRampPalette(rev(grad.green))

# applied to lattice
MyLatticeFont <- list(font = 1, cex = 1)
MyLabelFont <- list(font = 1, cex = 0.9)
MyStripFont <- list(font = 1, cex = 0.8)
MyLatticeStrip = strip.custom(par.strip.text = MyStripFont)
MyLatticeTheme <- list(
  par.main.text = MyLatticeFont
  , par.xlab.text = MyLabelFont
  , par.ylab.text = MyLabelFont
  , add.text = MyLabelFont
  , axis.text = MyStripFont
  , fontsize = list(text = 11, points = 7)
  , plot.symbol = list(col = k.purple, pch = 19, alpha = 1, cex = 1)
  , plot.polygon = list(col = k.purple)
  , superpose.symbol = list(pch = c(20, 18, 16, 10, 12, 3, 4, 6, 2)
                            , col = myPalBrand)
  , superpose.polygon = list(col = myPalBrand2)
  , box.umbrella = list(col = k.purple, lty = 2, lwd = 1.25)
  , box.rectangle = list(col = k.purple, lwd = 1.5)
  , box.dot = list(col = k.purple, pch = 16, cex = 1.5, alpha = 1)
  , plot.line = list(myPalBrand)
  , add.line = list(col = myPalBrand2, lty = 2)
  , strip.background = list(col = myPalNeut[c(2,4,6)])
  , strip.shingle = list(col = myPalNeut[c(1,3,5)])
)

# fourfold colour scheme
fourfold_k <- c(k.darkred, k.purple, k.pink, k.brightblue, k.darkred, k.purple)

# strucplot grapcons
# TO DO implement interpolation
shading_k <- function(observed = NULL, residuals = NULL, expected = NULL,
                      df = NULL, col = c(k.purple, k.pink)) {
  if (length(col) != 2) { stop("Need exactly two colors!") }
  function(res) gpar(fill = ifelse(res > 0, col[1], col[2]))
}
class(shading_k) <- "grapcon_generator"

# applied to ggplot2
myGgTheme <- theme(plot.title = element_text(colour = myPalNeut[7], size = 10)
                   , axis.title = element_text(colour = myPalNeut[7], size = 10)
                   , axis.text = element_text(colour = myPalNeut[7], size = 7)
                   , axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 15, b = 0, l = 0))
                   , axis.line = element_line(colour = myPalNeut[7], size = 0.5)
                   , axis.ticks = element_line(colour = myPalNeut[7], size = 0.5)
                   , panel.background = element_blank()
                   , panel.border = element_blank()
                   , legend.title = element_text(colour = myPalNeut[7], size = 8)
                   , legend.text = element_text(colour = myPalNeut[7], size = 6)
                   , panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()
                   , legend.background = element_blank()
                   , legend.box.background = element_blank()
                   , legend.key = element_blank()
)

myGgTheme_facets <- theme(plot.title = element_text(colour = myPalNeut[7], size = 10)
                          , axis.title = element_text(colour = myPalNeut[7], size = 10)
                          , axis.text = element_text(colour = myPalNeut[7], size = 7)
                          , axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 15, b = 0, l = 0))
                          , axis.line = element_line(colour = myPalNeut[7], size = 0.5)
                          , axis.ticks = element_line(colour = myPalNeut[7], size = 0.5)
                          , panel.background = element_rect(colour = myPalNeut[7], fill = "#EEEEEEEE")
                          , legend.title = element_text(colour = myPalNeut[7], size = 8)
                          , legend.text = element_text(colour = myPalNeut[7], size = 6)
                          , panel.grid.major = element_blank()
                          , panel.grid.minor = element_blank()
                          , legend.background = element_blank()
                          , legend.box.background = element_blank()
                          , legend.key = element_blank()
)

myGgFillScale <- scale_fill_manual(values = c(myPalBrand[8], myPalBrand[2]))
myGgFillScaleBlueOrange <- scale_fill_manual(values = c(k.brightblue, k.brightorange))
myGgColourScale <- scale_colour_manual(values = c(myPalBrand[8], myPalBrand[2]))
myGgColourScaleBlueOrange <- scale_colour_manual(values = c(k.brightblue, k.brightorange))
myGgPinkGradient <- scale_color_gradient(low = myPalBrand[1], high = myPalBrand[1])
myGgSeaGradient <- scale_color_gradient(low = myPalBrand[2], high = myPalBrand[2])
myGgSapphireGradient <- scale_color_gradient(low = myPalBrand[3], high = myPalBrand[3])
myGgHeatGradient <- scale_color_gradient2(low = myPalBrand[3]
                                          , mid = myPalBrand[2]
                                          , high = myPalBrand[1]
                                          , midpoint = 3000)
