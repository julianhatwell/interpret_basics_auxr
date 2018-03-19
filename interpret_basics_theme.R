# Theme for TC report
# a nice colour palette
myPal <- c("#8DD3C7", "#B0A8B3", "#9FB5D6", "#9EC0FA", "#DB8072")
myPalDark <- c("#4D8377", "#504853", "#3F5576", "#3E40DA", "#AB2012")
myPalContrasts <- c(myPalDark[1], myPalDark[5], myPal[2], myPal[4]
                    ,"#999999"
                    , myPal[1], myPal[5], myPalDark[2], myPalDark[4])
myPal.range <- colorRampPalette(c("#FFFFFF", myPal[3:1]))
myPal.rangeDark <- colorRampPalette(c("#FFFFFF", myPalDark[3:1]))
myPal.rangeDiv <- colorRampPalette(c(myPal[1], "#FFFFFF", myPal[5]))

# fourfold 6
myPalFourFold <- c(myPalContrasts[6], myPalContrasts[1]
                   , myPalContrasts[7], myPalContrasts[4]
                   , myPalContrasts[2], myPalContrasts[9])

# applied to lattice
MyLatticeFont <- list(font = 1, fontsize = 12, col = myPalDark[2])
MyStripFont <- list(font = 1, cex = 0.8, col = myPalDark[2])

MyLatticeTheme <- list(
  par.main.text = MyLatticeFont
  , par.xlab.text = MyLatticeFont
  , par.ylab.text = MyLatticeFont
  , axis.text = MyLatticeFont
  , add.text = MyLatticeFont
  , fontsize = list(text = 11, points = 7)
  , plot.symbol = list(col = myPalDark[5], pch = 19, alpha = 0.5, cex = 0.75)
  , box.umbrella = list(col = myPalDark[5], lty = 2, lwd = 1.25)
  , box.rectangle = list(fill = myPal.range(100)[3], col = myPalDark[5], lwd = 1.5)
  , box.dot = list(col = myPalDark[5], pch = 15, cex = 0.8, alpha = 0.5)
  , superpose.line = list(col = myPalContrasts)
  , superpose.symbol = list(col = myPalContrasts)
  , axis.line = list(col = myPal[2])
  , strip.background = list(col=myPal.range(100)[3])
  , strip.shingle = list(col=myPal.range(100)[12])
  , strip.border = list(col = myPal[2])
)
MyLatticeStrip = strip.custom(par.strip.text = MyLatticeFont)

# applied to ggplot2
myGgTheme <- theme(plot.title = element_text(colour = myPalDark[2], size = 10)
                   , axis.title = element_text(colour = myPalDark[2], size = 10)
                   , axis.text = element_text(colour = myPalDark[2], size = 10)
                   , legend.title = element_text(colour = myPalDark[2])
                   , legend.text = element_text(colour = myPalDark[2]))

myGgThemeSilentX <- theme(plot.title = element_text(colour = myPalDark[2], size = 10)
                          , axis.title = element_text(colour = myPalDark[2], size = 10)
                          , axis.text.y = element_text(colour = myPalDark[2], size = 10)
                          , axis.text.x = element_blank()
                          , legend.title = element_text(colour = myPalDark[2])
                          , legend.text = element_text(colour = myPalDark[2]))

myGgThemeSilentY <- theme(plot.title = element_text(colour = myPalDark[2], size = 10)
                          , axis.title = element_text(colour = myPalDark[2], size = 10)
                          , axis.text.x = element_text(colour = myPalDark[2], size = 10)
                          , axis.text.y = element_blank()
                          , legend.title = element_text(colour = myPalDark[2])
                          , legend.text = element_text(colour = myPalDark[2]))
