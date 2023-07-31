library(data.table)
library(ggplot2)

# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 22, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18, color = "black"),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


getPairedCutSiteDistances = function(mismatchesByRead, filteredAndFormattedInput = TRUE) {

  if (filteredAndFormattedInput) {
    colnames(mismatchesByRead) = c("Chromosome", "Mismatch_Pos_0", "Mismatch_Pos_1", "Three_Prime_Relative_Mismatch_Pos",
                                   "Mismatch_Type", "Mismatch_Strand", "Read_Length")
  } else {
    colnames(mismatchesByRead) = c("Chromosome", "Read_Pos_0", "Read_Pos_1", "Three_Prime_Relative_Mismatch_Pos",
                                   "Mismatch_Type", "Mismatch_Strand", "Read_Sequence")

    mismatchesByRead[,Read_Length := nchar(Read_Sequence)]
  }

  mismatchesByRead[,Three_Prime_Cut_Site_Distance := abs(Three_Prime_Relative_Mismatch_Pos)]
  mismatchesByRead[,Five_Prime_Cut_Site_Distance := Read_Length-abs(Three_Prime_Relative_Mismatch_Pos)+1]

  return(mismatchesByRead[,.(Three_Prime_Cut_Site_Distance, Five_Prime_Cut_Site_Distance)])

}

getCorrelationResults = function(pairedCutSiteDistances) {

  return(cor.test(pairedCutSiteDistances$Three_Prime_Cut_Site_Distance,
                  pairedCutSiteDistances$Five_Prime_Cut_Site_Distance, method = "pearson"))

}

plotPairedCutSiteDistances = function(pairedCutSiteDistances,
                                       xAxisLabel = "3' Cut Site Distance", yAxisLabel = "5' CutSite Distance",
                                       title = "Paired Cut Site Distances") {

  plot = ggplot(pairedCutSiteDistances, aes(x = Three_Prime_Cut_Site_Distance, y = Five_Prime_Cut_Site_Distance)) +
    geom_count() + blankBackground + defaultTextScaling +
    labs(title = title, x = xAxisLabel, y = yAxisLabel)

  print(plot)

}
