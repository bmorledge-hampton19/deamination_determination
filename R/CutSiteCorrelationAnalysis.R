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

THREE_PRIME_CUT_SITE_DISTANCE = "Three_Prime_Cut_Site_Distance"
READ_LENGTH = "Read_Length"

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

  return(mismatchesByRead[,.(Three_Prime_Cut_Site_Distance, Five_Prime_Cut_Site_Distance, Read_Length)])

}

getShufflePairedCutSiteDistances = function(pairedCutSiteDistances) {
  pairedCutSiteDistances = copy(pairedCutSiteDistances)
  pairedCutSiteDistances[,Three_Prime_Cut_Site_Distance :=
                           sample(Three_Prime_Cut_Site_Distance, nrow(pairedCutSiteDistances))]
  pairedCutSiteDistances[,Read_Length := Three_Prime_Cut_Site_Distance + Five_Prime_Cut_Site_Distance - 1]
  return(pairedCutSiteDistances)
}

getCorrelationResults = function(pairedCutSiteDistances, dependentVariable) {

  return(cor.test(pairedCutSiteDistances[[dependentVariable]],
                  pairedCutSiteDistances$Five_Prime_Cut_Site_Distance, method = "pearson"))

}

getCouplingRate = function(pairedCutSiteDistances) {
  linearModel = lm(Three_Prime_Cut_Site_Distance ~ Five_Prime_Cut_Site_Distance, pairedCutSiteDistances)
  return(linearModel$coefficients[2])
}

getReadLengthSlope = function(pairedCutSiteDistances) {
  linearModel = lm(Read_Length ~ Five_Prime_Cut_Site_Distance, pairedCutSiteDistances)
  return(linearModel$coefficients[2])
}

plotPairedCutSiteDistances = function(pairedCutSiteDistances, dependentVariable,
                                      xAxisLabel = "5' Cut Site Distance", yAxisLabel = NULL,
                                      title = "Paired Cut Site Distances", sizeRange = c(0.1,5),
                                      xlim = NULL, ylim = NULL, xAxisBreaks = waiver(), yAxisBreaks = waiver(),
                                      breaks = waiver()) {

  if (is.null(yAxisLabel)) {
    if (dependentVariable == THREE_PRIME_CUT_SITE_DISTANCE) {
      yAxisLabel = "3' Cut Site Distance"
    } else if (dependentVariable == READ_LENGTH) {
      yAxisLabel = "Read Length"
    }
  }

  plot = ggplot(pairedCutSiteDistances, aes(x = Five_Prime_Cut_Site_Distance, y = !!sym(dependentVariable))) +
    geom_count() + scale_size(range = sizeRange, breaks = breaks) +
    geom_smooth(method = "lm", formula = y~x, color = "red2", se = FALSE) +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_x_continuous(breaks = xAxisBreaks) + scale_y_continuous(breaks = yAxisBreaks) + blankBackground + defaultTextScaling

  print(plot)

}
