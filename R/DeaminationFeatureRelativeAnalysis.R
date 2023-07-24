library(data.table)
library(stringr)
library(ggplot2)

THREE_PRIME = "three_prime"
FIVE_PRIME = "five_prime"
HALF_MAX = "half_max"

FEATURE_RELATIVE_POSITION = "Feature_Relative_Mismatch_Pos"
FEATURE_RELATIVE_THREE_PRIME_CUT_SITE = "Feature_Relative_Three_Prime_Cut_Site"
FEATURE_RELATIVE_FIVE_PRIME_CUT_SITE = "Feature_Relative_Five_Prime_Cut_Site"
MEAN_THREE_PRIME_CUT_SITE_DISTANCE = "Mean_Three_Prime_Cut_Site_Distance"
MEAN_FIVE_PRIME_CUT_SITE_DISTANCE = "Mean_Five_Prime_Cut_Site_Distance"



# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 22, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18, color = "black"),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Takes a file of mismatches related to genome features and extracts information on relative mismatch and cut site location.
parseMismatchesByFeature = function(table) {
  table = copy(table)
  colnames(table) = c("Chromosome", "Mismatch_Pos_0", "Mismatch_Pos_1", "Read_Relative_Mismatch_Pos", "Mismatch_Type",
                      "Mismatch_Strand", "Read_Length", "Feature_Pos", "Feature_Strand")

  # Record mismatch positions relative to the transcription factor binding site
  table[Feature_Strand == '+', Feature_Relative_Mismatch_Pos := Mismatch_Pos_0 - Feature_Pos]
  table[Feature_Strand == '-', Feature_Relative_Mismatch_Pos := Feature_Pos - Mismatch_Pos_0]

  # Record cut sites relative to the transcription factor binding site
  table[Mismatch_Strand == '+', Three_Prime_Cut_Site := Mismatch_Pos_0 - Read_Relative_Mismatch_Pos - 0.5]
  table[Mismatch_Strand == '-', Three_Prime_Cut_Site := Mismatch_Pos_0 + Read_Relative_Mismatch_Pos + 0.5]

  table[Mismatch_Strand == '+', Five_Prime_Cut_Site := Mismatch_Pos_0 - (Read_Length + Read_Relative_Mismatch_Pos) - 0.5]
  table[Mismatch_Strand == '-', Five_Prime_Cut_Site := Mismatch_Pos_0 + Read_Length + Read_Relative_Mismatch_Pos + 0.5]

  table[Feature_Strand == '+', Feature_Relative_Three_Prime_Cut_Site := Three_Prime_Cut_Site - Feature_Pos]
  table[Feature_Strand == '-', Feature_Relative_Three_Prime_Cut_Site := Feature_Pos - Three_Prime_Cut_Site]

  table[Feature_Strand == '+', Feature_Relative_Five_Prime_Cut_Site := Five_Prime_Cut_Site - Feature_Pos]
  table[Feature_Strand == '-', Feature_Relative_Five_Prime_Cut_Site := Feature_Pos - Five_Prime_Cut_Site]

  table[,Same_Strand := Mismatch_Strand == Feature_Strand]

}


# This functions take a parsed table of feature-relative mismatch data and return a table of
# the relevant feature-relative data, as specified by the "feature" parameter.
getFeatureRelativeCounts = function(table, dataType, strandAlign = FALSE) {

  if (dataType == FEATURE_RELATIVE_POSITION) {
    if (strandAlign) {
      table[Same_Strand == FALSE, Feature_Relative_Mismatch_Pos := -Feature_Relative_Mismatch_Pos]
      table = table[,.N, by = list(Feature_Relative_Mismatch_Pos)]
    } else {
      table = table[,.N, by = list(Feature_Relative_Mismatch_Pos, Same_Strand)]
    }
  } else if (dataType == FEATURE_RELATIVE_THREE_PRIME_CUT_SITE){
    if (strandAlign) {
      table[Same_Strand == FALSE, Feature_Relative_Three_Prime_Cut_Site := -Feature_Relative_Three_Prime_Cut_Site]
      table = table[,.N, by = list(Feature_Relative_Three_Prime_Cut_Site)]
    } else {
      table = table[,.N, by = list(Feature_Relative_Three_Prime_Cut_Site, Same_Strand)]
    }
  } else if (dataType == FEATURE_RELATIVE_FIVE_PRIME_CUT_SITE) {
    if (strandAlign) {
      table[Same_Strand == FALSE, Feature_Relative_Five_Prime_Cut_Site := -Feature_Relative_Five_Prime_Cut_Site]
      table = table[,.N, by = list(Feature_Relative_Five_Prime_Cut_Site)]
    } else {
      table = table[,.N, by = list(Feature_Relative_Five_Prime_Cut_Site, Same_Strand)]
    }
  } else if (dataType == MEAN_THREE_PRIME_CUT_SITE_DISTANCE) {
    if (strandAlign) {
      table[Same_Strand == FALSE, Feature_Relative_Mismatch_Pos := -Feature_Relative_Mismatch_Pos]
      table = table[,.(Mean_Three_Prime_Cut_Site_Distance = mean(abs(Read_Relative_Mismatch_Pos))),
                    by = list(Feature_Relative_Mismatch_Pos)]
    } else {
      table = table[,.(Mean_Three_Prime_Cut_Site_Distance = mean(abs(Read_Relative_Mismatch_Pos))),
                    by = list(Feature_Relative_Mismatch_Pos, Same_Strand)]
    }
  } else if (dataType == MEAN_FIVE_PRIME_CUT_SITE_DISTANCE) {
    if (strandAlign) {
      table[Same_Strand == FALSE, Feature_Relative_Mismatch_Pos := -Feature_Relative_Mismatch_Pos]
      table = table[,.(Mean_Five_Prime_Cut_Site_Distance = mean(Read_Length-abs(Read_Relative_Mismatch_Pos)+1)),
                    by = list(Feature_Relative_Mismatch_Pos)]
    } else {
      table = table[,.(Mean_Five_Prime_Cut_Site_Distance = mean(Read_Length-abs(Read_Relative_Mismatch_Pos)+1)),
                    by = list(Feature_Relative_Mismatch_Pos, Same_Strand)]
    }
  } else stop("Unrecognized value for feature parameter")

  if (strandAlign) {

    if (dataType == MEAN_THREE_PRIME_CUT_SITE_DISTANCE || dataType == MEAN_FIVE_PRIME_CUT_SITE_DISTANCE) {
      positionCol = FEATURE_RELATIVE_POSITION
      dataCol = dataType
    } else {
      positionCol = dataType
      dataCol = "N"
    }

    return(table[,c(..positionCol, ..dataCol)])

  } else {
    return(table)
  }

}


# This function plots counts of feature-relative features, as specified by the "feature" parameter.
plotFeatureRelativeCounts = function(countsTable, countsDataType,
                                     meanCutSiteDistanceTable = NULL, cutSiteStrandPolarity = NULL,
                                     title = NULL, sameStrand = NULL, xlim = NULL, xAxisBreaks = waiver(),
                                     relativeFreqMax = NULL, meanCutSiteDistanceMin = NULL, meanCutSiteDistanceMax = NULL,
                                     secondaryAxisBreaks = waiver()) {

  countsTable = copy(countsTable)
  if (!is.null(sameStrand)) {
    countsTable = countsTable[Same_Strand == sameStrand]
  }
  countsTable[, Relative_Frequency := N/sum(N)]

  if (is.null(relativeFreqMax)) {relativeFreqMax = max(countsTable$Relative_Frequency) * 1.1}

  if (countsDataType == FEATURE_RELATIVE_POSITION) {
    xAxisLabel = "Feature Relative Position"
    if (is.null(title)) title = "Feature Relative Mismatch Position Frequencies"
  } else if (countsDataType == FEATURE_RELATIVE_THREE_PRIME_CUT_SITE){
    xAxisLabel = "Feature Relative 3' Cut Site"
    if (is.null(title)) title = "Feature Relative 3' Cut Site Position Frequencies"
  } else if (countsDataType == FEATURE_RELATIVE_FIVE_PRIME_CUT_SITE) {
    xAxisLabel = "Feature Relative 5' Cut Site"
    if (is.null(title)) title = "Feature Relative 5' Cut Site Position Frequencies"
  } else stop("Unrecognized value for data type parameter")

  if (!is.null(meanCutSiteDistanceTable)) {

    if (cutSiteStrandPolarity == THREE_PRIME) {
      meanCutSiteDistanceCol = MEAN_THREE_PRIME_CUT_SITE_DISTANCE
      secondaryAxisLabel = "Mean 3' Cut Site Distance"
      lineColor = "red2"
    } else if (cutSiteStrandPolarity == FIVE_PRIME) {
      meanCutSiteDistanceCol = MEAN_FIVE_PRIME_CUT_SITE_DISTANCE
      secondaryAxisLabel = "Mean 5' Cut Site Distance"
      lineColor = "deepskyblue3"
    } else {
      stop("Unrecognized value for cut site strand polarity")
    }

    countsTable = merge.data.table(countsTable, meanCutSiteDistanceTable, by = FEATURE_RELATIVE_POSITION)

    adjustmentFactor = (max(countsTable[[meanCutSiteDistanceCol]]) - min(countsTable[[meanCutSiteDistanceCol]]))*0.05
    if (is.null(meanCutSiteDistanceMin)) {meanCutSiteDistanceMin = min(countsTable[[meanCutSiteDistanceCol]]) - adjustmentFactor}
    if (is.null(meanCutSiteDistanceMax)) {meanCutSiteDistanceMax = max(countsTable[[meanCutSiteDistanceCol]]) + adjustmentFactor}
    yAxisRatio = relativeFreqMax/(meanCutSiteDistanceMax-meanCutSiteDistanceMin)

  }

  plot = ggplot(countsTable, aes(x = !!sym(countsDataType), y = Relative_Frequency)) +
                geom_bar(stat = "identity") + coord_cartesian(xlim = xlim, expand = FALSE) +
                scale_x_continuous(breaks = xAxisBreaks) +
                labs(title = title, x = xAxisLabel, y = "Relative Frequency") +
                blankBackground + defaultTextScaling

  if (!is.null(meanCutSiteDistanceTable)) {
    plot = plot +
      geom_line(aes(x = !!sym(FEATURE_RELATIVE_POSITION),
                    y = (countsTable[[meanCutSiteDistanceCol]] - meanCutSiteDistanceMin) * yAxisRatio ),
                linewidth = 1, color = lineColor) +
      scale_y_continuous(breaks = secondaryAxisBreaks, limits = c(0,relativeFreqMax), expand = expansion(),
                         sec.axis = sec_axis(~. / yAxisRatio + meanCutSiteDistanceMin, name = secondaryAxisLabel)) +
      theme(axis.text.y.right = element_text(size = 18, color = lineColor),
            axis.title.y.right = element_text(size = 22, color = lineColor))

  } else {

    plot = plot + scale_y_continuous(limits = c(0,relativeFreqMax), expand = expansion())

  }

  print(plot)

}
