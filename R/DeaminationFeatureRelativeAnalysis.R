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
COUNTS = "N"
RELATIVE_FREQUENCY = "Relative_Frequency"

LINE = "line"
BAR = "bar"



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

  if (dataType == MEAN_THREE_PRIME_CUT_SITE_DISTANCE || dataType == MEAN_FIVE_PRIME_CUT_SITE_DISTANCE) {
    positionCol = FEATURE_RELATIVE_POSITION
    dataCol = dataType
  } else {
    positionCol = dataType
    dataCol = "N"
  }

  fullRange = data.table(Position = min(table[[positionCol]]):max(table[[positionCol]]))
  table = merge.data.table(table, fullRange, by.x = positionCol, by.y = "Position", all = TRUE)
  if (dataCol == "N") {table[is.na(N), N := 0]}

  if (strandAlign) {
    return(table[order(table[[positionCol]]),c(..positionCol, ..dataCol)])
  } else {
    return(table[order(table[[positionCol]])])
  }

}


# This function plots counts of feature-relative features, as specified by the "feature" parameter.
plotFeatureRelativeCounts = function(countsTable, countsDataType, countsYVar = COUNTS,
                                     plotCountsMirror = FALSE, smoothCountsData = FALSE,
                                     countsPlotType = LINE, countsPlotColor = "black", countsMirrorColor = "gray",
                                     meanCutSiteDistanceTable = NULL, cutSiteStrandPolarity = NULL,
                                     smoothMeanCutSiteDistance = FALSE,
                                     title = NULL, sameStrand = NULL, xlim = NULL, xAxisBreaks = waiver(),
                                     mainYAxisMin = NULL, mainYAxisMax = NULL, mainYAxisLabel = NULL,
                                     meanCutSiteDistanceMin = NULL, meanCutSiteDistanceMax = NULL,
                                     secondaryAxisBreaks = waiver()) {

  countsTable = copy(countsTable)
  if (!is.null(sameStrand)) { countsTable = countsTable[Same_Strand == sameStrand | is.na(sameStrand)] }
  if (countsYVar == RELATIVE_FREQUENCY) { countsTable[, Relative_Frequency := N/sum(N)] }

  if (smoothCountsData) {
    countsTable[,c(paste0(countsYVar,"_Smoothed")) := frollmean(countsTable[[countsYVar]], 5, align = "center", na.rm = TRUE)]
    countsYVar = paste0(countsYVar,"_Smoothed")
  }
  countsTable = countsTable[complete.cases(countsTable)]

  adjustmentFactor = (max(countsTable[[countsYVar]]) - min(countsTable[[countsYVar]]))*0.05
  if (is.null(mainYAxisMax)) {mainYAxisMax = max(countsTable[[countsYVar]]) + adjustmentFactor}
  if (is.null(mainYAxisMin)) {
    mainYAxisMin = min(countsTable[[countsYVar]])
    if (mainYAxisMin != 0) mainYAxisMin = mainYAxisMin - adjustmentFactor
  }

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

    if (!is.null(sameStrand)) { meanCutSiteDistanceTable = meanCutSiteDistanceTable[Same_Strand == sameStrand | is.na(sameStrand)] }

    if (smoothMeanCutSiteDistance) {
      meanCutSiteDistanceTable[,c(paste0(meanCutSiteDistanceCol,"_Smoothed")) :=
                                 frollmean(meanCutSiteDistanceTable[[meanCutSiteDistanceCol]], 5, align = "center", na.rm = TRUE)]
      meanCutSiteDistanceCol = paste0(meanCutSiteDistanceCol,"_Smoothed")
    }
    meanCutSiteDistanceTable = meanCutSiteDistanceTable[complete.cases(meanCutSiteDistanceTable)]

    countsTable = merge.data.table(countsTable, meanCutSiteDistanceTable, by = FEATURE_RELATIVE_POSITION)

    adjustmentFactor = (max(countsTable[[meanCutSiteDistanceCol]]) - min(countsTable[[meanCutSiteDistanceCol]]))*0.05
    if (is.null(meanCutSiteDistanceMin)) {meanCutSiteDistanceMin = min(countsTable[[meanCutSiteDistanceCol]]) - adjustmentFactor}
    if (is.null(meanCutSiteDistanceMax)) {meanCutSiteDistanceMax = max(countsTable[[meanCutSiteDistanceCol]]) + adjustmentFactor}
    yAxisRatio = (mainYAxisMax-mainYAxisMin)/(meanCutSiteDistanceMax-meanCutSiteDistanceMin)

  }

  plot = ggplot(countsTable, aes(x = !!sym(countsDataType), y = !!sym(countsYVar)))

  if (countsPlotType == BAR) {
    plot = plot + geom_bar(stat = "identity") + coord_cartesian(xlim = xlim, expand = FALSE, fill = countsPlotColor)
  } else if (countsPlotType == LINE) {
    plot = plot + geom_line(linewidth = 1, color = countsPlotColor)
  }
  if (is.null(mainYAxisLabel)) {mainYAxisLabel = countsYVar}
  plot = plot +
    scale_x_continuous(breaks = xAxisBreaks) +
    labs(title = title, x = xAxisLabel, y = mainYAxisLabel) +
    blankBackground + defaultTextScaling +
    theme(axis.text.y.left = element_text(size = 18, color = countsPlotColor),
          axis.title.y.left = element_text(size = 22, color = countsPlotColor))

  if (plotCountsMirror) {
    mirroredCountsTable = countsTable[countsTable[[countsDataType]] <= 0, c(..countsDataType, ..countsYVar)]
    mirroredCountsTable[,c(countsDataType) := -mirroredCountsTable[[countsDataType]]]
    plot = plot + geom_line(aes(x = !!sym(countsDataType), y = !!sym(countsYVar)), data = mirroredCountsTable,
                            linewidth = 1, color = countsMirrorColor)
  }

  if (!is.null(meanCutSiteDistanceTable)) {
    plot = plot +
      geom_line(aes(x = !!sym(FEATURE_RELATIVE_POSITION),
                    y = countsTable[[meanCutSiteDistanceCol]] * yAxisRatio - (meanCutSiteDistanceMax*yAxisRatio-mainYAxisMax)),
                linewidth = 1, color = lineColor) +
      scale_y_continuous(breaks = secondaryAxisBreaks, limits = c(mainYAxisMin,mainYAxisMax), expand = expansion(),
                         sec.axis = sec_axis(~. / yAxisRatio + (meanCutSiteDistanceMax-mainYAxisMax/yAxisRatio),
                                             name = secondaryAxisLabel)) +
      theme(axis.text.y.right = element_text(size = 18, color = lineColor),
            axis.title.y.right = element_text(size = 22, color = lineColor))

  } else {

    plot = plot + scale_y_continuous(limits = c(mainYAxisMin,mainYAxisMax), expand = expansion())

  }

  print(plot)

}