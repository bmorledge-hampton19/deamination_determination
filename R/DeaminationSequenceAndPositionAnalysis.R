library(data.table)
library(stringr)
library(ggplot2)

THREE_PRIME = "three_prime"
FIVE_PRIME = "five_prime"

# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 22, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Computes and returns sequence length and GC content for each read in a given table.
getReadLengthAndGCContent = function(table) {
  return(table[,.(Length = nchar(V7),
                  GC_Content = str_count(V7, "C|G") / (nchar(V7) - str_count(V7, 'N')) )])
}


# Displays the frequencies of each mismatch type.
plotMismatchTypeFrequencies = function(mismatchTable, title = "Mismatch Type Frequency") {

  print(
    ggplot(mismatchTable[, .N, by = Mismatch][, Freq := N/sum(N)], aes(x = Mismatch, y = Freq)) +
      geom_bar(stat = "identity") +
      labs(title = title, x = "Mismatch Type", y = "Frequency") +
      blankBackground + defaultTextScaling + theme(axis.text = element_text(size = 15))
  )

}

# Displays the frequencies of positions for the chosen mismatch types.
# (Types can be chosen by inclusion or omission)
plotMismatchPositionFrequencies = function(mismatchTable, includedTypes = list(), omittedTypes = list(),
                                           title = "Position Frequency", posType = THREE_PRIME) {

  if (posType == THREE_PRIME) {
    xAxisLabel = "3' Relative Position"
  } else if (posType == FIVE_PRIME) {
    xAxisLabel = "5' Relative Position"
    mismatchTable = copy(mismatchTable)
    mismatchTable[, Position := Read_Length + Position + 1]
  } else stop("Unrecognized value for posType parameter.")

  if ( length(includedTypes) > 0 && length(omittedTypes) > 0) {
    stop("Included types and omitted types given simultaneously")
  } else if (length(includedTypes) > 0) {
    frequencyData = mismatchTable[Mismatch %in% includedTypes, .N, by = Position][, Freq := N/sum(N)]
  } else if (length(omittedTypes) > 0) {
    frequencyData = mismatchTable[!(Mismatch %in% omittedTypes), .N, by = Position][, Freq := N/sum(N)]
  } else {
    frequencyData = mismatchTable[, .N, by = Position][, Freq := N/sum(N)]
  }

  print(
    ggplot(frequencyData, aes(x = Position, y = Freq)) +
      geom_bar(stat = "identity") +
      labs(title = title, x = xAxisLabel, y = "Frequency") +
      blankBackground + defaultTextScaling
  )

}

# Plot mismatch position using a facet plot with timepoint (optional) on one dimension and read length on the other
# Input should be a list of data.tables with names as timepoint information, or a single data.table
# to construct a plot without timepoint information.
plotPositionAcrossTimepointAndReadLength = function(simplifiedTables, includedTypes = list(), omittedTypes = list(),
                                                    title = "Mismatch Position Frequencies", posType = THREE_PRIME,
                                                    zScoreTables = NULL, zScoreCutoff = 4) {

  #If passed a single data.table, wrap it in a list.
  if (is.data.table(simplifiedTables)) {
    simplifiedTables = list(NONE = simplifiedTables)
  }
  simplifiedTables = lapply(simplifiedTables, copy)
  noTimepointInfo = all(names(simplifiedTables) == "NONE")

  if (!is.null(zScoreTables)) {
    zScoreTables = lapply(zScoreTables, copy)
    if (is.data.table(zScoreTables)) {
      zScoreTables = list(None = zScoreTables)
    }
  }

  if (posType == THREE_PRIME) {
    xAxisLabel = "3' Relative Position"
    xAxisBreaks = c(0, -10, -20)
  } else if (posType == FIVE_PRIME) {
    xAxisLabel = "5' Relative Position"
    xAxisBreaks = c(0, 10, 20)
    simplifiedTables = lapply(simplifiedTables, copy)
    lapply(simplifiedTables, function(x) x[, Position := Read_Length + Position + 1])
    if (!is.null(zScoreTables)) {
      lapply(zScoreTables, function(x) x[, Position := Read_Length + Position + 1])
    }
  } else stop("Unrecognized value for posType parameter.")

  aggregateMismatchesTable = rbindlist(lapply(seq_along(simplifiedTables),
                                       function(i) simplifiedTables[[i]][,Timepoint := names(simplifiedTables)[i]]))

  if ( length(includedTypes) > 0 && length(omittedTypes) > 0) {
    stop("Included types and omitted types given simultaneously")
  } else if (length(includedTypes) > 0) {
    aggregateMismatchesTable = aggregateMismatchesTable[Mismatch %in% includedTypes]
  } else if (length(omittedTypes) > 0) {
    aggregateMismatchesTable = aggregateMismatchesTable[!(Mismatch %in% omittedTypes)]
  }

  groupedPositionFrequencies = (aggregateMismatchesTable[, .N, by = list(Position,Read_Length,Timepoint)]
                                [, Frequency := N/sum(N), by = list(Read_Length, Timepoint)])
  if (!is.null(zScoreTables)) {
    aggregateZScoreTable = rbindlist(lapply(seq_along(zScoreTables),
                                            function(i) zScoreTables[[i]][,Timepoint := names(zScoreTables)[i]]))
    setkey(groupedPositionFrequencies,Position,Read_Length,Timepoint)
    setkey(aggregateZScoreTable,Position,Read_Length,Timepoint)
    groupedPositionFrequencies = groupedPositionFrequencies[aggregateZScoreTable]
    groupedPositionFrequencies[,Meets_Cutoff := Z_Score >= zScoreCutoff]
  }

  plot = ggplot(groupedPositionFrequencies, aes(Position, Frequency)) +
    labs(title = title, x = xAxisLabel, y = "Relative Mismatch Frequency") +
    blankBackground + defaultTextScaling

  if (!is.null(zScoreTables)) {
    plot = plot +
      geom_bar(aes(fill = Meets_Cutoff), stat = "identity") +
      scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none")
  } else {
    plot = plot + geom_bar(stat = "identity")
  }

  if (noTimepointInfo) {
    plot = plot + facet_grid(rows = vars(Read_Length))
  } else {
    plot = plot + facet_grid(Read_Length~factor(Timepoint, levels = names(simplifiedTables)))
  }


  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          strip.background = element_rect(color = "black", size = 1),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y = element_text(size = 16),
          panel.grid.major.x = element_line(color = "black", size = 0.5, linetype = 2)) +
    scale_y_continuous(sec.axis = dup_axis(~., name = "Read Length")) +
    scale_x_continuous(breaks = xAxisBreaks)

  print(plot)

}


# Plot read length frequencies across timepoint (optional).
# Input should be a list of data.tables with names as timepoint information, or a single data.table
# to construct a plot without timepoint information.
plotReadLengthFrequencies = function(simplifiedTables, title = "Read Length Frequencies") {

  #If passed a single data.table, wrap it in a list.
  if (is.data.table(simplifiedTables)) {
    simplifiedTables = list(None = simplifiedTables)
    noTimepointInfo = TRUE
  } else noTimepointInfo = FALSE

  aggregateTable = rbindlist(lapply(seq_along(simplifiedTables),
                                    function(i) simplifiedTables[[i]][,Timepoint := names(simplifiedTables)[i]]))

  readLengthFrequencies = (aggregateTable[, .N, by = list(Read_Length,Timepoint)]
                                         [, Freq := N/sum(N), by = list(Timepoint)])
  maxFrequency = round(max(readLengthFrequencies$Freq), digits = 2)
  yAxisBreaks = c(0, maxFrequency/2, maxFrequency)

  plot = ggplot(readLengthFrequencies, aes(Read_Length, Freq)) +
    geom_bar(stat = "identity") +
    labs(title = title, x = "Read_Length", y = "Frequency") +
    blankBackground + defaultTextScaling

  if (!noTimepointInfo) {
    plot = plot + facet_grid(cols = vars(Timepoint))
  }

  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          strip.background = element_rect(color = "black", size = 1),
          axis.text.y = element_text(size = 16), strip.text.y = element_text(size = 16)) +
    scale_y_continuous(breaks = yAxisBreaks)

  print(plot)

}


# Gets the mode. (Why isn't there a base R function for this?)
getMode = function(values) {
  frequencyTable = table(values)
  return(as.numeric(names(frequencyTable)[which.max(frequencyTable)]))
}


# Returns a table of the mean, median, and mode positions for each read length at each time point.
getGroupedPositionStats = function(simplifiedTables, includedTypes = list(),
                                   omittedTypes = list(), posType = THREE_PRIME) {

  #If passed a single data.table, wrap it in a list.
  if (is.data.table(simplifiedTables)) {
    simplifiedTables = list(None = simplifiedTables)
  }

  if (posType == FIVE_PRIME) {
    simplifiedTables = lapply(simplifiedTables, copy)
    lapply(simplifiedTables, function(x) x[, Position := Read_Length + Position + 1])
  } else if (posType != THREE_PRIME) stop("Unrecognized value for posType parameter.")

  aggregateTable = rbindlist(lapply(seq_along(simplifiedTables),
                                    function(i) simplifiedTables[[i]][,Timepoint := names(simplifiedTables)[i]]))

  if ( length(includedTypes) > 0 && length(omittedTypes) > 0) {
    stop("Included types and omitted types given simultaneously")
  } else if (length(includedTypes) > 0) {
    aggregateTable = aggregateTable[Mismatch %in% includedTypes]
  } else if (length(omittedTypes) > 0) {
    aggregateTable = aggregateTable[!(Mismatch %in% omittedTypes)]
  }

  groupedPositionStats = aggregateTable[order(Timepoint, Read_Length),
                                        .(mean = mean(Position), median = median(Position),
                                          mode = getMode(Position), IQR = IQR(Position),
                                          Standard_Deviation = sd(Position)),
                                        by = list(Read_Length,Timepoint)]

  # Get the absolute change in position relative to the 3' end for each timepoint
  shortestReadMeanPos = unlist(lapply(unique(groupedPositionStats$Timepoint), function(x) {
    rep(groupedPositionStats[Timepoint == x, mean][1],length(groupedPositionStats[Timepoint==x,mean]))
  }))
  groupedPositionStats[,Absolute_Pos_Change := abs(mean-shortestReadMeanPos)]

  return(groupedPositionStats)

}


POS_DIFF = "Position Difference"
IQR = "IQR"
STDEV = "Standard Deviation"
# Plots grouped stats across timepoints and positioning type (3' vs. 5')
plotGroupedPositionStats = function(threePrimeGroupedStats, fivePrimeGroupedStats, stat,
                                    title = paste(stat,"Over Time"), xAxisBreaks = waiver()) {

  # Combine the two data sets with an extra column for the position type
  threePrimeGroupedStats = copy(threePrimeGroupedStats)
  threePrimeGroupedStats[,Position_Type := "3' Relative Position"]

  fivePrimeGroupedStats = copy(fivePrimeGroupedStats)
  fivePrimeGroupedStats[,Position_Type := "5' Relative Position"]

  aggregateData = rbindlist(list(threePrimeGroupedStats,fivePrimeGroupedStats))

  if (stat == POS_DIFF) {
    groupedStatsPlot =
      ggplot(aggregateData, aes(Read_Length, Absolute_Pos_Change, color = Timepoint, linetype = Position_Type)) +
      geom_line(size = 1.25)
  } else if (stat == IQR) {
    groupedStatsPlot =
      ggplot(aggregateData, aes(Read_Length, IQR, color = Timepoint, shape = Position_Type)) +
      geom_jitter(height = 0, size = 3)
  } else if (stat == STDEV) {
    groupedStatsPlot =
      ggplot(aggregateData, aes(Read_Length, Standard_Deviation, color = Timepoint, shape = Position_Type)) +
      geom_jitter(height = 0, size = 3)
  } else {
    stop(paste("Unrecognized stat:",stat))
  }

  groupedStatsPlot = groupedStatsPlot + blankBackground + defaultTextScaling +
    labs(title = title, x = "Read Length", y = stat) +
    scale_x_continuous(breaks = xAxisBreaks)

  if (all(aggregateData$Timepoint == "None")) {
    groupedStatsPlot = groupedStatsPlot +
      scale_color_grey() + theme(legend.position = "none")
  } else groupedStatsPlot = groupedStatsPlot + scale_color_brewer(palette = "Set1")

  print(groupedStatsPlot)

}


# Plots the trinucleotide context for a set of mismatches
plotTrinucleotideContext = function(simplifiedMismatchTable, includedTypes = list(), omittedTypes = list(),
                                    title = "Trinucleotide Context") {

  if ( length(includedTypes) > 0 && length(omittedTypes) > 0) {
    stop("Included types and omitted types given simultaneously")
  } else if (length(includedTypes) > 0) {
    simplifiedMismatchTable = simplifiedMismatchTable[Mismatch %in% includedTypes]
  } else if (length(omittedTypes) > 0) {
    simplifiedMismatchTable = simplifiedMismatchTable[!(Mismatch %in% omittedTypes)]
  }

  trinucContextFrequencies = (simplifiedMismatchTable[!is.na(Trinuc_Context), .N, by = Trinuc_Context]
                              [, Freq := N/sum(N)])

  print(
    ggplot(trinucContextFrequencies, aes(Trinuc_Context, Freq, fill = str_sub(Trinuc_Context, 1, 1))) +
      geom_bar(stat = "identity") +
      labs(title = title, x = "Trinucleotide Context", y = "Frequency") +
      guides(x = guide_axis(angle = 45)) +
      theme(legend.position = "none", axis.text.x = element_text(size = 12)) +
      blankBackground + defaultTextScaling + scale_color_brewer(palette = "Set1")
  )

}