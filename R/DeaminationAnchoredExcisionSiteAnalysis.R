library(data.table)
library(stringr)
library(ggplot2)

THREE_PRIME = "three_prime"
FIVE_PRIME = "five_prime"

# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 22, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18, color = "black"),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Calculates nucleotide frequencies at individual positions in a set of reads (relative to the specified end).
# Takes optional arguments to add columns with additional, repeating information (e.g. read length, timepoint).
tabulateNucleotideFrequenciesByPosition = function(sequences, posType = THREE_PRIME,
                                                   paddingColNames = list(), paddingInfo = list()) {

  maxStringLength = max(nchar(sequences))
  if (posType == THREE_PRIME) {
    positions = -1:-maxStringLength
  } else if (posType == FIVE_PRIME) {
    positions = 1:maxStringLength
  } else stop("Unrecognized value for posType parameter.")

  nucFreqByPosition = data.table(Position = rep(positions, 4),
                                 Nucleotide = unlist(lapply(c('A','C','G','T'),
                                                            function(x) rep(x, maxStringLength))),
                                 Frequency = unlist(lapply(c('A','C','G','T'), function(x) {
                                   lapply(positions, function(y) {
                                     nucleotides = str_sub(sequences, y, y)
                                     return(sum(nucleotides == x) / sum(nucleotides != ''))
                                   })
                                 }))
  )

  mapply(function(x,y) nucFreqByPosition[, eval(x) := y], paddingColNames, paddingInfo)

  return(nucFreqByPosition)

}


tabulateNucFreqByPosTimepointAndLength = function(simplifiedTablesByTimepoint, posType = THREE_PRIME,
                                                  combineNucleotides = list(), combinedNames = list()) {

  #If passed a single data.table, wrap it in a list.
  if (is.data.table(simplifiedTablesByTimepoint)){
    simplifiedTablesByTimepoint = list(NONE = simplifiedTables)
  }
  noTimepoints = all(names(simplifiedTablesByTimepoint) == "NONE")

  # First, stratify the data by length and timepoint (optional),
  # calculating relative nucleotide frequencies within each category.
  nucFreqTable = rbindlist(lapply(seq_along(simplifiedTablesByTimepoint), function(i) {
    rbindlist(lapply(unique(simplifiedTablesByTimepoint[[i]]$Read_Length), function(x) {

      if (noTimepoints) {
        paddingInfo = list(x, "NONE")
      } else {
        paddingInfo = list(x, names(simplifiedTablesByTimepoint)[i])
      }

      tabulateNucleotideFrequenciesByPosition(simplifiedTablesByTimepoint[[i]][Read_Length == x]$Read_Sequence,
                                              posType, paddingColNames = c("Read_Length", "Timepoint"),
                                              paddingInfo = paddingInfo)

    }))
  }))

  # Next, perform any transformations on the frequencies as specified by the "combine" parameters.
  # For example, convert to counts of Purines and pyrimidines.
  if (length(combineNucleotides) > 0) {
    nucFreqTable = combineNucleotideFrequencies(nucFreqTable, combineNucleotides, combinedNames)
  }

  if (dim(nucFreqTable)[1] > 0) {
    return(nucFreqTable)
  } else {
    return(data.table(Position = numeric(), Nucleotide = character(),
                      Read_Length = numeric(), Timepoint = character()))
  }

}


# Calculates frequencies as a combination of individual nucleotide frequencies.
# For example, calculates purine frequencies by combining 'A' and 'G' nucleotides.
combineNucleotideFrequencies = function(frequencyData, combineNucleotides = list(), combinedNames = list()) {

  return(rbindlist(mapply(function(x,y) {
    frequencyData[sapply(Nucleotide, grepl, x),
                  .(Frequency = sum(Frequency), Nucleotide = y),
                  by = list(Position, Read_Length, Timepoint)]
  }, combineNucleotides, combinedNames, SIMPLIFY = FALSE)))

}


# Plot nucleotide frequencies for a series of read lengths.
plotNucFreqVsReadLengthBarPlot = function(nucFreqTable, posType = THREE_PRIME,
                                          title = "Nuc Freq by Length and Timepoint",
                                          yAxisLabel = "Nucleotide Frequency", secondaryYAxisLabel = "Read Length",
                                          yStripFontSize = 16, ylim = c(0,1.2), yAxisBreaks = c(0,0.5,1),
                                          combineNucleotides = list(), combinedNames = list(),
                                          showThreePrimeCutSite = FALSE, showFivePrimeCutSite = FALSE,
                                          expansionOffset = 0, xlim = NULL, xAxisBreaks = -3:3*10, xAxisLabel = NULL,
                                          minReadLength = NULL, maxReadLength = NULL, positionOffset = 0) {

  if (is.null(xAxisLabel)) {
    if (posType == THREE_PRIME) {
      xAxisLabel = "3' Relative Position"
    } else if (posType == FIVE_PRIME) {
      xAxisLabel = "5' Relative Position"
    } else stop("Unrecognized value for posType parameter.")
  }

  if (length(combineNucleotides) > 0) {
    nucFreqTable = combineNucleotideFrequencies(nucFreqTable, combineNucleotides, combinedNames)
  }

  combinedReadLengths = all(nucFreqTable$Read_Length == "combined")

  if (!combinedReadLengths) {
    if (!is.null(minReadLength)) nucFreqTable = nucFreqTable[Read_Length >= minReadLength]
    if (!is.null(maxReadLength)) nucFreqTable = nucFreqTable[Read_Length <= maxReadLength]
  }

  if (is.null(xlim)) {
    xlim = c(min(nucFreqTable$Position)-1 + positionOffset,
             max(nucFreqTable$Position)+1 + positionOffset)
  }

  plot = ggplot(nucFreqTable, aes(Position + positionOffset, Frequency, fill = Nucleotide)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    blankBackground + defaultTextScaling

  if (any(nucFreqTable$Timepoint != "NONE") && !combinedReadLengths) {
    plot = plot + facet_grid(Read_Length~factor(Timepoint, levels = unique(nucFreqTable$Timepoint)))
  } else if (any(nucFreqTable$Timepoint != "NONE") && combinedReadLengths) {
    plot = plot + facet_grid(cols = vars(factor(Timepoint, levels = unique(nucFreqTable$Timepoint))))
  } else if (all(nucFreqTable$Timepoint == "NONE") && !combinedReadLengths) {
    plot = plot + facet_grid(rows = vars(Read_Length))
  }

  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          strip.background = element_rect(color = "black", linewidth = 1),
          axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          axis.text.y.left = element_text(size = 10),
          strip.text.y = element_text(size = yStripFontSize)) +
    coord_cartesian(ylim = ylim, xlim = xlim, expand = FALSE) +
    scale_x_continuous(breaks = xAxisBreaks)

  if (combinedReadLengths) {
    plot = plot + scale_y_continuous(breaks = yAxisBreaks)
  } else {
    plot = plot + scale_y_continuous(sec.axis = dup_axis(~., name = secondaryYAxisLabel), breaks = yAxisBreaks)
  }

  if (length(unique(nucFreqTable$Nucleotide)) == 1) {
    plot = plot + theme(legend.position = "none") + scale_fill_manual(values = "grey35")
  } else {
    plot = plot + scale_fill_brewer(palette = "Set1")
  }

  if (showThreePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos + positionOffset),
                             nucFreqTable[,.(Cut_Site_Pos = max(Position) - expansionOffset + 0.5),
                                          by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }
  if (showFivePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos + positionOffset),
                             nucFreqTable[,.(Cut_Site_Pos = min(Position) + expansionOffset - 0.5),
                                          by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }

  print(plot)

}


# This function takes a table of mismatch data for whole reads and returns a table of
# that same data for one side of that data, anchored on the mismatch
orientDataToMismatch = function(simplifiedTable, posType, expansionOffset = 0, halfPositions = FALSE) {

  # Check to see if the input table is empty. If it is, return an empty table.
  if (dim(simplifiedTable)[1] == 0) return(data.table(Position = numeric(), Mismatch = character(),
                                                      Read_Length = numeric(), Read_Sequence = character()))

  if (posType == THREE_PRIME) {
    if (halfPositions) {
      sequences = str_sub(simplifiedTable$Read_Sequence,
                          -nchar(simplifiedTable$Read_Sequence), simplifiedTable$Position - 0.5)
    } else {
      sequences = str_sub(simplifiedTable$Read_Sequence,
                          -nchar(simplifiedTable$Read_Sequence), simplifiedTable$Position)
    }
    position = -1
  } else if (posType == FIVE_PRIME) {
    if (halfPositions) {
      sequences = str_sub(simplifiedTable$Read_Sequence, simplifiedTable$Position + 0.5, -1)
    } else {
      sequences = str_sub(simplifiedTable$Read_Sequence, simplifiedTable$Position, -1)
    }
    position = 1
  } else stop("Unrecognized value for posType parameter.")

  return(data.table(Position = position, Mismatch = simplifiedTable$Mismatch,
                    Read_Length = nchar(sequences)-expansionOffset, Read_Sequence = sequences))

}


# Given one or more nucleotide sequences, determine their frequency in a set of reads
# across position and read length (may need to be subset by timepoint).
# Optionally, a table of nucleotide frequencies may be provided to normalize by.
# If queryPattern is given (not null), querySequences is ignored, and query pattern is
# used in its place, interpreted as a regular expression (MUCH faster for more complex patterns).
# It is recommended that these patterns are all preceded by the '^' symbol avoid ambiguity on
# where the pattern started.
getSequenceFreqByReadLengthAndPos = function(readSequencesTable, expectedNucleotideFrequencies = NULL,
                                             querySequences = c("TGG"), posType = THREE_PRIME,
                                             queryPattern = NULL, maxQueryLength = NULL) {

  if (is.null(maxQueryLength)) maxQueryLength = max(nchar(querySequences))

  # Iterate through all the different read lengths present in the given table(s),
  # Calculating expected and observed frequencies across all the query sequences at each position
  # NOTE: In order to determine what positions to iterate through, the following code assumes that
  #       all read sequences under a specific read length actually have that length.
  seqFreqTable = rbindlist(lapply(unique(readSequencesTable$Read_Length), function(readLength) {

    sequences = readSequencesTable[Read_Length == readLength]$Read_Sequence
    if (!is.null(expectedNucleotideFrequencies)) {
      if (!is.null(queryPattern)) stop("Background calculation for queryPattern not implemented.")
      readLengthSpecificNucFreq = expectedNucleotideFrequencies[Read_Length == readLength]
    }
    if (posType == THREE_PRIME) {
      positions = -max(nchar(sequences)):(-maxQueryLength)
    } else if (posType == FIVE_PRIME) {
      positions = 1:(max(nchar(sequences)) - maxQueryLength + 1)
    } else stop("Unrecognized value for posType parameter.")

    frequencies = sapply(positions, function(position) {

      if (is.null(queryPattern)) {
        observedFrequency = sum(sapply(querySequences, function(querySequence) {
          sum(str_sub(sequences, position, position + nchar(querySequence) - 1) == querySequence)/length(sequences)
        }))
      } else {
        observedFrequency = sum(
          grepl(queryPattern, str_sub(sequences, position, position + maxQueryLength - 1))
        )/length(sequences)
      }

      if (!is.null(expectedNucleotideFrequencies)) {
        expectedFrequency = sum(sapply(querySequences, function(querySequence) {
          prod(sapply(1:nchar(querySequence), function(i) {
            readLengthSpecificNucFreq[Position == position + i - 1 &
                                        Nucleotide == str_sub(querySequence, i, i)]$Frequency
          }))
        }))

        # Add pseudo counts to avoid dividing by or taking the log (during plotting) of 0.
        pseudoCount = 1/length(sequences)
        observedFrequency = observedFrequency + pseudoCount
        expectedFrequency = expectedFrequency + pseudoCount

        return(observedFrequency/expectedFrequency)
      } else {
        return(observedFrequency)
      }

    })

    if (!is.null(expectedNucleotideFrequencies)) {
      return(data.table(Position = positions, Read_Length = readLength, Enrichment = frequencies))
    } else {
      return(data.table(Position = positions, Read_Length = readLength, Frequency = frequencies))
    }

  }))

  if (dim(seqFreqTable)[1] > 0) {
    return(seqFreqTable)
  } else {
    if (!is.null(expectedNucleotideFrequencies)) {
      return(data.table(Position = numeric(), Read_Length = numeric(), Enrichment = numeric()))
    } else {
      return(data.table(Position = numeric(), Read_Length = numeric(), Frequency = numeric()))
    }
  }

}


# Plots sequence enrichment values in a barplot, similar to nucleotide frequencies
plotSequenceEnrichment = function(enrichmentTablesByTimepoint, posType = THREE_PRIME,
                                  title = "Sequence Enrichment by Length and Timepoint",
                                  yAxisLabel = expression("log"[2]*"(Enrichment)"),
                                  secondaryYAxisLabel = "Read Length", yStripFontSize = 16, yAxisTickTextSize = 8,
                                  showThreePrimeCutSite = FALSE, showFivePrimeCutSite = FALSE,
                                  querySequences = c("TGG"), expansionOffset = 0, xAxisBreaks = -3:3*10,
                                  minReadLength = NULL, maxReadLength = NULL) {

  if (posType == THREE_PRIME) {
    xAxisLabel = "3' Relative Position"
  } else if (posType == FIVE_PRIME) {
    xAxisLabel = "5' Relative Position"
  } else stop("Unrecognized value for posType parameter.")

  # If passed a single data.table, no timepoint information is present.
  # Otherwise, combine the given tables into one agggregate table.
  if (is.data.table(enrichmentTablesByTimepoint)) {
    fullEnrichmentTable = copy(enrichmentTablesByTimepoint)[,Timepoint := "NONE"]
  } else {
    fullEnrichmentTable = rbindlist(lapply(seq_along(enrichmentTablesByTimepoint), function(i) {
      copy(enrichmentTablesByTimepoint[[i]])[,Timepoint := names(enrichmentTablesByTimepoint)[i]]
    }))
  }

  maxQueryLength = max(nchar(querySequences))

  if (!is.null(minReadLength)) fullEnrichmentTable = fullEnrichmentTable[Read_Length >= minReadLength]
  if (!is.null(maxReadLength)) fullEnrichmentTable = fullEnrichmentTable[Read_Length <= maxReadLength]

  plot = ggplot(fullEnrichmentTable, aes(Position, log2(Enrichment))) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    blankBackground + defaultTextScaling

  if (all(fullEnrichmentTable$Timepoint == "NONE")) {
    plot = plot + facet_grid(rows = vars(Read_Length))
  } else {
    plot = plot + facet_grid(Read_Length~factor(Timepoint, levels = unique(fullEnrichmentTable$Timepoint)))
  }

  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          strip.background = element_rect(color = "black", linewidth = 1),
          axis.text.y = element_text(size = yAxisTickTextSize),
          axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
          strip.text.y = element_text(size = yStripFontSize)) +
    #coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(sec.axis = dup_axis(~., name = secondaryYAxisLabel)) +
    scale_x_continuous(breaks = xAxisBreaks) + geom_hline(yintercept = 0)

  if (showThreePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos),
                             fullEnrichmentTable[,.(Cut_Site_Pos = max(Position) - expansionOffset +
                                                      maxQueryLength - 0.5),
                                                 by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }
  if (showFivePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos),
                             fullEnrichmentTable[,.(Cut_Site_Pos = min(Position) + expansionOffset - 0.5),
                                                 by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }

  print(plot)

}


# Plot nucleotide frequencies for a series of read lengths.
# Notes on statistics:
#   A background mean and standard deviation are calculated from a given number of bases at the
#   start (near 0) or end of the x-axis. The bar plot is colored based on whether the frequency at each position
#   meets a z-score threshold with respect to this background. The positions which are checked can be specified
#   using the same start- or end-based notation.
# Notes on xAxisOffset and positionOffset:
#   Both offsets are applied after statistics and serves purely to change the visuals on the graph.
#   The position values for the x-axis will be incremented by the given values, but cut-site positions will only be
#   changed by positionOffset, not xAxisOffset. (This is useful for adjusting which nucleotide in the sequence
#   is used to orient the bars in the plot.)
# Combining read lengths:
#   If all of the entries in the "Read_Length" columns are the string "combined", no secondary y-axis will be plotted.)
plotSequenceFrequencies = function(seqFreqTablesByTimepoint, posType = THREE_PRIME,
                                   title = "Seq Freq by Length and Timepoint", yAxisLabel = "Sequence Frequency",
                                   secondaryYAxisLabel = "Read Length", yStripFontSize = 16,
                                   showThreePrimeCutSite = FALSE, showFivePrimeCutSite = FALSE,
                                   expansionOffset = 0, querySequences = c("TGG"), xAxisBreaks = -3:3*10,
                                   minReadLength = NULL, maxReadLength = NULL, xlim = NULL, ylimMax = NULL, ylimMaxModifier = NULL,
                                   startBasedBackgroundNum = NULL, endBasedBackgroundNum = NULL, backgroundPositions = NULL,
                                   startBasedCheckPositions = NULL, endBasedCheckPositions = NULL, checkPositions = NULL,
                                   displayPValue = FALSE, zScoreCutoff = 4, significanceAsteriskBreakpoints = NULL,
                                   xAxisLabel = NULL, xAxisOffset = 0, positionOffset = 0, readLengthOffset = 0,
                                   defaultColor = "grey35") {

  if (is.null(xAxisLabel)) {
    if (posType == THREE_PRIME) {
      xAxisLabel = "3' Relative Position"
    } else if (posType == FIVE_PRIME) {
      xAxisLabel = "5' Relative Position"
    } else stop("Unrecognized value for posType parameter.")
  }

  if (sum(c(!is.null(startBasedBackgroundNum), !is.null(endBasedBackgroundNum)), !is.null(backgroundPositions)) > 1) {
    stop("More than 1 background position determinant given")
  }

  # If passed a single data.table, no timepoint information is present.
  # Otherwise, combine the given tables into one agggregate table.
  if (is.data.table(seqFreqTablesByTimepoint)) {
    fullFrequencyTable = copy(seqFreqTablesByTimepoint)[,Timepoint := "NONE"]
  } else {
    fullFrequencyTable = rbindlist(lapply(seq_along(seqFreqTablesByTimepoint), function(i) {
      copy(seqFreqTablesByTimepoint[[i]])[,Timepoint := names(seqFreqTablesByTimepoint)[i]]
    }))
  }

  combinedReadLengths = all(fullFrequencyTable$Read_Length == "combined")

  maxFrequency = signif(max(fullFrequencyTable$Frequency),digits = 2)
  if (is.null(ylimMax)) {
    yAxisBreaks = c(0, maxFrequency/2, maxFrequency)
  } else {
    yAxisBreaks = c(0, ylimMax/2, ylimMax)
  }

  maxQueryLength = max(nchar(querySequences))

  if (!combinedReadLengths) {
    fullFrequencyTable[,Read_Length := Read_Length + readLengthOffset]
    if (!is.null(minReadLength)) fullFrequencyTable = fullFrequencyTable[Read_Length >= minReadLength]
    if (!is.null(maxReadLength)) fullFrequencyTable = fullFrequencyTable[Read_Length <= maxReadLength]
  }

  if (combinedReadLengths) fullFrequencyTable[,Read_Length := "combined"]

  # Determine zScores from background positions, if specified.
  if (!is.null(startBasedBackgroundNum) || !is.null(endBasedBackgroundNum) || !is.null(backgroundPositions)) {

    if (!is.null(startBasedBackgroundNum)) {
      backgroundRows = fullFrequencyTable[abs(Position) <= abs(startBasedBackgroundNum)]
    } else if (!is.null(endBasedBackgroundNum)) {
      endPositionsByTimepointAndRL = fullFrequencyTable[, .(End = max(abs(Position))),
                                                        by = list(Timepoint, Read_Length)]
      setkey(endPositionsByTimepointAndRL, Timepoint, Read_Length)
      endPosByFullFreqRow = endPositionsByTimepointAndRL[list(fullFrequencyTable$Timepoint,
                                                              fullFrequencyTable$Read_Length)]$End
      backgroundRows = fullFrequencyTable[abs(Position) + abs(endBasedBackgroundNum) > endPosByFullFreqRow]
    } else if (!is.null(backgroundPositions)) {
      backgroundRows = fullFrequencyTable[Position %in% backgroundPositions]
    }

    backgroundStats = backgroundRows[,.(Mean = mean(Frequency), SD = sd(Frequency)), by = list(Timepoint, Read_Length)]
    setkey(backgroundStats, Timepoint, Read_Length)
    statsByFullFreqRow = backgroundStats[list(fullFrequencyTable$Timepoint, fullFrequencyTable$Read_Length)]
    fullFrequencyTable[,Z_Score := (Frequency-statsByFullFreqRow$Mean) / statsByFullFreqRow$SD]
    fullFrequencyTable[,P_Value := 2*pnorm(abs(Z_Score), lower.tail = FALSE)]
    fullFrequencyTable[,Significant := abs(Z_Score) > zScoreCutoff]

    if (!is.null(startBasedCheckPositions) || !is.null(endBasedCheckPositions) || !is.null(checkPositions)) {
      fullFrequencyTable[,Relevant := FALSE]
      if (!is.null(startBasedCheckPositions)) {
        fullFrequencyTable[abs(Position) %in% abs(startBasedCheckPositions), Relevant := TRUE]
      }

      if (!is.null(endBasedCheckPositions)) {
        endPositionsByTimepointAndRL = fullFrequencyTable[, .(End = max(abs(Position))),
                                                          by = list(Timepoint, Read_Length)]
        setkey(endPositionsByTimepointAndRL, Timepoint, Read_Length)
        endPosByFullFreqRow = endPositionsByTimepointAndRL[list(fullFrequencyTable$Timepoint,
                                                                fullFrequencyTable$Read_Length)]$End
        fullFrequencyTable[(endPosByFullFreqRow - abs(Position) + 1) %in% abs(endBasedCheckPositions), Relevant := TRUE]
      }

      if (!is.null(checkPositions)) {
        fullFrequencyTable[Position %in% checkPositions, Relevant := TRUE]
      }

    } else {
      fullFrequencyTable[,Relevant := TRUE]
    }
  }

  if (is.null(xlim)) {
    xlim = c(min(fullFrequencyTable$Position)-1 + xAxisOffset + positionOffset,
             max(fullFrequencyTable$Position)+1 + xAxisOffset + positionOffset)
  }

  plot = ggplot(fullFrequencyTable, aes(Position + xAxisOffset + positionOffset, Frequency)) +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    blankBackground + defaultTextScaling

  if (!is.null(startBasedBackgroundNum) || !is.null(endBasedBackgroundNum) || !is.null(backgroundPositions)) {
    plot = plot +
      geom_bar(aes(fill = Significant & Relevant), stat = "identity") +
      scale_fill_manual(values = c("TRUE" = "red", "FALSE" = defaultColor), guide = "none")
  } else {
    plot = plot + geom_bar(stat = "identity", fill = defaultColor)
  }

  if (displayPValue) {
    plot = plot + geom_text(aes(label = ifelse(Relevant & Significant, paste0("p=",format(signif(P_Value, 3))), '')),
                            size = 3, hjust = 1, vjust = 0, nudge_x = -0.05, nudge_y = 0.05)
  } else if (!is.null(significanceAsteriskBreakpoints)) {
    fullFrequencyTable[,Asterisks := ""]
    for (i in seq_along(significanceAsteriskBreakpoints)) {
      fullFrequencyTable[Relevant & P_Value < significanceAsteriskBreakpoints[i],
                         Asterisks := paste0(rep('*',i), collapse = '')]
    }
    plot = plot + geom_text(aes(label = Asterisks), size = 4, hjust = 0.5, vjust = 0, nudge_x = 0, nudge_y = -0.05)
  }

  if (any(fullFrequencyTable$Timepoint != "NONE") && !combinedReadLengths) {
    plot = plot + facet_grid(Read_Length~factor(Timepoint, levels = unique(fullFrequencyTable$Timepoint)))
  } else if (any(fullFrequencyTable$Timepoint != "NONE") && combinedReadLengths) {
    plot = plot + facet_grid(cols = vars(factor(Timepoint, levels = unique(fullFrequencyTable$Timepoint))))
  } else if (all(fullFrequencyTable$Timepoint == "NONE") && !combinedReadLengths) {
    plot = plot + facet_grid(rows = vars(Read_Length))
  }

  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          strip.background = element_rect(color = "black", linewidth = 1),
          axis.text.y = element_text(size = 10),
          axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          strip.text.y = element_text(size = yStripFontSize)) +
    scale_x_continuous(breaks = xAxisBreaks)

  if (combinedReadLengths) {
    plot = plot + scale_y_continuous(breaks = yAxisBreaks)
  } else {
    plot = plot + scale_y_continuous(sec.axis = dup_axis(~., name = secondaryYAxisLabel), breaks = yAxisBreaks)
  }

  if (is.null(ylimMax)) ylimMax = maxFrequency

  if (is.null(ylimMaxModifier)) {
    if (displayPValue) {
      ylimMaxModifier = 1.4
    } else {
      ylimMaxModifier = 1.2
    }
  }

  plot = plot + coord_cartesian(ylim = c(0, ylimMax*ylimMaxModifier), xlim = xlim, expand = FALSE)

  if (showThreePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos + positionOffset),
                             fullFrequencyTable[,.(Cut_Site_Pos = max(Position) - expansionOffset +
                                                     maxQueryLength - 0.5),
                                                by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }
  if (showFivePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos + positionOffset),
                             fullFrequencyTable[,.(Cut_Site_Pos = min(Position) + expansionOffset - 0.5),
                                                by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }

  print(plot)

}
