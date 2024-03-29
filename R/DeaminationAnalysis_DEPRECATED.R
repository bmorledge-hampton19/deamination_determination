library(data.table)
library(stringr)
library(ggplot2)

THREE_PRIME = 3
FIVE_PRIME = 5

# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 26, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Filters the given mismatch data based the presence of 'N', the number of mismatches per read,
# and read length.
filterResults = function(mismatchData, removeRowsWithN = TRUE, maxMismatchesAllowed = NA,
                         minReadLength = NA, maxReadLength = NA) {

  if (removeRowsWithN) mismatchData = mismatchData[!grepl('N', V5)]

  if (!is.na(maxMismatchesAllowed)) {
    mismatchData = mismatchData[str_count(V4, ':') < maxMismatchesAllowed]
  }

  if (!is.na(minReadLength)) {
    mismatchData = mismatchData[V3-V2 >= minReadLength]
  }

  if (!is.na(maxReadLength)) {
    mismatchData = mismatchData[V3-V2 <= maxReadLength]
  }

  return(mismatchData)

}


# Simplifies the given table of bed formatted mismatches into a smaller
# three-column table of positions, sequence context, and read length and optionally,
# read sequence and trinuc context.
# The expansionOffset parameter is used to change the Position and Read_Length values
# relative to the number of expanded bases for the read. As a result, Read_Length should still reflect the
# original read length, but Position will change to reflect the mismatch's position in the new sequence
# NOTE: In the event that the read contains expanded sequence context information, that information
#       no longer reflects the read sequence and instead reflects an expansion of the corresponding reference
#       genome sequence. This also means that if the read sequence contained any indels, the Read_Length field will
#       be inaccurate and the Position field MAY be inaccurate, even with the expansionOffset parameter.
simplifyTable = function(table, includeReadSequence = FALSE, includeTrinucleotideContext = FALSE,
                         expansionOffset = 0) {

  returnTable = table[,list(Position = as.numeric(unlist(strsplit(V4, ':'))) - expansionOffset,
                            Mismatch = unlist(strsplit(V5, ':')),
                            Read_Length = unlist(mapply(rep, nchar(V7)-2*expansionOffset,
                                                        lapply(strsplit(V5,':'), length))))]

  if (includeReadSequence || includeTrinucleotideContext) {
    returnTable[, Read_Sequence := unlist(mapply(rep, table$V7, lapply(strsplit(table$V5,':'), length)))]
  }

  if (includeTrinucleotideContext) {
    returnTable[abs(Position) > 1 & abs(Position) < Read_Length,
                Trinuc_Context := paste0(str_sub(Read_Sequence, Position - 1, Position - 1),
                                         str_sub(Mismatch, 1, 1),
                                         str_sub(Read_Sequence, Position + 1, Position + 1))]
    if (!includeReadSequence) returnTable[, Read_Sequence := NULL]
  }

  return(returnTable)

}


# Takes a simplified table of mismatches and filters them based on the position of the mismatch
filterMismatchesByPosition = function(mismatchTable, minPos, maxPos, posType, expansionOffset = 0) {

  if (posType == THREE_PRIME) {
    positions = mismatchTable$Position + expansionOffset
  } else if (posType == FIVE_PRIME) {
    positions = mismatchTable$Read_Length + mismatchTable$Position + 1 - expansionOffset
  } else stop("Unrecognized value for posType parameter.")

  return(mismatchTable[positions >= minPos & positions <= maxPos])

}


# Takes a simplified table of mismatches and filters them based on the position of the mismatch.
# Unlike the filterMismatchesByPosition function, min and max positions are dynamic based on read length.
# (This function is actually still used under the hood.)
# See constant tables below for examples of how to format the posConstraintsByReadLength table.
filterMismatchesByPositionAndReadLength = function(mismatchTable, posConstraintsByReadLength,
                                                   expansionOffset = 0) {

  return(rbindlist(mapply( function(readLength, minPos, maxPos) {
    filterMismatchesByPosition(mismatchTable[Read_Length == readLength],
                               minPos, maxPos, THREE_PRIME, expansionOffset)
  },
  posConstraintsByReadLength$Read_Length,
  posConstraintsByReadLength$Min_Pos,
  posConstraintsByReadLength$Max_Pos,
  SIMPLIFY = FALSE)))

}


HUMAN_THREE_PRIME_POS_CONSTRAINTS = data.table(Read_Length = c(23, 24, 25, 26, 27, 28, 29, 30, 31),
                                               Max_Pos = c(-3, -3, -4, -5, -5, -6, -6, -6, -6),
                                               Min_Pos = c(-9, -9, -10, -10, -10, -11, -11, -12, -12))
# HUMAN_FIVE_PRIME_POS_CONSTRAINTS = copy(HUMAN_THREE_PRIME_POS_CONSTRAINTS)
# HUMAN_FIVE_PRIME_POS_CONSTRAINTS[, Max_Pos := Read_Length + Max_Pos + 1]
# HUMAN_FIVE_PRIME_POS_CONSTRAINTS[, Min_Pos := Read_Length + Min_Pos + 1]


# The following functions get the first or last quartile of a data.table sorted by position.
getFirstPositionQuartile = function(table) {
  return(table[order(Position)][1:(dim(table)[1]/4)])
}
getLastPositionQuartile = function(table) {
  return(table[order(Position)][(dim(table)[1]*3/4):dim(table)[1]])
}


# Determines if the given mismatch position and type data for a single read
# contains a tandem C>T mutation.
hasTandemCTMismatch = function(mismatchPositions, mismatchTypes) {

  if (length(mismatchPositions) < 2) return(NA)

  for (i in 2:length(mismatchPositions)) {
    if (abs(mismatchPositions[i] - mismatchPositions[i-1]) == 1 &&
        mismatchTypes[i] == "C>T" && mismatchTypes[i-1] == "C>T") return(TRUE)
  }

  return(FALSE)

}


# Finds adjacent C>T mismatches in the given table of reads.
# Can also filter results by boundary conditions (e.g. a threePrimeBoundary of -4 will filter out
# tandem mismatches that extend into the -3, -2, or -1 positions.)
getTandemCTMismatches = function(table, threePrimeBoundary = NA, fivePrimeBoundary = NA) {

  multiMismatchReads = table[grepl(':', V4)]
  tandemCTMismatches = multiMismatchReads[mapply(hasTandemCTMismatch,
                                                 lapply(V4, function(x) as.numeric(unlist(strsplit(x,':')))),
                                                 strsplit(V5, ':'))]
  if (!is.na(threePrimeBoundary)) {
    tandemCTMismatches = tandemCTMismatches[sapply(strsplit(V4, ':'),
                                                   function(x) min(as.numeric(x)) < threePrimeBoundary)]
  }
  if (!is.na(fivePrimeBoundary)) {
    tandemCTMismatches = tandemCTMismatches[sapply(strsplit(V4, ':'),
                                                   function(x) min(as.numeric(x)) > fivePrimeBoundary)]
  }
  return(tandemCTMismatches)

}

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
                                                    title = "Mismatch Position Frequencies", posType = THREE_PRIME) {

  #If passed a single data.table, wrap it in a list.
  if (is.data.table(simplifiedTablesByTimepoint)) {
    simplifiedTablesByTimepoint = list(None = simplifiedTablesByTimepoint)
    noTimepointInfo = TRUE
  } else noTimepointInfo = FALSE

  if (posType == THREE_PRIME) {
    xAxisLabel = "3' Relative Position"
    xAxisBreaks = c(0, -10, -20)
  } else if (posType == FIVE_PRIME) {
    xAxisLabel = "5' Relative Position"
    xAxisBreaks = c(0, 10, 20)
    simplifiedTables = lapply(simplifiedTables, copy)
    lapply(simplifiedTables, function(x) x[, Position := Read_Length + Position + 1])
  } else stop("Unrecognized value for posType parameter.")

  aggregateTable = rbindlist(lapply(seq_along(simplifiedTables),
                                    function(i) simplifiedTables[[i]][,Timepoint := names(simplifiedTables)[i]]))

  if ( length(includedTypes) > 0 && length(omittedTypes) > 0) {
    stop("Included types and omitted types given simultaneously")
  } else if (length(includedTypes) > 0) {
    aggregateTable = aggregateTable[Mismatch %in% includedTypes]
  } else if (length(omittedTypes) > 0) {
    aggregateTable = aggregateTable[!(Mismatch %in% omittedTypes)]
  }

  groupedPositionFrequencies = (aggregateTable[, .N, by = list(Position,Read_Length,Timepoint)]
                                [, Freq := N/sum(N), by = list(Read_Length, Timepoint)])

  plot = ggplot(groupedPositionFrequencies, aes(Position, Freq)) +
    geom_bar(stat = "identity") +
    labs(title = title, x = xAxisLabel, y = "Relative Mismatch Frequency") +
    blankBackground + defaultTextScaling

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
          panel.grid.major.x = element_line(color = "red", size = 0.5, linetype = 2)) +
    scale_y_continuous(sec.axis = dup_axis(~., name = "Read Length")) +
    scale_x_continuous(breaks = xAxisBreaks)

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
  if (is.data.table(simplifiedTablesByTimepoint)) {
    simplifiedTablesByTimepoint = list(None = simplifiedTablesByTimepoint)
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

  if (all(aggregateData$Timepoint == None)) {
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
    simplifiedTablesByTimepoint = list(simplifiedTablesByTimepoint)
    noTimepoints = TRUE
  } else noTimepoints = FALSE

  # First, stratify the data by length and timepoint (optional),
  # calculating relative nucleotide frequencies within each category.
  nucFreqTable = rbindlist(lapply(seq_along(simplifiedTablesByTimepoint), function(i) {
    rbindlist(lapply(unique(simplifiedTablesByTimepoint[[i]]$Read_Length), function(x) {

      if (noTimepoints) {
        paddingInfo = list(x, NA)
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

  return(nucFreqTable)

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
                                          yStripFontSize = 16,
                                          combineNucleotides = list(), combinedNames = list(),
                                          showThreePrimeCutSite = FALSE, showFivePrimeCutSite = FALSE,
                                          expansionOffset = 0) {

  if (posType == THREE_PRIME) {
    xAxisLabel = "3' Relative Position"
    xAxisBreaks = c(0, -10, -20, -30)
  } else if (posType == FIVE_PRIME) {
    xAxisLabel = "5' Relative Position"
    xAxisBreaks = c(0, 10, 20, 30)
  } else stop("Unrecognized value for posType parameter.")

  if (length(combineNucleotides) > 0) {
    nucFreqTable = combineNucleotideFrequencies(nucFreqTable, combineNucleotides, combinedNames)
  }

  plot = ggplot(nucFreqTable, aes(Position, Frequency, fill = Nucleotide)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    blankBackground + defaultTextScaling

  if (all(is.na(nucFreqTable$Timepoint))) {
    plot = plot + facet_grid(rows = vars(Read_Length))
  } else {
    plot = plot + facet_grid(Read_Length~factor(Timepoint, levels = unique(nucFreqTable$Timepoint)))
  }

  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          strip.background = element_rect(color = "black", size = 1),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y = element_text(size = yStripFontSize)) +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(sec.axis = dup_axis(~., name = secondaryYAxisLabel)) +
    scale_x_continuous(breaks = xAxisBreaks)

  if (length(unique(nucFreqTable$Nucleotide)) == 1) {
    plot = plot + theme(legend.position = "none") + scale_fill_grey()
  } else {
    plot = plot + scale_fill_brewer(palette = "Set1")
  }

  if (showThreePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos),
                             nucFreqTable[,.(Cut_Site_Pos = max(Position) - expansionOffset + 0.5),
                                          by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }
  if (showFivePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos),
                             nucFreqTable[,.(Cut_Site_Pos = min(Position) + expansionOffset - 0.5),
                                          by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }

  print(plot)

}


# This function takes a table of mismatch data for whole reads and returns a table of
# that same data for one side of that data, anchored on the mismatch
orientDataToMismatch = function(simplifiedTable, posType, expansionOffset = 0) {

  if (posType == THREE_PRIME) {
    sequences = str_sub(simplifiedTable$Read_Sequence,
                        -nchar(simplifiedTable$Read_Sequence), simplifiedTable$Position)
    position = -1
  } else if (posType == FIVE_PRIME) {
    sequences = str_sub(simplifiedTable$Read_Sequence, simplifiedTable$Position, -1)
    position = 1
  } else stop("Unrecognized value for posType parameter.")

  return(data.table(Position = position, Mismatch = simplifiedTable$Mismatch,
                    Read_Length = nchar(sequences)-expansionOffset, Read_Sequence = sequences))

}


# Given one or more nucleotide sequences, determine their frequency in a set of reads
# across position and read length (may need to be subset by timepoint).
# Optionally, a table of nucleotide frequencies may be provided to normalize by.
getSequenceFreqByReadLengthAndPos = function(readSequencesTable, expectedNucleotideFrequencies = NULL,
                                             querySequences = c("TGG"), posType = THREE_PRIME) {

  maxQueryLength = max(nchar(querySequences))

  # Iterate through all the different read lengths present in the given table(s),
  # Calculating expected and observed frequencies across all the query sequences at each position
  # NOTE: In order to determine what positions to iterate through, the following code assumes that
  #       all read sequences under a specific read length actually have that length.
  return(rbindlist(lapply(unique(readSequencesTable$Read_Length), function(readLength) {

    sequences = readSequencesTable[Read_Length == readLength]$Read_Sequence
    if (!is.null(expectedNucleotideFrequencies)) {
      readLengthSpecificNucFreq = expectedNucleotideFrequencies[Read_Length == readLength]
    }
    if (posType == THREE_PRIME) {
      positions = -nchar(sequences[1]):(-1 - maxQueryLength + 1)
    } else if (posType == FIVE_PRIME) {
      positions = 1:(nchar(sequences[1]) - maxQueryLength + 1)
    } else stop("Unrecognized value for posType parameter.")

    frequencies = sapply(positions, function(position) {

      observedFrequency = sum(sapply(querySequences, function(querySequence) {
        sum(str_sub(sequences, position, position + nchar(querySequence) - 1) == querySequence)/length(sequences)
      }))

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

  })))

}


# Plots sequence enrichment values in a barplot, similar to nucleotide frequencies
plotSequenceEnrichment = function(enrichmentTablesByTimepoint, posType = THREE_PRIME,
                                  title = "Sequence Enrichment by Length and Timepoint",
                                  yAxisLabel = expression("log"[2]*"(Enrichment)"),
                                  secondaryYAxisLabel = "Read Length", yStripFontSize = 16, yAxisTickTextSize = 8,
                                  showThreePrimeCutSite = FALSE, showFivePrimeCutSite = FALSE,
                                  querySequences = c("TGG"), expansionOffset = 0) {

  if (posType == THREE_PRIME) {
    xAxisLabel = "3' Relative Position"
    xAxisBreaks = c(0, -10, -20, -30, -40)
  } else if (posType == FIVE_PRIME) {
    xAxisLabel = "5' Relative Position"
    xAxisBreaks = c(0, 10, 20, 30, 40)
  } else stop("Unrecognized value for posType parameter.")

  # If passed a single data.table, no timepoint information is present.
  # Otherwise, combine the given tables into one agggregate table.
  if (is.data.table(enrichmentTablesByTimepoint)) {
    fullEnrichmentTable = copy(enrichmentTablesByTimepoint)[,Timepoint := NA]
  } else {
    fullEnrichmentTable = rbindlist(lapply(seq_along(enrichmentTablesByTimepoint), function(i) {
      copy(enrichmentTablesByTimepoint[[i]])[,Timepoint := names(enrichmentTablesByTimepoint)[i]]
    }))
  }

  maxQueryLength = max(nchar(querySequences))

  plot = ggplot(fullEnrichmentTable, aes(Position, log2(Enrichment))) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    blankBackground + defaultTextScaling

  if (all(is.na(fullEnrichmentTable$Timepoint))) {
    plot = plot + facet_grid(rows = vars(Read_Length))
  } else {
    plot = plot + facet_grid(Read_Length~factor(Timepoint, levels = unique(fullEnrichmentTable$Timepoint)))
  }

  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          strip.background = element_rect(color = "black", size = 1),
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
plotSequenceFrequencies = function(seqFreqTablesByTimepoint, posType = THREE_PRIME,
                                   title = "Seq Freq by Length and Timepoint", yAxisLabel = "Sequence Frequency",
                                   secondaryYAxisLabel = "Read Length", yStripFontSize = 16,
                                   showThreePrimeCutSite = FALSE, showFivePrimeCutSite = FALSE,
                                   expansionOffset = 0, querySequences = c("TGG")) {

  if (posType == THREE_PRIME) {
    xAxisLabel = "3' Relative Position"
    xAxisBreaks = c(0, -10, -20, -30)
  } else if (posType == FIVE_PRIME) {
    xAxisLabel = "5' Relative Position"
    xAxisBreaks = c(0, 10, 20, 30)
  } else stop("Unrecognized value for posType parameter.")

  # If passed a single data.table, no timepoint information is present.
  # Otherwise, combine the given tables into one agggregate table.
  if (is.data.table(seqFreqTablesByTimepoint)) {
    fullFrequencyTable = copy(seqFreqTablesByTimepoint)[,Timepoint := NA]
  } else {
    fullFrequencyTable = rbindlist(lapply(seq_along(seqFreqTablesByTimepoint), function(i) {
      copy(seqFreqTablesByTimepoint[[i]])[,Timepoint := names(seqFreqTablesByTimepoint)[i]]
    }))
  }

  maxQueryLength = max(nchar(querySequences))

  plot = ggplot(fullFrequencyTable, aes(Position, Frequency)) +
    geom_bar(stat = "identity") +
    labs(title = title, x = xAxisLabel, y = yAxisLabel) +
    blankBackground + defaultTextScaling

  if (all(is.na(fullFrequencyTable$Timepoint))) {
    plot = plot + facet_grid(rows = vars(Read_Length))
  } else {
    plot = plot + facet_grid(Read_Length~factor(Timepoint, levels = unique(fullFrequencyTable$Timepoint)))
  }

  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          strip.background = element_rect(color = "black", size = 1),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y = element_text(size = yStripFontSize)) +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(sec.axis = dup_axis(~., name = secondaryYAxisLabel)) +
    scale_x_continuous(breaks = xAxisBreaks)

  if (showThreePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos),
                             fullFrequencyTable[,.(Cut_Site_Pos = max(Position) - expansionOffset +
                                                                  maxQueryLength - 0.5),
                                                by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }
  if (showFivePrimeCutSite) {
    plot = plot + geom_vline(aes(xintercept = Cut_Site_Pos),
                             fullFrequencyTable[,.(Cut_Site_Pos = min(Position) + expansionOffset - 0.5),
                                                by = list(Read_Length, Timepoint)],
                             linetype = "dashed")
  }

  print(plot)

}


# Plot a scatter and accompanying trend line for the frequency of given nucleotides vs. read length.
plotNucFreqVsReadLengthLinearRegression = function(data, columns, title = "Nucleotides Vs. Read Length") {

  plotData = data.table(Nuc_Freq = data[,rowSums(.SD), .SDcols = columns] /
                          data[,rowSums(.SD), .SDcols = c("A_Counts", "C_Counts", "G_Counts", "T_Counts")],
                        Read_Length = data$Read_Length)

  linearRegressionSummary = summary(lm(plotData$Read_Length ~ plotData$Nuc_Freq))

  statsText = substitute(atop(italic(y) == m*italic(x) + b, italic(r)^2~"="~r2),
                         list(m = format(round(linearRegressionSummary$coefficients[2], 2), nsmall = 2),
                              b = format(round(linearRegressionSummary$coefficients[1], 2), nsmall = 2),
                              r2 = format(round(linearRegressionSummary$r.squared, 2), nsmall = 2)))
  statsText = as.expression(statsText)

  print(
    ggplot(data, aes(GC_Content, Read_Length)) +
      geom_point(size = 2) +
      geom_smooth(method = lm, se = FALSE) +
      annotate("text", label = "stats go here", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size = 18/.pt) +
      labs(title = title, x = "GC Content", y = "Read Length") +
      blankBackground + defaultTextScaling
  )

}


# Testing resized points and displaying expression text
testFun = function(data, myExpression, title = "testing") {
  test = CPDReadLengthAndGCData$rep1$hr1[, .(Size = log(.N) + 1), by = list(Read_Length, GC_Content)]
  print(
    ggplot(data, aes(test$GC_Content, test$Read_Length, size = test$Size)) +
      geom_point() +
      geom_smooth(method = lm, se = FALSE) +
      annotate("text", label = myExpression,
               x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size = 18/.pt) +
      labs(title = title, x = "GC Content", y = "Read Length") +
      blankBackground + defaultTextScaling
  )
}
