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

  returnTable = table[,list(Position = as.numeric(unlist(strsplit(as.character(V4), ':'))) - expansionOffset,
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


# Takes a simplified table of mismatches and filters them based on the position of the mismatch.
# Valid positions are defined using a zscore table that is stratified by both position and read length.
# NOTE: Position is assumed to be three-prime oriented (negative position values)
filterMismatchesByPositionAndReadLength = function(mismatchTable, posZScoreByReadLength,
                                                   zScoreCutoff = 4, expansionOffset = 0) {

  return(rbindlist(lapply(unique(posZScoreByReadLength$Read_Length), function(readLength) {
    relevantMismatches = mismatchTable[Read_Length == readLength]
    validPositions = posZScoreByReadLength[Read_Length == readLength & Z_Score >= zScoreCutoff, Position]
    return(relevantMismatches[(Position + expansionOffset) %in% validPositions])
  })))

}


# This function gets a z-score for mismatch frequency at every read position, stratified by read length.
# Z scores are calculated relative to a background defined as the first 10 positions on the 5' end (which
# should just be noise).
# NOTE: Positions are assumed to be given relative to the three prime end (negative values)
getMismatchFrequencyZScoreByPosAndReadLength = function(simplifiedMismatchData, expansionOffset,
                                                        noisePositions = 10) {

  # First, get position frequencies for each read length.
  mismatchPositionFrequencies = simplifiedMismatchData[, .N, by = list(Position,Read_Length)]
  mismatchPositionFrequencies[, Frequency := N/sum(N), by = list(Read_Length)]
  mismatchPositionFrequencies[,Position := Position + expansionOffset] # Adjust position for expansion offset

  # Next, iterate through read lengths, determining peak regions for each.
  return(rbindlist(lapply(unique(mismatchPositionFrequencies$Read_Length), function(readLength) {

    # Determine the cutoff value based on the background.
    relevantPositionFrequencies = mismatchPositionFrequencies[Read_Length == readLength]
    background = relevantPositionFrequencies[Position >= -readLength & Position < -readLength + noisePositions,
                                             Frequency]
    backgroundMean = mean(background)
    backgroundSD = sd(background)

    return(data.table(Read_Length = readLength,
                      Position = relevantPositionFrequencies$Position,
                      Frequency = relevantPositionFrequencies$Frequency,
                      Z_Score = (relevantPositionFrequencies$Frequency - backgroundMean) / backgroundSD))

  })))

}


# Plot mismatch frequency z-score using a facet plot with timepoint (optional) on one dimension and read length on the other
# Input should be a list of data.tables with names as timepoint information, or a single data.table
# to construct a plot without timepoint information.
# DON'T USE THIS IT LOOKS BAD :(
plotZScoreAcrossTimepointAndReadLength = function(zScoreTables, title = "Mismatch Position Frequencies",
                                                  zScoreCutoff = 4) {

  #If passed a single data.table, wrap it in a list.
  if (is.data.table(zScoreTables)) {
    zScoreTables = list(NONE = zScoreTables)
  }
  noTimepointInfo = all(names(zScoreTables) == "NONE")

  xAxisLabel = "3' Relative Position"
  xAxisBreaks = c(0, -10, -20, -30)

  aggregateZScoreTable = rbindlist(lapply(seq_along(zScoreTables),
                                          function(i) zScoreTables[[i]][,Timepoint := names(zScoreTables)[i]]))
  # aggregateZScoreTable[,Log_10_Z_Score := log10(Z_Score)]
  # maxLog10ZScore = round(max(aggregateZScoreTable$Log_10_Z_Score), digits = 2)
  # yAxisBreaks = c(0, maxLog10ZScore/2, maxLog10ZScore)
  # zScoreCutoff = log10(zScoreCutoff)
  maxZScore = round(max(aggregateZScoreTable$Z_Score), digits = 2)
  yAxisBreaks = c(0, maxZScore/2, maxZScore)

  plot = ggplot(aggregateZScoreTable, aes(Position, Z_Score, fill = Z_Score >= zScoreCutoff)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
    labs(title = title, x = xAxisLabel, y = expression("log"[10]*"(Mismatch Frequency Z-Score)")) +
    blankBackground + defaultTextScaling

  if (noTimepointInfo) {
    plot = plot + facet_grid(rows = vars(Read_Length)) +
      scale_y_continuous(breaks = yAxisBreaks)
  } else {
    plot = plot + facet_grid(Read_Length~factor(Timepoint, levels = names(zScoreTables))) +
      scale_y_continuous(sec.axis = dup_axis(~., name = "Read Length"), breaks = yAxisBreaks)
  }

  plot = plot +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          strip.background = element_rect(color = "black", size = 1),
          axis.text.y = element_text(size = 12), strip.text.y = element_text(size = 16),
          axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank()) +
    scale_x_continuous(breaks = xAxisBreaks) +
    geom_hline(yintercept = zScoreCutoff, linetype = 2, size = 0.5)

  print(plot)

}


# HUMAN_THREE_PRIME_POS_CONSTRAINTS = data.table(Read_Length = c(23, 24, 25, 26, 27, 28, 29, 30, 31),
#                                                Max_Pos = c(-3, -3, -4, -5, -5, -6, -6, -6, -6),
#                                                Min_Pos = c(-9, -9, -10, -10, -10, -11, -11, -12, -12))
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
