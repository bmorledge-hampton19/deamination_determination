library(data.table)
library(stringr)
library(ggplot2)

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
# three-column table of positions, sequence context, and read length.
simplifyTable = function(table) {
  return(table[,list(Position = as.numeric(unlist(strsplit(V4, ':'))),
                     Mismatch = unlist(strsplit(V5, ':')),
                     Read_Length = unlist(mapply(rep, V3-V2, lapply(strsplit(V5,':'), length))))])
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
                                           title = "Position Frequency") {

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
      labs(title = title, x = "3` Relative Position", y = "Frequency") +
      blankBackground + defaultTextScaling
  )

}
