library(data.table)
library(stringr)
library(ggplot)

# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 26, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Filters the given mismatch data based the presence of 'N' and the number
filterResults = function(mismatchData, removeRowsWithN = TRUE, maxMismatchesAllowed = 3) {
  if (removeRowsWithN) mismatchData = mismatchData[!grepl('N', V5)]
  mismatchData = mismatchData[str_count(V4, ':') <= maxMismatchesAllowed]
}

# Condenses the given table of bed formatted mismatches into a smaller
# two-column table of positions and sequence context.
condenseTable = function(table) {
  return(table[,list(Position = unlist(strsplit(V4, ':')),
                     Mismatch = unlist(strsplit(V5, ':')))])
}

# Displays the frequencies of each mismatch type.
plotMismatchTypeFrequencies = function(mismatchTable) {

}

# Displays the frequencies of positions for the chosen mismatch types.
# (Types can be chosen by inclusion or omission)
plotMismatchPositionFrequencies = function(mismatchTable, includedTypes = list(), omittedTypes = list()) {

}
