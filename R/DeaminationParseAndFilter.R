library(data.table)
library(stringr)

THREE_PRIME = 3
FIVE_PRIME = 5


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