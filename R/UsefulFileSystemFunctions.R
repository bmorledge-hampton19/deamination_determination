library(stringr)

# A list containing the strings designated to files containing a specific type of data.
dataTypeStrings =
  list(mismatchesByRead = "mismatches_by_read",
       simplifiedMismatches = "simplified_mismatches",
       mismatchPositionStats = "mismatch_position_stats",
       nucleosomePeriodicity = "deamination_periodicity_results",
       nucleotideFrequencyData = "nucleotide_frequencies",
       sequenceEnrichment = "sequence_enrichment",
       sequenceFrequency = "sequence_frequency")

generateFilePath = function(directory, dataTypeString, fileExtension,
                            cellType = NA, lesion = NA, timepoint = NA, repitition = NA,
                            includedMismatches = character(), omittedMismatches = character(),
                            strandPolarity = NA, anchored = F, sequence = NA, expansionNum = NA,
                            additionalInformation = character()) {

  if ( length(includedMismatches) > 0 && length(omittedMismatches) > 0) {
    stop("Included mismatches and omitted mismatches given simultaneously")
  }

  includedMismatches = lapply(includedMismatches, function(x) str_replace(x, '>', "_to_"))
  omittedMismatches = lapply(omittedMismatches, function(x) str_replace(x, '>', "_to_"))

  if (!is.na(expansionNum)) expansionText = paste0(expansionNum,"bp_expanded")
  else expansionText = NA

  if (anchored) anchoredText = "anchored"
  else anchoredText = NA

  basenameParts = c(cellType, lesion, timepoint, repitition, includedMismatches, omittedMismatches, "omitted",
                    strandPolarity, anchoredText, sequence, additionalInformation, dataTypeString, expansionText)
  if (length(omittedMismatches) == 0) basenameParts[5+length(includedMismatches)] = NA
  basenameParts = basenameParts[!is.na(basenameParts)]

  return(file.path(directory,paste0(paste(basenameParts,collapse = '_'),fileExtension)))

}
