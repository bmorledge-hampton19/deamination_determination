library(stringr)
library(data.table)

# A list containing the strings designated to files containing a specific type of data.
dataTypeStrings =
  list(mismatchesByRead = "mismatches_by_read",
       mismatchesByReadRelationFormatted = "mismatches_by_read_relation_formatted",
       simplifiedMismatches = "simplified_mismatches",
       mismatchPositionStats = "mismatch_position_stats",
       nucleosomePeriodicity = "deamination_periodicity_results",
       nucleotideFrequencyData = "nucleotide_frequencies",
       sequenceEnrichment = "sequence_enrichment",
       sequenceFrequency = "sequence_frequency",
       alignedReads = "aligned_reads",
       mismatchFrequencyZScores = "mismatch_frequency_z-scores",
       readLengthcounts = "read_length_counts",
       tandemMismatchCountsComparison = "tandem_mismatch_counts_comparison",
       sequenceLogoInput = "sequence_logo_input",
       featureRelativeData = "feature_relative_data",
       featureRelativePositionCounts = "feature_relative_position_counts",
       featureRelativeCutSiteCounts = "feature_relative_cut_site_counts",
       featureRelativeMeanCutSiteDistance = "feature_relative_mean_cut_site_distance",
       pairedCutSiteDistances = "paired_cut_site_distances"
       )

generateFilePath = function(directory, dataTypeString, fileExtension,
                            cellType = NA, lesion = NA, timepoint = NA, repitition = NA,
                            includedMismatches = character(), omittedMismatches = character(),
                            strandPolarity = NA, anchored = F, sequence = NA, expansionNum = NA,
                            additionalInformation = character(), filtering = NA) {

  if (!is.na(timepoint) && timepoint == "NONE") timepoint = NA

  if ( length(includedMismatches) > 0 && length(omittedMismatches) > 0) {
    stop("Included mismatches and omitted mismatches given simultaneously")
  }

  includedMismatches = lapply(includedMismatches, function(x) str_replace_all(x, '>', "_to_"))
  omittedMismatches = lapply(omittedMismatches, function(x) str_replace_all(x, '>', "_to_"))

  if (!is.na(expansionNum) && expansionNum > 0) expansionText = paste0(expansionNum,"bp_expanded")
  else expansionText = NA

  if (anchored) anchoredText = "anchored"
  else anchoredText = NA

  basenameParts = c(cellType, lesion, timepoint, repitition, includedMismatches, omittedMismatches, "omitted",
                    strandPolarity, anchoredText, sequence, dataTypeString, expansionText, filtering,
                    additionalInformation)
  if (length(omittedMismatches) == 0) basenameParts[5+length(includedMismatches)] = NA
  basenameParts = basenameParts[!is.na(basenameParts)]

  return(file.path(directory,paste0(paste(basenameParts,collapse = '_'),fileExtension)))

}


# Reads in data using fread and the generateFilePath function, but with a twist:
# Parameters can be passed as lists, and when they are, the returned data object
# will be a list (potentially a list of lists) with names/tiers corresponding to the different parameters.
# NOTE: Currently only supports stratification by timepoint and repitition.
# NOTE NOTE: This really would be easier with a class structure for data objects.
readData = function(directory, dataTypeString, fileExtension,
                    cellType = NA, lesion = NA, timepoint = NA, repitition = NA,
                    includedMismatches = character(), omittedMismatches = character(),
                    strandPolarity = NA, anchored = F, sequence = NA, expansionNum = NA,
                    additionalInformation = character(), filterText = NA) {

}


# Writes data using the given parameters, but with a fun and exciting twist:
# If the input data is in list format with names corresponding to parameters,
# ?????? How will this work ??????

