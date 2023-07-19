library(data.table)
library(stringr)
library(ggplot2)

THREE_PRIME = "three_prime"
FIVE_PRIME = "five_prime"
HALF_MAX = "half_max"

TFBS_RELATIVE_POSITION = "TFBS_Relative_Mismatch_Pos"
TFBS_RELATIVE_THREE_PRIME_CUT_SITE = "TFBS_Relative_Three_Prime_Cut_Site"
TFBS_RELATIVE_FIVE_PRIME_CUT_SITE = "TFBS_Relative_Five_Prime_Cut_Site"


# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 22, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18, color = "black"),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Takes a file of mismatches related to TFBSs and extracts information on relative mismatch and cut site location.
parseMismatchesByTFBS = function(table) {
  table = copy(table)
  colnames(table) = c("Chromosome", "Mismatch_Pos_0", "Mismatch_Pos_1", "Read_Relative_Mismatch_Pos", "Mismatch_Type",
                      "Mismatch_Strand", "Read_Length", "TFBS_Pos", "TFBS_Strand")

  # Record mismatch positions relative to the transcription factor binding site
  table[TFBS_Strand == '+', TFBS_Relative_Mismatch_Pos := Mismatch_Pos_0 - TFBS_Pos]
  table[TFBS_Strand == '-', TFBS_Relative_Mismatch_Pos := TFBS_Pos - Mismatch_Pos_0]

  # Record cut sites relative to the transcription factor binding site
  table[Mismatch_Strand == '+', Three_Prime_Cut_Site := Mismatch_Pos_0 - Read_Relative_Mismatch_Pos - 0.5]
  table[Mismatch_Strand == '-', Three_Prime_Cut_Site := Mismatch_Pos_0 + Read_Relative_Mismatch_Pos + 0.5]

  table[Mismatch_Strand == '+', Five_Prime_Cut_Site := Mismatch_Pos_0 - (Read_Length + Read_Relative_Mismatch_Pos) - 0.5]
  table[Mismatch_Strand == '-', Five_Prime_Cut_Site := Mismatch_Pos_0 + Read_Length + Read_Relative_Mismatch_Pos + 0.5]

  table[TFBS_Strand == '+', TFBS_Relative_Three_Prime_Cut_Site := Three_Prime_Cut_Site - TFBS_Pos]
  table[TFBS_Strand == '-', TFBS_Relative_Three_Prime_Cut_Site := TFBS_Pos - Three_Prime_Cut_Site]

  table[TFBS_Strand == '+', TFBS_Relative_Five_Prime_Cut_Site := Five_Prime_Cut_Site - TFBS_Pos]
  table[TFBS_Strand == '-', TFBS_Relative_Five_Prime_Cut_Site := TFBS_Pos - Five_Prime_Cut_Site]

  table[,Same_Strand := Mismatch_Strand == TFBS_Strand]

}


# This functions take a parsed table of TFBS-relative mismatch data and return a table of
# counts for the relevant TFBS-relative data, as specified by the "feature" parameter.
getTfbsRelativeCounts = function(table, feature) {

  if (feature == TFBS_RELATIVE_POSITION) {
    return(table[,.N, by = list(TFBS_Relative_Mismatch_Pos, Same_Strand)])
  } else if (feature == TFBS_RELATIVE_THREE_PRIME_CUT_SITE){
    return(table[,.N, by = list(TFBS_Relative_Three_Prime_Cut_Site, Same_Strand)])
  } else if (feature == TFBS_RELATIVE_FIVE_PRIME_CUT_SITE) {
    return(table[,.N, by = list(TFBS_Relative_Five_Prime_Cut_Site, Same_Strand)])
  } else stop("Unrecognized value for feature parameter")

}


# This function plots counts of TFBS-relative features, as specified by the "feature" parameter.
plotTFBSRelativeCounts = function(countsTable, feature, title = NULL, sameStrand = NULL,
                                  ylim = NULL, xlim = NULL, xAxisBreaks = waiver()) {

  countsTable = copy(countsTable)
  if (!is.null(sameStrand)) {
    countsTable = countsTable[Same_Strand == sameStrand]
  }
  countsTable[, Relative_Frequency := N/sum(N)]

  if (feature == TFBS_RELATIVE_POSITION) {
    xAxisLabel = "TFBS Relative Position"
    if (is.null(title)) title = "TFBS Relative Mismatch Position Frequencies"
  } else if (feature == TFBS_RELATIVE_THREE_PRIME_CUT_SITE){
    xAxisLabel = "TFBS Relative 3' Cut Site"
    if (is.null(title)) title = "TFBS Relative 3' Cut Site Position Frequencies"
  } else if (feature == TFBS_RELATIVE_FIVE_PRIME_CUT_SITE) {
    xAxisLabel = "TFBS Relative 5' Cut Site"
    if (is.null(title)) title = "TFBS Relative 5' Cut Site Position Frequencies"
  } else stop("Unrecognized value for feature parameter")

  print(
    ggplot(countsTable, aes(x = !!sym(feature), y = Relative_Frequency)) +
      geom_bar(stat = "identity") + coord_cartesian(ylim = ylim, xlim = xlim) +
      scale_x_continuous(breaks = xAxisBreaks) +
      labs(title = title, x = xAxisLabel, y = "Relative Frequency") +
      blankBackground + defaultTextScaling
  )

}
