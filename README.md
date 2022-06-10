# deamination_determination
Using mismatches in XR-seq reads to pinpoint CPD positions
***

## Table of Contents
1. [Project Overview](#project-overview)
2. [Dependencies and File System Naming Conventions](#dependencies-and-file-system-naming-conventions)
3. [Aligning Reads](#aligning-reads)
4. [Identifying Mismatches](#identifying-mismatches)
5. [Expanding and Filtering Mismatch Data](#expanding-and-filtering-mismatch-data)
6. [Mismatch Analysis and Figure Generation](#mismatch-analysis-and-figure-generation)
7. [Data Availability](#data-availability)
8. [Acknowledgements](#acknowledgements)
***

## Project Overview

#### Accompanying Literature
The code in this repository was used to produce the findings in [NO LINK YET](). This paper contains additional insight into the analysis using this code.

#### Background
Cytosine nucleotides deaminate much faster in cyclobutane pyrimidine dimers (CPDs) than normal. Because of this, XR-sequencing reads, which represent fragments from nucleotide excision repair of various helix-distorting lesions such as CPDs (more information [here](http://genesdev.cshlp.org/content/29/9/948)) may contain C>T mismatches where the CPD was present. Normally, the position of these lesions in the original repair fragment can only be estimated, but by identifying these C>T mismatches, the lesion positions can be pinpointed more accurately.

#### The Methodology
For single-nucleotide mismatches, reads are aligned using default alignment parameters in [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), which are somewhat permissive of mismatches. The exact locations of the mismatches are extracted post-alignment from the resulting sam file.

For tandem mismatches (specifically, CC>TT) an alternative method must be used. Aligners such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) do not allow multiple mismatches during seeding, so depending on the read length and the position of the tandem mismatch, the read may not align at all. However, given the specific nature of the tandem mismatch, the problem can be coereced into an exact matching situation. For every TT sequence in a given read, a new read is produced with the TT changed to a CC to represent a reversion of the potential mismatch. These reverted reads are check against a reference genome for exact matches, and if exactly one read matches (and the original read does NOT match), the tandem mismatch in that read is confirmed. Exact matching can be performed using bowtie2 with specific parameters outlined in the [Aligning Reads](#aligning-reads) section.
***

## Dependencies and File System Naming Conventions
Admittedly, the software in this repository was largely not designed with other users in mind. Despite this, the core analysis is still very much reproducible after taking a few steps to set up a specific environment and following a few conventions pertaining to the project's file system structure. In the event that the following steps are insufficient to reproduce the desired results, you are more than welcome to [post an issue](../../issues) or email me directly at b.morledge-hampton@wsu.edu.

#### Dependencies
Besides the code in this repository, the analysis requires that the following packages are installed (The most up-to-date version of each package is recommended):
Python:
  [benbiohelpers](https://github.com/bmorledge-hampton19/benbiohelpers)
R:
  [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
  [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
Must be available through command line (Recommend installing through apt, where possible):
  [bedtools](https://bedtools.readthedocs.io/en/latest/) 
  [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  [samtools](http://www.htslib.org/)
  [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  
Although not required, the scripts in [this repository](https://github.com/bmorledge-hampton19/XR-seq_Analysis) may be helpful for [trimming](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/TrimAdaptorSequences.py) and [aligning](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/TrimmedFastqToSam.py) sequencing reads and [parsing them to bed](NEED_TO_MAKE_THIS_SCRIPT_ALSO_DID_YOU_REMOVE_MUTPERIOD_DEPENDENCIES).
  
The [NHF1 periodicity analysis](R/NHF1DeaminationPeriodicityAnalysis.rmd) is largely deprecated, but if for some reason you want to run it anyway, you will need to install [mutperiod](https://github.com/bmorledge-hampton19/mutperiod) for the periodicity analysis and have [this repository](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis) downloaded in the same directory as deamination determination.

#### Directory structure naming conventions
After cloning the deamination_determination repository, you will need to create a directory at the top level called "data". Within this data directory, an additional directory will be required for each species or cell type being analyzed. The names supported by the R notebooks in this repository are:
- Arabidopsis
- CS-B
- Dm
- E_coli_WT
- GM12878
- NHF1
- XP-C
- yeast

Within each of these directories, up to three other directories should be created:
- "mismatches_by_read" for storing single-base mismatch data by read (e.g. The result from [SamMismatchesToBed.py](python/SamMismatchesToBed.py))
- "CC_to_TT_tandem_deaminations" for storing tandem mismatch data by read (e.g. The result from [PotentialMismatchesSamToBed.py](python/PotentialMismatchesSamToBed.py))
- "general_read_sequences" for storing bed files containing general non-mismatch-associated sequencing results

When running the R notebooks contained in this repository, they will ask you to provide a "bioinformatics directory". This is merely the directory containing the cloned repository (i.e. one directory up from "deamination_deterimnation/"). Once this directory is selected, the underlying file system is inferred using the naming conventions in this section.

#### Data file naming conventions

***

## Aligning Reads
***

## Identifying Mismatches
***

## Expanding and Filtering Mismatch Data
***

## Mismatch Analysis and Figure Generation
***

## Data Availability
***

## Acknowledgements
