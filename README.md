# deamination_determination
Using mismatches in XR-seq reads to pinpoint CPD positions
***

## Table of Contents
1. [Project Overview](#project-overview)
2. [Dependencies and File System Naming Conventions](#dependencies-and-file-system-naming-conventions)
3. [Preparing Data for Analysis](#preparing-data-for-analysis)
4. [Analysis and Figure Generation](#analysis-and-figure-generation)
5. [Data Availability](#data-availability)
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

#### A Disclaimer...
Admittedly, the software in this repository was largely not designed with other users in mind. Despite this, the core analysis is still very much reproducible after taking a few steps to set up a specific environment and following a few conventions pertaining to the project's file system structure. In the event that the following steps are insufficient to reproduce the desired results, you are more than welcome to [post an issue](../../issues) or email me directly at b.morledge-hampton@wsu.edu. It is also worth noting that we plan to produce a more user-friendly software package with the same functionality in the near future.

#### Dependencies
Besides the code in this repository, the analysis requires that the following packages are installed (The most up-to-date version of each package is recommended):
- Python:
  - [benbiohelpers](https://github.com/bmorledge-hampton19/benbiohelpers)
- R:
  - [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
  - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- Must be available through command line (Recommend installing through apt, where possible):
  - [bedtools](https://bedtools.readthedocs.io/en/latest/)
  - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - [samtools](http://www.htslib.org/)
  - [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  
Although not required, the scripts in [this repository](https://github.com/bmorledge-hampton19/XR-seq_Analysis) may be helpful for [trimming](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/TrimAdaptorSequences.py) and [aligning](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/TrimmedFastqToSam.py) sequencing reads and [parsing them to bed](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/SamToBed.py). If desired, these three operations can be chained together using the [AlignXRSeqReads.py](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/AlignXRSeqReads.py) script.

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

The corresponding data sets for the above names can be found in the [Data Availability](#data-availability) section

Within each of these directories, up to two other directories should be created:
- "mismatches_by_read" for storing mismatch data by read (e.g. The results from [SamMismatchesToBed.py](python/SamMismatchesToBed.py) or [PotentialMismatchesSamToBed.py](python/PotentialMismatchesSamToBed.py))
- "general_read_sequences" for storing bed files containing general non-mismatch-associated sequencing results

When running the R notebooks contained in this repository, they will ask you to provide a "bioinformatics directory". This is merely the directory containing the cloned repository (i.e. one directory up from "deamination_determination/"). Once this directory is selected, the underlying file system is inferred using the naming conventions in this section.

#### Data file naming conventions
When obtaining fastq sequencing files (or sam alignment files), a strict file naming system needs to be maintained in order for the rest of the pipeline to run without errors. The file names should be made up of the following identifying information, separated by underscores:
- Organism or cell type name
- Lesion type ("CPD", "6-4", or "cisplatin")
- Timepoint (in shorthand format, such as "10min" or "1h")
- Repitition number (In the format "rep#" where '#' is the repitition number or "all_reps" when combined)

Here is an example of a valid sam file name: "NHF1_CPD_1h_rep1.sam"

#### Deprecated Periodicity Analysis
The [NHF1 periodicity analysis](R/NHF1DeaminationPeriodicityAnalysis.rmd) is largely deprecated, but if for some reason you want to run it anyway, you will need to install [mutperiod](https://github.com/bmorledge-hampton19/mutperiod) for the periodicity analysis and have [this repository](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis) downloaded in the same directory as deamination determination. The analysis requires that you use mutperiod to generate the periodicity data and stratify it by sequence context. Seriously, probably don't try to run this. The analysis turned up very little in terms of a clear periodicity.
***

## Preparing Data for Analysis

#### Aligning Reads
For this analysis, reads were aligned using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) after trimming adaptor sequences with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Other aligners have not been tested, but if they alter the sam file format significantly, they may not function properly with the rest of the pipeline.

For finding single-nucleotide mismatches, default bowtie2.4.5 parameters were used. For finding exact matches with potential tandem mismatches, the following parameters were used: `--reorder --no-1mm-upfront --no-unal --score-min C,0,0 -k 2`.

#### Identifying Single Nucleotide Mismatches
In order to identify single nucleotide mismatches from an aligned sam file, use the [SamMismatchesToBed.py](python/SamMismatchesToBed.py) script. When running the script directly (recommended), you should see a short dialog asking for a few input parameters. Sam files can be given by selecting individual files or by giving a directory which will be recursively searched for files ending in ".sam" or ".sam.gz". The "omit indels" box should be checked, and the "specify single output dir" box can be used to output the results directly to the relevant "mismatches_by_read" directory, if desired.

#### Identifying Tandem Mismatches
In order to identify tandem nucleotide mismatches, reads were first trimmed, and then the [FindPotentialTandemMismatches.py](python/FindPotentialTandemMismatches.py) script was used to split reads based on potential tandem mismatch sites. Fastq files can be given by selecting individual files or by giving a directory which will be recursively searched for files ending in "trimmed.fastq" or "trimmed.fastq.gz". Next, the resulting reads should be aligned with bowtie using the parameters outlined in the [Aligning Reads](#Aligning-Reads) section. Finally, true tandem mismatches should be parsed from the resulting sam file using the [PotentialTandemMismatchesSamToBed.py](python/PotentialTandemMismatchesSamToBed.py) script.

#### Producing General Read Data
Some steps of the analysis require data in bed format from all reads, independent of mismatches. Aligned reads reads can be coerced to bed format using a combination of samtools and bedtools. A script for this conversion is available [here](https://github.com/bmorledge-hampton19/XR-seq_Analysis/blob/main/SamToBed.py).

#### Expanding and Filtering Data
Once the data from the above steps has been generated and is present in the [relevant directories](#directory-structure-naming-conventions), the sequence context of the data should be expanded by at least 3 base pairs using [this script](https://github.com/bmorledge-hampton19/benbiohelpers/blob/main/python/benbiohelpers/FileSystemHandling/ExpandSequenceContext.py). For mismatch data files, the expanded sequence context should replace the 7th column (index 6), and for general read files, the new sequence should be appended. For general reads, the "output file suffix" should be in the format "\_aligned_reads_#bp_expanded", whereas for mismatch data, the suffix should be in the format "\_#bp_expanded". After being expanded, data can be optionally TGG filtered using the [FilterOnExpandedContext.py](python/FilterOnExpandedContext.py) script.

#### Combining Data Repititions
Lastly, the individual data repititions need to be combined. This can be achieved by running [this script](https://github.com/bmorledge-hampton19/benbiohelpers/blob/main/python/benbiohelpers/FileSystemHandling/CombineReps.py) on the "mismatches_by_read" and "general_read_sequences" directories. This step IS REQUIRED to run the analyses within the R notebooks
***

## Analysis and Figure Generation
If the data has been prepared correctly, running the analyses and generating figures is easy: Simply run/knit the R notebooks provided in the [R/](R/) directory of the repository for the desired organism and analysis. However, there are a few caveats here...

#### Order matters
Any "DataGeneration" notebooks must be run BEFORE their "FigureGeneration" counterparts. If a "FigureGeneration" notebook is run with specific parameters (i.e. TGG filtering), then its corresponding "DataGeneration" notebook must have also been run with those parameters at least once.  
The "DeaminationSequenceAndPositionAnalysis" notebooks must be run BEFORE their "DeaminationAnchoredExcisionSite" counterparts.

#### Defining the "Bioinformatics_Projects" Directory
Each notebook requires a path to the "Bioinformatics_Projects" directory. (More information on this directory is available at the end of the [directory structure naming conventions section](#directory-structure-naming-conventions). In brief, this is the parent directory to the cloned deamination_determination repository.  
If you are running the R notebooks in interactive mode, you will be able to select this directory through a file browser dialog. If you are knitting the notebook, you will select this directory from a list of choices in the notebook's Shiny UI. If the desired directory is not listed, you will need to add it manually to the R notebook (generally on line 11).

#### Figure Output
When knitting figure generation notebooks, the figures are automatically output to the R/deamination_analysis_output directory, which will be created if it doesn't already exist.
***

## Data Availability
The sequencing data used for the analyses supported by this repository can be found at the following locations:
- Arabidopsis: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108932)
- CS-B: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67941](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67941)
- Dm: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138846)
- E_coli_WT: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92734](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92734)
- GM12878: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82213)
- NHF1: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76391)
- XP-C: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67941](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67941)
- yeast: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110621)
