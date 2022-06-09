# deamination_determination
Using mismatches in XR-seq reads to pinpoint CPD positions
***

## Table of Contents
1. [Project Overview](#project-overview)
2. [Directory Structure and Naming Conventions](#directory-structure-and-naming-conventions)
3. [Aligning Reads](#aligning-reads)
4. [Identifying Mismatches](#identifying-mismatches)
5. [Mismatch Analysis and Figure Generation](#mismatch-analysis-and-figure-generation)
6. [Data Availability](#data-availability)
7. [Acknowledgements](#acknowledgements)
***

## Project Overview

#### Accompanying Literature
The code in this repository was used to produce the findings in [NO LINK YET](). This paper contains additional insight into the analysis using this code.

#### Background
Cytosine nucleotides deaminate much faster in cyclobutane pyrimidine dimers (CPDs) than normal. Because of this, XR-sequencing reads, which represent fragments from nucleotide excision repair of various lesions (including CPDs) may contain C>T mismatches where the CPD was present. Normally, the position of these lesions in the original repair fragment can only be estimated, but by identifying these C>T mismatches, the lesion positions can be pinpointed more accurately.

#### The Methodology
For single-nucleotide mismatches, reads are aligned using default alignment parameters in [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), which are somewhat permissive of mismatches. The exact locations of the mismatches are extracted post-alignment from the resulting sam file.

For tandem mismatches (specifically, CC>TT) an alternative method must be used. Aligners such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) do not allow multiple mismatches during seeding, so depending on the read length and the position of the tandem mismatch, the read may not align at all. However, given the specific nature of the tandem mismatch, the problem can be coereced into an exact matching situation. For every TT sequence in a given read, a new read is produced with the TT changed to a CC to represent a reversion of the potential mismatch. These reverted reads are check against a reference genome for exact matches, and if exactly one read matches (and the original read does NOT match), the tandem mismatch in that read is confirmed. Exact matching can be performed using bowtie2 with specific parameters outlined in the (Aligning Reads)[#aligning-reads] section.
***

## Directory Structure and Naming Conventions
***

## Aligning Reads
***

## Identifying Mismatches
***

## Mismatch Analysis and Figure Generation
***

## Data Availability
***

## Acknowledgements
