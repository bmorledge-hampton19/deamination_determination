# This script takes a sam file and locates all the mismatches (potential deamination events) in it.
# The results are written to a bed file along with a metadata file containing the total number of
# reads analyzed, and the number of mismatches found.
import os, gzip
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.DNA_SequenceHandling import reverseCompliment

strandFromIsReverseComplement = {True:'-', False:'+'}

def samMismatchesToBed(samFilePaths: List[str], verbose = False):
    
    for samFilePath in samFilePaths:

        print(f"\nWorking in {os.path.basename(samFilePath)}")

        gzipped = samFilePath.endswith(".gz")

        # Create output file paths (bed file + metadata)
        if gzipped:
            outputBedFilePath = samFilePath.rsplit('.',2)[0] + "_mismatches_by_read.bed"
            openFunction = gzip.open
        else:
            outputBedFilePath = samFilePath.rsplit('.',1)[0] + "_mismatches_by_read.bed"
            openFunction = open
        metadataFilePath = outputBedFilePath.rsplit('.',1)[0] + ".metadata"



        # Read through the sam file line by line, looking for mismatches and recording them.
        alignedReadsCounter = 0
        mismatchesCounter = 0
        with openFunction(samFilePath, "rt") as samFile:
            with open(outputBedFilePath, 'w') as outputBedFile:
                for line in samFile:

                    # Skip header lines.
                    if line.startswith('@'):
                        if verbose: print("Skipping header")
                        continue

                    splitLine = line.split()

                    # Skip lines that didn't align.
                    if splitLine[2] == '*': 
                        if verbose: print("Skipping unaligned read")
                        continue

                    # Increment the aligned reads counter
                    alignedReadsCounter += 1

                    # Using the "MD" optional field, determine if there were any mismatches.
                    # For those those reads that don't contain mismatches (only numeric characters), skip them.
                    mismatchDesignations = splitLine[17].rsplit(':',1)[1]
                    if mismatchDesignations.isnumeric(): 
                        if verbose: print("Skipping read with no mismatches")
                        continue

                    # Get the sequence (and determine if it should actually be the reverse compliment).
                    readSequence = splitLine[9]
                    isReverseCompliment = bool(int(splitLine[1]) & 0b10000)

                    # Next, produce a reference-relative sequence, removing inserted bases and
                    # putting in placeholders for deleted bases.
                    referenceRelativeSequence = ''
                    cigarString = splitLine[5]
                    readPos = 0
                    alphaPositions = [i for i,char in enumerate(cigarString) if char.isalpha()]
                    lastAlphaPosition = -1
                    for alphaPosition in alphaPositions:
                        numeric = int(cigarString[lastAlphaPosition+1:alphaPosition])
                        alpha = cigarString[alphaPosition]
                        if alpha == 'M':
                            referenceRelativeSequence += readSequence[readPos:readPos+numeric]
                            readPos += numeric
                        elif alpha == 'I':
                            readPos += numeric
                        elif alpha == 'D':
                            referenceRelativeSequence += '-'*numeric
                        else: raise ValueError(f"Unexpected cigar string character {alpha} in {cigarString}")
                        lastAlphaPosition = alphaPosition
                    if len(alphaPositions) == 0: referenceRelativeSequence = readSequence
                    refSeqLength = len(referenceRelativeSequence)

                    # Catalogue all the mismatches based on their position in the reference sequence.
                    i = 0
                    refSeqPos = 0
                    mismatches = dict()
                    while i < len(mismatchDesignations):
                        
                        # First, determine how many bases to advance on the reference sequence (matches).
                        numericString = ''
                        while mismatchDesignations[i:i+1].isnumeric():
                            numericString += mismatchDesignations[i]
                            i += 1
                        refSeqPos += int(numericString)

                        # If the next character is a '^', advance the reference position for every directly trailing alpha character.
                        if mismatchDesignations[i:i+1] == '^':
                            i += 1
                            while mismatchDesignations[i:i+1].isalpha():
                                i += 1
                                refSeqPos += 1
                        # Otherwise, directly trailing alpha characters represent mismatches!
                        else:
                            while mismatchDesignations[i:i+1].isalpha():
                                mismatchesCounter += 1
                                mismatches[refSeqPos] = mismatchDesignations[i] + '>' + referenceRelativeSequence[refSeqPos]
                                i += 1
                                refSeqPos += 1

                    # If no mismatches were found (i.e. the MD string only referenced insertions and deletions), skip this line.
                    # NOTE: It turns out you can also do this with the "XM" string, which just counts mismatches.
                    if len(mismatches) == 0: continue

                    # Write the pertinent data for this read to the bed output file.
                    # Adjust values as necessary if the read is on the '-' strand.
                    threePrimeMismatchPositions = list()
                    mismatchSequences = list()
                    for pos in mismatches:
                        if isReverseCompliment:
                            threePrimeMismatchPositions.append(str(-pos - 1))
                            mismatchSequences.append(reverseCompliment(mismatches[pos][0]) + '>' + reverseCompliment(mismatches[pos][2]))
                        else:
                            threePrimeMismatchPositions.append(str(pos - refSeqLength))
                            mismatchSequences.append(mismatches[pos])

                    outputBedFile.write('\t'.join((splitLine[2], str(int(splitLine[3]) - 1), str(int(splitLine[3]) - 1 + refSeqLength), 
                                                   ':'.join(threePrimeMismatchPositions), ':'.join(mismatchSequences),
                                                   strandFromIsReverseComplement[isReverseCompliment])) + '\n')

                    if verbose:
                        print(f"Byte flag: {splitLine[1]}")
                        print(f"Read sequence: {readSequence}")
                        print(f"CIGAR string: {cigarString}")
                        print(f"Derived reference sequence: {referenceRelativeSequence}")
                        print(f"MD string: {splitLine[17]}")
                        print("Bed line:", '\t'.join((splitLine[2], str(int(splitLine[3]) - 1), str(int(splitLine[3]) - 1 + refSeqLength),
                                                      ':'.join(threePrimeMismatchPositions), ':'.join(mismatchSequences),
                                                      strandFromIsReverseComplement[isReverseCompliment])) + '\n')
                        pass


        # Write the metadata
        with open(metadataFilePath, 'w') as metadataFile:
            metadataFile.write(f"Total_Aligned_Reads:\t{alignedReadsCounter}\n")
            metadataFile.write(f"Mismatches:\t{mismatchesCounter}\n")


def main():
    # Create the Tkinter dialog.
    dialog = TkinterDialog(workingDirectory=os.path.join(__file__,"..",".."))
    dialog.createMultipleFileSelector("Sam Read Files:",0,".sam.gz",("Sam Files",(".sam.gz",".sam")), 
                                      additionalFileEndings = [".sam"])

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections = dialog.selections
    samMismatchesToBed(selections.getFilePathGroups()[0])


if __name__ == "__main__": main()