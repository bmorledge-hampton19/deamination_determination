# This script takes a sam file with exact matches (aligned from a potential mismatches fastq file) and
# records uniquely mapped read positions to a new bed file.

import os, gzip
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.DNA_SequenceHandling import reverseCompliment
from benbiohelpers.CustomErrors import UserInputError


strandFromIsReverseComplement = {True:'-', False:'+'}

def potentialTandemMismatchesSamToBed(samFilePaths: List[str], outputDir = None):
    
    for samFilePath in samFilePaths:

        print(f"\nWorking in {os.path.basename(samFilePath)}")

        if samFilePath.endswith(".gz"): openFunction = gzip.open
        else: openFunction = open

        # Create output file paths (bed file + metadata)
        if outputDir is None: thisOutputDir = os.path.dirname(samFilePath)
        else: thisOutputDir = outputDir
        mismatchType = os.path.basename(samFilePath).rsplit("_mismatches.sam", 1)[0].rsplit("potential_", 1)[1]
        outputBedFileBasename = (os.path.basename(samFilePath).rsplit(f"potential_{mismatchType}_mismatches.sam",1)[0] +
                                 f"{mismatchType}_mismatches_by_read.bed")
        outputBedFilePath = os.path.join(thisOutputDir, outputBedFileBasename)

        # Read through the sam file line by line, looking for mismatches and recording them.
        currentReadCategory = None
        readsInCurrentCategory = None
        observedReadCategories = set()
        lastLine = None
        writeLastLine = False
        with openFunction(samFilePath, "rt") as samFile:
            with open(outputBedFilePath, 'w') as outputBedFile:
                for line in samFile:

                    # Skip header lines.
                    if line.startswith('@'):
                        continue

                    splitLine = line.split()

                    # Look for the unaligned flag and raise an error if it is found.
                    # (only aligned reads should be reported.)
                    if int(splitLine[1]) & 0b100: 
                        raise UserInputError("Unaligned flag found. Only aligned reads should be present.")

                    # Check the current read category to determine if it has changed, if the last category needs
                    # to be written, etc.
                    thisLineReadCategory = splitLine[0].split('$')[0]

                    if currentReadCategory is None:
                        currentReadCategory = thisLineReadCategory
                        observedReadCategories.add(currentReadCategory)
                        readsInCurrentCategory = 1

                    elif thisLineReadCategory != currentReadCategory:
                        if readsInCurrentCategory == 1 and lastLine.split()[0].split('$')[1]:
                            writeLastLine = True
                        currentReadCategory = thisLineReadCategory
                        if currentReadCategory in observedReadCategories:
                            raise UserInputError("Sam file is not clustered by read category. "
                                                 f"Category {currentReadCategory} found in separate clusters.")
                        else: observedReadCategories.add(currentReadCategory)
                        readsInCurrentCategory = 1

                    else: readsInCurrentCategory += 1

                    # If the last line was flagged for writing, write it!
                    if writeLastLine:
                        lastSplitLine = lastLine.split()

                        readSequence = lastSplitLine[9]
                        if bool(int(lastSplitLine[1]) & 0b10000):
                            readSequence = reverseCompliment(readSequence)
                            strand = '-'
                        else: strand = '+'

                        mismatchstartPos, mismatchEndPos = lastSplitLine[0].split('$')[1].split(':')
                        mismatchCenterPos = (int(mismatchstartPos) + int(mismatchEndPos) - 1) / 2
                        mismatchCenterPos = mismatchCenterPos - len(readSequence) # Convert to 3' oriented (negative) position

                        outputBedFile.write('\t'.join((lastSplitLine[2],
                                                       str(int(lastSplitLine[3]) - 1),
                                                       str(int(lastSplitLine[3]) - 1 + len(readSequence)), 
                                                       str(mismatchCenterPos),
                                                       mismatchType.replace("_to_",'>'),
                                                       strand,
                                                       readSequence)) + '\n')
                        
                        writeLastLine = False

                    # Store this line for later.
                    lastLine = line

def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory=os.path.join(__file__,"..",".."),
                       title = "Potential Tandem Mismatches: Sam to Bed") as dialog:
        dialog.createMultipleFileSelector("Sam Read Files:",0,".sam.gz",("Sam Files",(".sam.gz",".sam")), 
                                        additionalFileEndings = [".sam"])
        with dialog.createDynamicSelector(1, 0) as outputDirDynSel:
                outputDirDynSel.initCheckboxController("Specify single output dir")
                outputDirDialog = outputDirDynSel.initDisplay(True, "outputDir")
                outputDirDialog.createFileSelector("Output Directory:", 0, directory = True)

    # Get the user's input from the dialog.
    selections = dialog.selections

    if outputDirDynSel.getControllerVar(): outputDir = selections.getIndividualFilePaths("outputDir")[0]
    else: outputDir = None

    potentialTandemMismatchesSamToBed(selections.getFilePathGroups()[0], outputDir)


if __name__ == "__main__": main()