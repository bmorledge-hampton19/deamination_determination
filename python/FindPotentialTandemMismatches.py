# This script takes fastq files (for XR-seq reads)
# and finds positions that could be tandem mismatch based on given mismatched and reference sequences.
# For every potential tandem mismatch event, the fastq read is rewritten with the event reversed.
# The new read is given a new "sub-identifier" so reads with multiple potential deamination events
# can have the resulting reads from those events grouped together.
import os, gzip
from typing import IO, List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CustomErrors import UserInputError


def writeNewFastqEntry(file: IO, readID, remainingHeaderLine, readSequence, qualityScoreLine):
    """
    This function writes an updated fastq read entry using the given information.
    NOTE: This function will NOT make changes to read ID or read sequence based on the presence of potential mismatches.
    """
    file.write('@' + readID + ' ' + remainingHeaderLine)
    file.write(readSequence + '\n')
    file.write('+' + readID + ' ' + remainingHeaderLine)
    file.write(qualityScoreLine)


def findPotentialTandemDeaminations(inputFastqFilePaths: List[str], referenceSequence, mismatchSequence, outputDir = None):

    for inputFastqFilePath in inputFastqFilePaths:

        print(f"\nWorking in {os.path.basename(inputFastqFilePath)}")
        if "trimmed.fastq" not in inputFastqFilePath:
            raise UserInputError("Given fastq file does not appear to be trimmed. (Expected \"trimmed.fastq\" in the name.)")

        isGzipped = inputFastqFilePath.endswith(".gz")

        # Create the output file path
        if outputDir is None: outputDir = os.path.dirname(inputFastqFilePath)

        if isGzipped:
            outputFileBasename = os.path.basename(inputFastqFilePath).rsplit("trimmed.fastq",1)[0]
            outputFileBasename += "_potential_tandem_deaminations.fastq.gz"
            openFunction = gzip.open
        else:
            outputFileBasename = os.path.basename(inputFastqFilePath).rsplit("trimmed.fastq",1)[0]
            outputFileBasename += "_potential_tandem_deaminations.fastq"
            openFunction = open

        outputFastqFilePath = os.path.join(outputDir, outputFileBasename)

        # Read through the fastq file line by line, looking for reads with positions where a CC>TT mismatch could have occurred
        # and rewriting them with the potential mismatch reversed.
        with openFunction(inputFastqFilePath, "rt") as inputFastqFile:
            with openFunction(outputFastqFilePath, 'wt') as outputFastqFile:
                for line in inputFastqFile:
                    
                    # Get information from the four lines that represent the current fastq read.
                    readID, remainingHeaderLine = line.split(' ', 1)
                    readID = readID[1:]
                    readSequence = inputFastqFile.readline().strip()
                    inputFastqFile.readline()
                    qualityScoreLine = inputFastqFile.readline()

                    # Write the unaltered fastq entry to the new file (with '$' delimiter added to the read ID).
                    writeNewFastqEntry(outputFastqFile, readID + '$', remainingHeaderLine, readSequence, qualityScoreLine)

                    # Search for the given mismatch sequence in the read, and wherever it is found, create a new read with the 
                    # mismatch sequence replaced by the reference sequence.
                    potentialMismatchLocation = readSequence.find(mismatchSequence)
                    while potentialMismatchLocation != -1:
                        thisReadID = f"{readID}${potentialMismatchLocation}:{potentialMismatchLocation+len(mismatchSequence)}"
                        thisReadSequence = (readSequence[:potentialMismatchLocation] +
                                            referenceSequence +
                                            readSequence[potentialMismatchLocation+2:])
                        writeNewFastqEntry(outputFastqFile, thisReadID, remainingHeaderLine,
                                           thisReadSequence, qualityScoreLine)
                        potentialMismatchLocation = readSequence.find(mismatchSequence, potentialMismatchLocation + 1)


def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory=os.path.join(__file__,"..","..")) as dialog:
        dialog.createMultipleFileSelector("Trimmed Fastq Files:",0,"trimmed.fastq.gz",
                                          ("Fastq Files",(".fastq.gz",".fastq")), 
                                        additionalFileEndings = ["trimmed.fastq"])
        dialog.createTextField("Reference Sequence: ", 1, 0, defaultText = "CC")
        dialog.createTextField("Mismatched Sequence: ", 2, 0, defaultText = "TT")
        with dialog.createDynamicSelector(3, 0) as outputDirDynSel:
            outputDirDynSel.initCheckboxController("Specify single output dir")
            outputDirDialog = outputDirDynSel.initDisplay(True, "outputDir")
            outputDirDialog.createFileSelector("Output Directory:", 0, directory = True)

    # Get the user's input from the dialog.
    selections = dialog.selections

    if outputDirDynSel.getControllerVar(): outputDir = selections.getIndividualFilePaths("outputDir")[0]
    else: outputDir = None

    findPotentialTandemDeaminations(selections.getFilePathGroups()[0], *selections.getTextEntries(), outputDir)


if __name__ == "__main__": main()