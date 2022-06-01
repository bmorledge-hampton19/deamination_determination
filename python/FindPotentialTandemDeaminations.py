# This script takes fastq files and a table of acceptable read lengths and CPD positions (for XR-seq reads)
# and finds positions that could be tandem CC>TT mismatches caused by deamination.
# For every potential tandem deamination event, it rewrites the fastq read with the event reversed.
# The new read is given a new "sub-identifier" so reads with multiple potential deamination events
# can have the resulting reads from those events grouped together.
import os, gzip
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.DNA_SequenceHandling import reverseCompliment

def findPotentialTandemDeaminations(inputFastqFilePaths: List[str], validLengthsAndPositionsFilePath: str,
                                    outputDir = None):
    
    for inputFastqFilePath in inputFastqFilePaths:

        print(f"\nWorking in {os.path.basename(inputFastqFilePath)}")

        isGzipped = inputFastqFilePath.endswith(".gz")

        # Create the output file path
        if outputDir is None: outputDir = os.path.dirname(inputFastqFilePath)

        if isGzipped:
            outputFileBasename = os.path.basename(inputFastqFilePath).rsplit('.',2)[0] + "_potential_tandem_deaminations.tsv.gz"
            openFunction = gzip.open
        else:
            outputFileBasename = os.path.basename(inputFastqFilePath).rsplit('.',1)[0] + "_potential_tandem_deaminations.tsv"
            openFunction = open

        outputFastqFilePath = os.path.join(outputDir, outputFileBasename)

        # Retrieve valid read lengths and lesion positions from the given file.
        validPositionsByValidReadLength = dict()
        with open(validLengthsAndPositionsFilePath, 'r') as validLengthsAndPositionsFile:
            for line in validLengthsAndPositionsFile:
                pass

        # Read through the fastq file line by line, looking for reads with positions where a CC>TT mismatch could have occurred and
        # rewriting them with the potential mismatch reversed.
        with openFunction(inputFastqFilePath, "rt") as inputFastqFile:
            with openFunction(outputFastqFilePath, 'wt') as outputFastqFile:
                for line in inputFastqFile:
                    pass

                    

def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory=os.path.join(__file__,"..","..")) as dialog:
        dialog.createMultipleFileSelector("Trimmed Fastq Files:",0,"trimmed.fastq.gz",
                                          ("Fastq Files",(".fastq.gz",".fastq")), 
                                        additionalFileEndings = ["trimmed.fastq"])
        dialog.createFileSelector("Valid read lengths and lesion positions:", 0,
                                  ("tsv files", ".tsv"))
        with dialog.createDynamicSelector(2, 0) as outputDirDynSel:
            outputDirDynSel.initCheckboxController("Specify single output dir")
            outputDirDialog = outputDirDynSel.initDisplay(True, "outputDir")
            outputDirDialog.createFileSelector("Output Directory:", 0, directory = True)

    # Get the user's input from the dialog.
    selections = dialog.selections

    if outputDirDynSel.getControllerVar(): outputDir = selections.getIndividualFilePaths("outputDir")[0]
    else: outputDir = None

    findPotentialTandemDeaminations(selections.getFilePathGroups()[0], selections.getIndividualFilePaths()[0], outputDir)


if __name__ == "__main__": main()