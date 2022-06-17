# This script takes a sam file and writes read length and nucleotide frequency for each read to an output file.
# Reads can be filtered on read length, number of mismatches, and the presence of 'N' nucleotides as well.
import os, gzip
from typing import List
from benbiohelpers.DNA_SequenceHandling import reverseCompliment
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.InputParsing.ParseToIterable import parseToIterable


def getReadLengthVsGCContent(samFilePaths: List[str], outputDir = None,
                             acceptableReadLengths = None, maxMismatches = None,
                             allowN = False):

    for samFilePath in samFilePaths:

        print(f"\nWorking in {os.path.basename(samFilePath)}")

        # Prep for file reading and writing
        gzipped = samFilePath.endswith(".gz")

        if outputDir is None: outputDir = os.path.dirname(samFilePath)

        if gzipped:
            outputFileBasename = os.path.basename(samFilePath).rsplit('.',2)[0] + "_read_length_vs_nuc_freq.tsv"
            openFunction = gzip.open
        else:
            outputFileBasename = os.path.basename(samFilePath).rsplit('.',1)[0] + "_read_length_vs_nuc_freq.tsv"
            openFunction = open
        
        outputFilePath = os.path.join(outputDir, outputFileBasename)

        # Iterate through the input file, checking against the filtering options, and outputting read
        # length and nucleotide frequencies if valid
        with openFunction(samFilePath, "rt") as samFile:
            with open(outputFilePath, 'w') as outputFile:

                # Write headers for the output file
                outputFile.write("Read_Length\tA_Counts\tC_Counts\tG_Counts\tT_Counts\n")

                for line in samFile:

                    # Skip header lines
                    if line.startswith('@'): continue

                    splitLine = line.split()

                    # Skip lines that didn't align.
                    if splitLine[2] == '*': continue

                    # Find the XM field and get the number of mismatches from it to see if the read passes filtering
                    if splitLine[13].startswith("XM"):
                        if int(splitLine[13].rsplit(':',1)[1]) > maxMismatches: continue

                    elif splitLine[14].startswith("XM"): 
                        if int(splitLine[14].rsplit(':',1)[1]) > maxMismatches: continue

                    else: raise ValueError(f"XM field not at expected position:\n{line}")

                    # Get the read sequence length and ensure that it passes filtering.
                    readSequence = splitLine[9].upper()
                    if len(readSequence) not in acceptableReadLengths: continue
                    if not allowN and 'N' in readSequence: continue
                    if int(splitLine[1]) & 0b10000: readSequence = reverseCompliment(readSequence)

                    aCount = readSequence.count('A')
                    cCount = readSequence.count('C')
                    gCount = readSequence.count('G')
                    tCount = readSequence.count('T')

                    # Write read length and nucleotide frequencies to the output file.
                    outputFile.write(f"{len(readSequence)}\t{aCount}\t{cCount}\t{gCount}\t{tCount}\n")
                    

def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory=os.path.join(os.path.dirname(os.path.dirname(__file__)),"data")) as dialog:
        dialog.createMultipleFileSelector("Sam Read Files:",0,".sam.gz",("Sam Files",(".sam.gz",".sam")), 
                                          additionalFileEndings = [".sam"])
        dialog.createTextField("Acceptable Read Lengths: ", 1, 0, defaultText = "23-31")
        dialog.createTextField("Max Mismatches: ", 2, 0, defaultText = "2")
        dialog.createCheckbox("Allow N's in reads: ", 3, 0)

        with dialog.createDynamicSelector(4, 0) as outputDirDynSel:
            outputDirDynSel.initCheckboxController("Specify single output dir")
            outputDirDialog = outputDirDynSel.initDisplay(True, "outputDir")
            outputDirDialog.createFileSelector("Output Directory:", 0, directory = True)

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections = dialog.selections
    
    samFilePaths = selections.getFilePathGroups()[0]
    
    acceptableReadLengths = parseToIterable(selections.getTextEntries()[0], castType = int)
    maxMismatches = int(selections.getTextEntries()[1])
    allowN = selections.getToggleStates()[0]
    if outputDirDynSel.getControllerVar(): outputDir = selections.getIndividualFilePaths("outputDir")[0]
    else: outputDir = None

    getReadLengthVsGCContent(samFilePaths, outputDir, acceptableReadLengths, maxMismatches, allowN)


if __name__ == "__main__": main()