# This script takes a bed file resulting from SamMismatchesToBed.py and converts
# it to a more general bed format suitable for mutperiod.
# Also, data can be filtered based on read length, number of mismatches, the sequence of the mismatch,
# And the position of the mismatch within the read.
# Alternatively, this script can look for a specific type of tandem mutation at specific positions.
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs
from benbiohelpers.InputParsing.ParseToIterable import parseToIterable


# Given a bed file path containing information on sam read mismatches, filter the reads and mismatches
# based on the given parameters (all "acceptable..." parameters accept a list/range as input, and if
# any filtering parameters are passed NoneType, they will not be used in filtering at all.) and single
# mismatches and/or tandem mismatches to a custom-bed formatted file.
def samBedToCustomBed(samBedFilePaths: List[str], customSuffix: str = "mismatches_custom_input",
                      acceptableReadLengths = None, maxMismatches = None,
                      acceptableMismatchPositions = None, acceptableMismatchSequences = None):

    for samBedFilePath in samBedFilePaths:

        print(f"\nWorking in {os.path.basename(samBedFilePath)}")

        # Generate the required output file paths, each within its own directory (for convenience with mutperiod).
        inputFileDir = os.path.dirname(samBedFilePath)
        baseFileName = os.path.basename(samBedFilePath).rsplit('.')[0].rsplit("mismatches_by_read", 1)[0]
        outputDir = os.path.join(inputFileDir, baseFileName + customSuffix.rsplit("custom_input", 1)[0][:-1])
        customBedOutputFilePath = os.path.join(outputDir, baseFileName + customSuffix + ".bed")
        checkDirs(outputDir)

        # Loop through the sam bed file, filtering out lines that don't meet the requirements
        # and writing any remaining relevant data to the output file(s)
        with open(customBedOutputFilePath, 'w') as customBedOutputFile:
            with open(samBedFilePath, 'r') as samBedFile:
                for line in samBedFile:

                    (readChrom, readStart, readEnd, mismatchPositions, mismatchTypes, strand, _) = line.split()
                    readStart = int(readStart)
                    readEnd = int(readEnd)
                    
                    # Perform initial filtering on the read.
                    if ( (acceptableReadLengths is None or int(readEnd) - int(readStart) in acceptableReadLengths) and
                        (maxMismatches is None or mismatchPositions.count(':') < maxMismatches) ): pass
                    else: continue

                    # If the read passed filtering, write any mismatches that pass filtering.
                    for relativeMMPosition, mmType in zip(mismatchPositions.split(':'), mismatchTypes.split(':')):

                        relativeMMPosition = float(relativeMMPosition)
                        if round(relativeMMPosition) == relativeMMPosition: relativeMMPosition = int(relativeMMPosition)

                        if ( (acceptableMismatchPositions is None or relativeMMPosition in acceptableMismatchPositions) and
                             (acceptableMismatchSequences is None or mmType in acceptableMismatchSequences) ):

                            posOffset = ( (len(mmType)-1)/2 - 1 ) / 2
                            if strand == '+':
                                absoluteMMPosStart = int(readEnd + relativeMMPosition - posOffset)
                                absoluteMMPosEnd = int(readEnd + relativeMMPosition + posOffset)
                            else:
                                absoluteMMPosStart = int(readStart - relativeMMPosition - 1 - posOffset)
                                absoluteMMPosEnd = int(readStart - relativeMMPosition - 1 + posOffset)

                            customBedOutputFile.write('\t'.join((readChrom, str(absoluteMMPosStart), str(absoluteMMPosEnd+1),
                                                                 *mmType.split('>'), strand, str(relativeMMPosition))) + '\n')


def main():
    # Create the Tkinter dialog.
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(os.path.dirname(__file__)),"data"))
    dialog.createMultipleFileSelector("Mismatched Reads Bed Files:",0,"mismatches_by_read.bed",("Bed Files",".bed"))
    dialog.createTextField("Acceptable Read Lengths: ", 1, 0, defaultText = "23-31")
    dialog.createTextField("Max Mismatches: ", 2, 0, defaultText = "2")
    dialog.createTextField("Acceptable Mismatch Positions: ", 3, 0, defaultText = "-12:-3")
    dialog.createTextField("Acceptable Mismatch Sequences: ", 4, 0, defaultText = "C>T")
    dialog.createTextField("File suffix (before .bed, with \"mismatches_by_read\" omitted):",
                           5, 0, defaultText = "mismatches_custom_input")

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections = dialog.selections
    
    samBedFilePaths = selections.getFilePathGroups()[0]
    
    acceptableReadLengths = parseToIterable(selections.getTextEntries()[0], castType = int)
    if len(acceptableReadLengths) == 0: acceptableReadLengths = None

    maxMismatches = selections.getTextEntries()[1]
    if maxMismatches == '': maxMismatches = None
    else: maxMismatches = int(maxMismatches)

    acceptableMismatchPositions = parseToIterable(selections.getTextEntries()[2],rangeChar = ':', castType = float)
    if len(acceptableMismatchPositions) == 0: acceptableMismatchPositions = None

    acceptableMismatchSequences = parseToIterable(selections.getTextEntries()[3])
    if len(acceptableMismatchSequences) == 0: acceptableMismatchSequences = None

    samBedToCustomBed(samBedFilePaths, selections.getTextEntries()[4],
                      acceptableReadLengths, maxMismatches,
                      acceptableMismatchPositions, acceptableMismatchSequences)


if __name__ == "__main__": main()