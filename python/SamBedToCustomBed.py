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
def samBedToCustomBed(samBedFilePaths: List[str], 
                      outputSingleBaseMismatches = True, singleBaseSuffix = "_individual_mismatches",
                      outputTandemMismatches = True, tandemSuffix = "_tandem_mismatches",
                      acceptableReadLengths = None, maxMismatches = None, acceptableMismatchPositions = None, 
                      acceptableMismatchSequences = None):

    for samBedFilePath in samBedFilePaths:

        print(f"\nWorking in {os.path.basename(samBedFilePath)}")

        # Generate the required output file paths, each within its own directory (for convenience with mutperiod).
        inputFileDir = os.path.dirname(samBedFilePath)
        baseFileName = os.path.basename(samBedFilePath).rsplit('.')[0].rsplit("_mismatches_by_read")[0]

        if outputSingleBaseMismatches:
            singleMismatchOutputDir = os.path.join(inputFileDir, baseFileName + singleBaseSuffix)
            singleMismatchOutputFilePath = os.path.join(singleMismatchOutputDir, baseFileName + singleBaseSuffix + ".bed")
            checkDirs(singleMismatchOutputDir)
            singleMismatchOutputFile = open(singleMismatchOutputFilePath, 'w')
        if outputTandemMismatches:
            tandemMismatchOutputDir = os.path.join(inputFileDir, baseFileName + tandemSuffix)
            tandemMismatchOutputFilePath = os.path.join(tandemMismatchOutputDir, baseFileName + tandemSuffix + ".bed")
            checkDirs(tandemMismatchOutputDir)
            tandemMismatchOutputFile = open(tandemMismatchOutputFilePath, 'w')

        # Loop through the sam bed file, filtering out lines that don't meet the requirements
        # and writing any remaining relevant data to the output file(s)
        with open(samBedFilePath, 'r') as samBedFile:
            for line in samBedFile:

                (readChrom, readStart, readEnd, mismatchPositions, mismatchSequences, strand) = line.split()
                readStart = int(readStart)
                readEnd = int(readEnd)
                
                # Perform initial filtering on the read.
                if ( (acceptableReadLengths is None or int(readEnd) - int(readStart) in acceptableReadLengths) and
                     (maxMismatches is None or mismatchPositions.count(':') < maxMismatches) ):
                    pass
                else: continue

                # If the read passed filtering, write any single-base and tandem mismatches that pass filtering.
                lastPosition = None
                lastSequence = None
                for position, sequence in zip(mismatchPositions.split(':'), mismatchSequences.split(':')):
                    position = int(position)
                    if ( (acceptableMismatchPositions is None or position in acceptableMismatchPositions) and
                            (acceptableMismatchSequences is None or sequence in acceptableMismatchSequences) ):

                        if outputSingleBaseMismatches:
                            if strand == '+':
                                absolutePos = readEnd + position
                            else:
                                absolutePos = readStart - position - 1

                            singleMismatchOutputFile.write('\t'.join((readChrom, str(absolutePos), str(absolutePos+1),
                                                                      sequence[0], sequence[2], strand, str(position))) + '\n')

                        if outputTandemMismatches and lastPosition is not None and abs(position-lastPosition) == 1:
                            if strand == '+':
                                absolutePos = readEnd + lastPosition
                                refSeq = lastSequence[0] + sequence[0]
                                mutSeq = lastSequence[2] + sequence[2]
                            else:
                                absolutePos = readStart - lastPosition - 1
                                refSeq = sequence[0] + lastSequence[0]
                                mutSeq = sequence[2] + lastSequence[2]

                            tandemMismatchOutputFile.write('\t'.join((readChrom, str(absolutePos), str(absolutePos+2),
                                                                      refSeq, mutSeq, strand, 
                                                                      str((position+lastPosition)/2) )) + '\n')

                        lastPosition = position
                        lastSequence = sequence


def main():
    # Create the Tkinter dialog.
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(os.path.dirname(__file__)),"data"))
    dialog.createMultipleFileSelector("Mismatched Reads Bed Files:",0,"mismatches_by_read.bed",("Bed Files",".bed"))
    dialog.createTextField("Acceptable Read Lengths: ", 1, 0, defaultText = "23-31")
    dialog.createTextField("Max Mismatches: ", 2, 0, defaultText = "2")
    dialog.createTextField("Acceptable Mismatch Positions: ", 3, 0, defaultText = "-12$-3")
    dialog.createTextField("Acceptable Mismatch Sequences: ", 4, 0, defaultText = "C>T")

    with dialog.createDynamicSelector(5, 0) as singleBaseDynSel:
        singleBaseDynSel.initCheckboxController("Output single-base mismatches")
        singleBaseSuffixDialog = singleBaseDynSel.initDisplay(True, "singleBaseSuffix")
        singleBaseSuffixDialog.createTextField("File suffix: ", 0, 0, defaultText = "_individual_mismatches")

    with dialog.createDynamicSelector(6, 0) as tandemDynSel:
        tandemDynSel.initCheckboxController("Output tandem mismatches")
        tandemSuffixDialog = tandemDynSel.initDisplay(True, "tandemSuffix")
        tandemSuffixDialog.createTextField("File suffix: ", 0, 0, defaultText = "_tandem_mismatches")

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections = dialog.selections
    
    samBedFilePaths = selections.getFilePathGroups()[0]
    
    acceptableReadLengths = parseToIterable(selections.getTextEntries()[0])
    maxMismatches = int(selections.getTextEntries()[1])
    acceptableMismatchPositions = parseToIterable(selections.getTextEntries()[2],rangeChar = '$')
    acceptableMismatchSequences = [sequence.strip() for sequence in selections.getTextEntries()[3].split(',')]

    if singleBaseDynSel.getControllerVar():
        outputSingleBaseMismatches = True
        singleBaseSuffix = selections.getTextEntries("singleBaseSuffix")[0]
    else: 
        outputSingleBaseMismatches = False
        singleBaseSuffix = ''

    if tandemDynSel.getControllerVar():
        outputTandemMismatches = True
        tandemSuffix = selections.getTextEntries("tandemSuffix")[0]
    else: 
        outputTandemMismatches = False
        tandemSuffix = ''

    samBedToCustomBed(samBedFilePaths, 
                      outputSingleBaseMismatches, singleBaseSuffix,
                      outputTandemMismatches, tandemSuffix,
                      acceptableReadLengths, maxMismatches, acceptableMismatchPositions, 
                      acceptableMismatchSequences)


if __name__ == "__main__": main()