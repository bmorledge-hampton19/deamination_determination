# This script filters out reads based on the sequence after the cut site in their expanded context.
# For example, reads with a 3' post-cut-site TGG could be filtered out.
# Right now, this script also requires that a companion unexpanded file be present in 
# the same directory, which will be filtered in tandem. This may become optional in the
# future if the need arises.

import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.InputParsing.ParseToIterable import parseToIterable
from benbiohelpers.CustomErrors import UserInputError


def filterOnExpandedContext(expandedContextFilePaths: List[str], threePrimeBlacklist = list(), fivePrimeBlacklist = list(),
                            customSuffix = '', filterUnexpandedCompanionFiles = True):
    
    for expandedContextFilePath in expandedContextFilePaths:
        expandedBasename = os.path.basename(expandedContextFilePath)
        unexpandedBasename, expansionContext = expandedBasename.rsplit("bp_expanded", 1)[0].rsplit('_', 1)
        outputDirectory = os.path.dirname(expandedContextFilePath)
        print(f"\nWorking with {expandedBasename}")

        # Generate the relevant output file paths, finding the companion unexpanded filepath if necessary.
        # Also, perform some quality control checks.
        if not "bp_expanded" in expandedBasename:
            raise UserInputError(f"Given file: {expandedContextFilePath} does not appear to have expanded context.")

        expandedContextOutputFilePath = os.path.join(outputDirectory, expandedBasename.rsplit('.',1)[0] + customSuffix + ".bed")
        if filterUnexpandedCompanionFiles:
            unexpandedBasename += ".bed"
            unexpandedFilePath = os.path.join(os.path.dirname(expandedContextFilePath), unexpandedBasename)
            if not os.path.exists(unexpandedFilePath):
                raise UserInputError(f"Could not find unexpanded companion file: {unexpandedFilePath}")
            unexpandedOutputFilePath = os.path.join(outputDirectory, unexpandedBasename.rsplit('.',1)[0] + customSuffix + ".bed")

        expansionContext = int(expansionContext)
        maxBlacklistSequenceLength = max(len(sequence) for sequence in threePrimeBlacklist + fivePrimeBlacklist)
        if expansionContext < maxBlacklistSequenceLength:
            raise UserInputError(f"Expansion context is insufficient for blacklisted "
                                 f"sequence length: {maxBlacklistSequenceLength}")
        

        # Read through the files in tandem, searching for lines with the blacklisted sequences to be filtered and writing
        # the rest to the new output files.
        print("Searching for blacklisted sequences just outside cut-sites...")

        with open(expandedContextFilePath, 'r') as expandedContextFile:
            with open(expandedContextOutputFilePath, 'w') as expandedContextOutputFile:
                if filterUnexpandedCompanionFiles:
                    unexpandedFile = open(unexpandedFilePath, 'r')
                    unexpandedOutputFile = open(unexpandedOutputFilePath, 'w')
                
                for expandedLine in expandedContextFile:
                    splitExpandedLine = expandedLine.split()

                    # If filtering on unexpanded companion files as well, read in the next unexpanded line and ensure
                    # that both lines are referring to the same genomic location, as its possible
                    # that some lines were filtered out during the expansion process.
                    if filterUnexpandedCompanionFiles:
                        unexpandedLine = unexpandedFile.readline()
                        splitUnexpandedLine = unexpandedLine.split()

                        while splitExpandedLine[0:3] != splitUnexpandedLine[0:3]:
                            unexpandedLine = unexpandedFile.readline()
                            assert unexpandedLine, f"Reached EOF without finding match for {expandedLine}"
                            splitUnexpandedLine = unexpandedLine.split()

                    # Check for any blacklisted sequences.
                    expandedReadSequence = splitExpandedLine[6]
                    blacklist = False
                    for blacklistSequence in threePrimeBlacklist:
                        if len(blacklistSequence) == expansionContext: 
                            slicedReadSequence = expandedReadSequence[-expansionContext:]
                        else: 
                            slicedReadSequence =expandedReadSequence[-expansionContext:
                                                                        -expansionContext + len(blacklistSequence)]
                        if blacklistSequence == slicedReadSequence:
                            blacklist = True
                            break

                    if not blacklist:
                        for blacklistSequence in fivePrimeBlacklist:
                            slicedReadSequence = expandedReadSequence[expansionContext - len(blacklistSequence):
                                                                        expansionContext]
                            if blacklistSequence == slicedReadSequence:
                                blacklist = True
                                break

                    # If the sequence was not blacklisted, rewrite the relevant lines to their respective output files.
                    if not blacklist:
                        expandedContextOutputFile.write(expandedLine)
                        if filterUnexpandedCompanionFiles: unexpandedOutputFile.write(unexpandedLine)

                if filterUnexpandedCompanionFiles:
                    unexpandedFile.close()
                    unexpandedOutputFile.close()


def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory=os.path.join(__file__,"..","..")) as dialog:
        dialog.createMultipleFileSelector("Expanded context files:",0,"bp_expanded.bed",("Bed files",".bed"))
        dialog.createCheckbox("Filter companion unexpanded files:", 1, 0)
        dialog.createTextField("3' Blacklist sequences:", 2, 0, defaultText="TGG")
        dialog.createTextField("5' Blacklist sequences:", 3, 0, defaultText='')
        dialog.createTextField("Custom suffix:", 4, 0, defaultText = "_TGG_filtered")

    # Get the user's input from the dialog.
    selections = dialog.selections

    threePrimeBlacklist = parseToIterable(selections.getTextEntries()[0])
    fivePrimeBlacklist = parseToIterable(selections.getTextEntries()[1])

    filterOnExpandedContext(selections.getFilePathGroups()[0], threePrimeBlacklist, fivePrimeBlacklist,
                            selections.getTextEntries()[2], selections.getToggleStates()[0])


if __name__ == "__main__": main()