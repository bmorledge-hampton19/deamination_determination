# This script searches through a directory for "tandem mismatches by read" files and their
# "mismatches by read" counterparts and compares the number of tandem mismatches (of a specific type, 
# usually CC>TT) in each.

import os
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CustomErrors import UserInputError
from benbiohelpers.FileSystemHandling.GetLineCount import getLineCount


def compareTandemMismatchCounts(mismatchesByReadDir: str, referenceSequence, mismatchSequence):

    # Make sure the reference and mismatch sequences are the same size.
    if not len(referenceSequence) == len(mismatchSequence):
        raise UserInputError("Reference and mismatch sequences are not of equal length.")

    # Create and open an output file for the results
    resultsFilePath = os.path.join(mismatchesByReadDir,
                                   f"{referenceSequence}_to_{mismatchSequence}_tandem_mismatch_counts_comparison.tsv")
    with open(resultsFilePath, 'w') as resultsFile:
        # Write the header
        resultsFile.write('\t'.join(("Data_Group_Name", "Basic_Mismatches_File_Count",
                                     "Tandem_Mismatches_File_Count", "Proportion_Missed")) + '\n')

        # Iterate through the given directory lookingn for pairs of mismatches by read files to count tandem mismatches in.
        for item in os.listdir(mismatchesByReadDir):

            # Look for tandem mismatch files from combined repititions first.
            if item.endswith(f"all_reps_{referenceSequence}_to_{mismatchSequence}_mismatches_by_read.bed"):

                # Derive the basic non-tandem mismatches file from the tandem mismatches file.
                tandemMismatchesFilePath = os.path.join(mismatchesByReadDir, item)
                basicMismatchesFilePath = tandemMismatchesFilePath.replace(f"{referenceSequence}_to_{mismatchSequence}_",'')
                dataGroupName = os.path.basename(basicMismatchesFilePath).split(f"{referenceSequence}_to_{mismatchSequence}")[0]
                if not os.path.exists(basicMismatchesFilePath):
                    raise UserInputError(f"Companion file {basicMismatchesFilePath} does not exist for "
                                        f"{tandemMismatchesFilePath}.")
                print(f"Working with {os.path.basename(basicMismatchesFilePath)} and "
                      f"{os.path.basename(tandemMismatchesFilePath)}...")

                # First, count TRUE tandem mismatches in the basic mismatches file.
                basicMismatchesTandemCount = 0
                with open(basicMismatchesFilePath, 'r') as basicMismatchesFile:
                    for line in basicMismatchesFile:
                        splitLine = line.split()
                        mismatchPositions = splitLine[3].split(':')
                        mismatchTypes = splitLine[4].split(':')

                        # Make sure we have exactly two mismatches
                        if len(mismatchPositions) != len(referenceSequence): continue
                        # Make sure those mismatches are adjacent
                        if int(mismatchPositions[0]) + 1 != int(mismatchPositions[1]): continue
                        # Make sure the mismatch type matches the given reference and mismatch sequences.
                        if (mismatchTypes[0][0]+mismatchTypes[1][0]+mismatchTypes[0][2]+mismatchTypes[1][2] !=
                            referenceSequence+mismatchSequence): continue
                        
                        # If all the above checks passed, count it!
                        basicMismatchesTandemCount += 1

                tandemMismatchesTandemCount = getLineCount(tandemMismatchesFilePath)
                missedTandemsProportion = 1 - basicMismatchesTandemCount / tandemMismatchesTandemCount
                resultsFile.write('\t'.join((dataGroupName, str(basicMismatchesTandemCount),
                                             str(tandemMismatchesTandemCount), f"{missedTandemsProportion:.4f}")) + '\n')


def main():

    # Get the working directory from mutperiod if possible. Otherwise, just use this script's directory.
    try:
        from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
        workingDirectory = getDataDirectory()
    except ImportError:
        workingDirectory = os.path.dirname(__file__)

    #Create the Tkinter UI
    with TkinterDialog(workingDirectory=workingDirectory, title = "Compare Tandem Mismatch Counts") as dialog:
        dialog.createFileSelector("Parent Directory:", 0, directory = True)
        dialog.createTextField("Reference Sequence: ", 1, 0, defaultText = "CC")
        dialog.createTextField("Mismatched Sequence: ", 2, 0, defaultText = "TT")

    selections = dialog.selections
    compareTandemMismatchCounts(selections.getIndividualFilePaths()[0], *selections.getTextEntries())


if __name__ == "__main__": main()