# This script takes "mismatches by read" files and filters/formats them for use with the RelateMismatchesToFeatures.py script.
import os, subprocess, shutil
from typing import List
from benbiohelpers.FileSystemHandling.DirectoryHandling import getTempDir
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.InputParsing.CheckForNumber import checkForNumber
from benbiohelpers.InputParsing.ParseToIterable import parseToIterable
from benbiohelpers.CustomErrors import UserInputError


def formatMismatchesForRelation(mismatchesByReadFilePath: str, zScoresFilePath, zScoreCutoff = 4,
                                acceptableReadLengths = range(22,31),acceptableMismatchTypes: List[str] = ["C>T", "CC>TT"],
                                unacceptableMismatchTypes: List[str] = None, filterCompositeMismatches = True,
                                outputSuffix = "_relation_formatted", addMismatchTypesToName = True, outputToTmpDir = False,
                                sortOutput = True, verbose = False):
    """
    Filters mismatches by read to acceptable mismatch patterns in reads of valid length at positions with zscores greater than or equal to the given z-score cutoff.
    Also filter reads with sequences containing 'N' nucleotides (to keep things consistent with R pipeline).
    The bed entry is also formatted so that it is focused directly on the mismatched position(s) and read sequence is replaced with its length.
    The resulting file is sorted, first by chromosome, then by start position and end position.
    The output suffix is appended to the new file path just after "mismatches_by_read".
    Returns the new file path.
    """

    # Make sure invalid arguments were not given.
    if acceptableMismatchTypes is not None and unacceptableMismatchTypes is not None:
        raise UserInputError("Acceptable mismatch types and unacceptable mismatch types are incompatible and cannot both have values.")


    # Determine which positions are valid for each read length based on the given z-score cutoff.
    if zScoresFilePath is not None:
        # Initialize a nested dictionary with one key for read length, another key for position and a final boolean value stating whether the zscore meets the cutoff.
        if verbose: print("Determining positions with acceptable z-scores...")
        zscoreDictionary = dict()
        for readLength in acceptableReadLengths:
            zscoreDictionary[readLength] = dict()

        # Determine which positions have zscores which meet the cutoff.
        with open(zScoresFilePath, 'r') as zScoresFile:
            zScoresFile.readline() # Skip the headers
            for line in zScoresFile:
                readLength, position, frequency, zScore = line.split()
                readLength = int(readLength)
                if readLength not in acceptableReadLengths: continue
                else: zscoreDictionary[readLength][float(position)] = float(zScore) >= zScoreCutoff

    # Generate the output file path.
    if outputToTmpDir: outputDir = getTempDir(mismatchesByReadFilePath)
    else: outputDir = os.path.dirname(mismatchesByReadFilePath)

    outputBaseName = (os.path.basename(mismatchesByReadFilePath).rsplit("mismatches_by_read",1)[0] +
                      "mismatches_by_read" + outputSuffix +
                      os.path.basename(mismatchesByReadFilePath).rsplit("mismatches_by_read",1)[1])
    if addMismatchTypesToName:

        if acceptableMismatchTypes is not None:
            for mismatchType in acceptableMismatchTypes:
                mismatchName = mismatchType.replace('>', "_to_")
                if mismatchName in outputBaseName: pass
                else: outputBaseName = outputBaseName.replace("mismatches_by_read", f"{mismatchName}_mismatches_by_read")

        elif unacceptableMismatchTypes is not None:
            for mismatchType in unacceptableMismatchTypes:
                if "omitted_mismatches_by_read" in outputBaseName: pass
                else: outputBaseName = outputBaseName.replace("mismatches_by_read", "omitted_mismatches_by_read")
                mismatchName = mismatchType.replace('>', "_to_")
                if mismatchName in outputBaseName: pass
                else: outputBaseName = outputBaseName.replace("omitted_mismatches_by_read", f"{mismatchName}_omitted_mismatches_by_read")

    outputFilePath = os.path.join(outputDir, outputBaseName)

    # Iterate through the mismatches by read file, filtering and formatting each line as necessary.
    with open(mismatchesByReadFilePath, 'r') as mismatchesByReadFile, open(outputFilePath, 'w') as filteredAndFormattedFile:

        if verbose: print("Filtering and formatting mismatches...")

        for line in mismatchesByReadFile:

            # Store relevant variables and filter lines.
            splitLine = line.split()

            readLength = int(splitLine[2]) - int(splitLine[1])
            if readLength not in acceptableReadLengths: continue

            if filterCompositeMismatches and ':' in splitLine[4]: continue
            if acceptableMismatchTypes is not None and splitLine[4] not in acceptableMismatchTypes: continue
            if unacceptableMismatchTypes is not None and splitLine[4] in unacceptableMismatchTypes: continue

            position = float(splitLine[3])
            if zScoresFilePath is not None and not zscoreDictionary[readLength][position]: continue

            if 'N' in splitLine[6]: continue

            # Recenter the coordinates on the mismatch (0-based)
            if splitLine[5] == '+': mismatchCoordinate = int(splitLine[2]) + position
            else: mismatchCoordinate = int(splitLine[1]) - position - 1
            if int(mismatchCoordinate) == mismatchCoordinate: mismatchCoordinate = int(mismatchCoordinate)

            # Write the formatted line.
            filteredAndFormattedFile.write('\t'.join((splitLine[0], str(mismatchCoordinate), str(mismatchCoordinate+1),
                                                      splitLine[3], splitLine[4], splitLine[5], str(readLength))) + '\n')

    # Sort the output
    if sortOutput:
        if verbose: print("Sorting output...")
        subprocess.check_call(("sort", "-k1,1", "-k2,2n", "-k3,3n", "-s", "-o", outputFilePath, outputFilePath))

    return outputFilePath


def combineSingleAndTandemMismatches(singleMismatchFilePath: str, tandemMismatchFilePath: str,
                                     singleMismatchTypes: List[str] = ["C>T"], tandemMismatchTypes: List[str] = ["CC>TT"]):
    """
    Combines two mismatch files, one with single-base mismatches, and one with tandem (two-base) mismatches.
    This function is meant to work specifically with the output from the formatMismatchesForRelation function.
    The resulting file is named using the single-mismatch file path as a base, with both mismatch types included.
    Returns the output file path
    NOTE: It might be a good idea to make a variation of this into its own function in benbiohelpers at some point.
    """

    # Create the output file path
    singleMismatchTypes = '_'.join([singleMismatchType.replace('>', "_to_") for singleMismatchType in singleMismatchTypes])
    tandemMismatchTypes = '_'.join([tandemMismatchType.replace('>', "_to_") for tandemMismatchType in tandemMismatchTypes])
    outputFilePath = singleMismatchFilePath.replace(singleMismatchTypes, singleMismatchTypes + '_' + tandemMismatchTypes)

    # Copy over the contents of the single mismatch file to the output file. Then, split up and append the tandem mismatches.
    shutil.copy(singleMismatchFilePath, outputFilePath)
    with open(tandemMismatchFilePath, 'r') as tandemMismatchFile, open(outputFilePath, 'a') as outputFile:
        
        for line in tandemMismatchFile:

            splitLine = line.strip().split('\t')

            # Format for the first base in the lesion
            splitLine[1] = splitLine[1][:-2]
            splitLine[2] = splitLine[2][:-2]
            splitLine[3] = splitLine[3][:-2]
            outputFile.write('\t'.join(splitLine) + '\n')

            # Format for the second base in the lesion
            splitLine[1] = str(int(splitLine[1]) + 1)
            splitLine[2] = str(int(splitLine[2]) + 1)
            splitLine[3] = str(int(splitLine[3]) - 1)
            outputFile.write('\t'.join(splitLine) + '\n')


    subprocess.check_call(("sort", "-k1,1", "-k2,2n", "-k3,3n", "-s", "-o", outputFilePath, outputFilePath))
    return outputFilePath


def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory=os.path.abspath(os.path.join(__file__,"..","..","data")),
                       title = "Relate Mismatches to Feature") as dialog:
        dialog.createFileSelector("Mismatches by read:", 0, ("Bed file", ".bed"))
        dialog.createTextField("Acceptable Mismatch Types:", 1, 0, defaultText = "C>T")
        dialog.createFileSelector("Z-scores file:", 2, ("Tab-separated values file", ".tsv"))
        dialog.createTextField("Z-score cutoff:", 3, 0, defaultText = "4")
        dialog.createTextField("Acceptable read lengths:", 4, 0, defaultText="22-30")
        dialog.createCheckbox("Combine CC>TT tandem mismatch file:", 5, 0)
        

    # Get the user's input from the dialog.
    selections = dialog.selections

    formattedMismatchesByReadFilePath = formatMismatchesForRelation(selections.getIndividualFilePaths()[0], selections.getIndividualFilePaths()[1],
                                                                    checkForNumber(selections.getTextEntries()[1], True),
                                                                    parseToIterable(selections.getTextEntries()[2], castType = int),
                                                                    parseToIterable(selections.getTextEntries()[0]))
    
    if selections.getToggleStates()[0]:

        # NOTE: A lot of assumptions are being made here about the associated file paths.
        tandemFormattedMismatchesByReadFilePath = formatMismatchesForRelation(
            selections.getIndividualFilePaths()[0].replace("mismatches_by_read", "CC_to_TT_mismatches_by_read"),
            selections.getIndividualFilePaths()[1].replace("C_to_T_mismatch_frequency_z-scores","CC_to_TT_mismatch_frequency_z-scores"),
            checkForNumber(selections.getTextEntries()[1], True), parseToIterable(selections.getTextEntries()[2], castType = int), ["CC>TT"]
        )

        combineSingleAndTandemMismatches(formattedMismatchesByReadFilePath, tandemFormattedMismatchesByReadFilePath)


if __name__ == "__main__": main()