# This script takes a file of mismatch positions and a file of transcription factor binding sites (must both be sorted)
# and determines where those mismatches are relative to the TFBSs, within a given radius.

import os, subprocess
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import ENCOMPASSING_DATA, ENCOMPASSED_DATA
from benbiohelpers.CountThisInThat.SupplementalInformation import SimpleColumnSupInfoHandler
from benbiohelpers.InputParsing.CheckForNumber import checkForNumber
from benbiohelpers.InputParsing.ParseToIterable import parseToIterable
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs, getTempDir

 
class MismatchesAroundTFBSCounter(ThisInThatCounter):

    def setupOutputDataStratifiers(self):
        self.outputDataHandler.addEncompassedFeatureStratifier()
        self.outputDataHandler.addPlaceholderStratifier()
        self.outputDataHandler.addCustomSupplementalInformationHandler(SimpleColumnSupInfoHandler(
            outputName = "Motif_Midpoint", relevantData = ENCOMPASSING_DATA, dataCol = 1, removeDups = False
        ))
        self.outputDataHandler.addCustomSupplementalInformationHandler(SimpleColumnSupInfoHandler(
            outputName = "TFBS_Strand", relevantData = ENCOMPASSING_DATA, dataCol = 5, removeDups = False
        ))

    def setupOutputDataWriter(self):
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath, omitFinalStratificationCounts = True)


def filterAndFormatMismatches(mismatchesByReadFilePath: str, zScoresFilePath, zScoreCutoff = 4,
                              acceptableReadLengths = range(22,31), acceptableMismatchTypes = ("C>T", "CC>TT"),
                              outputSuffix = "_TFBS_formatted", addMismatchTypesToName = False):
    """
    Filters mismatches by read to acceptable mismatch patterns in reads of valid length at positions with zscores greater than or equal to the given z-score cutoff.
    The bed entry is also formatted so that it is focused directly on the mismatched position(s) and read sequence is replaced with its length.
    The resulting file is sorted, first by chromosome, then by start position and end position.
    The output suffix is appended to the new file path just after "mismatches_by_read".
    Returns the new file path.
    """

    # Determine which positions are valid for each read length based on the given z-score cutoff.
    # Initialize a nested dictionary with one key for read length, another key for position and a final boolean value stating whether the zscore meets the cutoff.
    print("Determining positions with acceptable z-scores...")
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
    outputDir = os.path.join(os.path.dirname(mismatchesByReadFilePath), ".tmp")
    checkDirs(outputDir)

    outputBaseName = (os.path.basename(mismatchesByReadFilePath).rsplit("mismatches_by_read",1)[0] +
                      "mismatches_by_read" + outputSuffix +
                      os.path.basename(mismatchesByReadFilePath).rsplit("mismatches_by_read",1)[1])
    if addMismatchTypesToName:

        for mismatchType in acceptableMismatchTypes:
            mismatchName = mismatchType.replace('>', "_to_")
            if mismatchName in outputBaseName: pass
            else: outputBaseName = outputBaseName.replace("mismatches_by_read", f"{mismatchName}_mismatches_by_read")

    outputFilePath = os.path.join(outputDir, outputBaseName)

    # Iterate through the mismatches by read file, filtering and formatting each line as necessary.
    with open(mismatchesByReadFilePath, 'r') as mismatchesByReadFile, open(outputFilePath, 'w') as filteredAndFormattedFile:

        print("Filtering and formatting mismatches...")

        for line in mismatchesByReadFile:

            # Store relevant variables and filter lines.
            splitLine = line.split()

            readLength = int(splitLine[2]) - int(splitLine[1])
            if readLength not in acceptableReadLengths: continue

            if splitLine[4] not in acceptableMismatchTypes: continue

            position = float(splitLine[3])
            if not zscoreDictionary[readLength][position]: continue

            # Recenter the coordinates on the mismatch (0-based)
            if splitLine[5] == '+': mismatchCoordinate = int(splitLine[2]) + position
            else: mismatchCoordinate = int(splitLine[1]) - position - 1
            if int(mismatchCoordinate) == mismatchCoordinate: mismatchCoordinate = int(mismatchCoordinate)

            # Write the formatted line.
            filteredAndFormattedFile.write('\t'.join((splitLine[0], str(mismatchCoordinate), str(mismatchCoordinate+1),
                                                      splitLine[3], splitLine[4], splitLine[5], str(readLength))) + '\n')
    
    # Sort the output
    print("Sorting output...")
    subprocess.check_call(("sort", "-k1,1", "-k2,2n", "-k3,3n", "-s", "-o", outputFilePath, outputFilePath))
    return outputFilePath

def relateMismatchesToTFBS(formattedMismatchesByReadFilePath: str, TFBS_FilePath: str,
                           originalMismatchesByReadFilePath: str = None, outputSuffix = "_TFBS_related_mismatches", midpointRadius = 100):
    """
    Relates mismatches to nearby transcription factor binding sites, recording the position of the TFBS midpoint as well as the strand it's on.
    Only mismatches that are nearby a TFBS, as defined by the midpointRadius argument, are tracked.
    Expects input in the format produced by the above filterAndFormatMismatches function.
    NOTE: The midpoint is assumed to be given directly by TFBS bed file. This bed file should NOT give coordinates for the full binding motif/site.
    """

    # Create the output file paths.
    
    if originalMismatchesByReadFilePath is None:
        outputBaseName = os.path.basename(formattedMismatchesByReadFilePath).rsplit('.',1)[0] + outputSuffix + ".bed"
        outputFilePath = os.path.join(os.path.dirname(os.path.dirname(formattedMismatchesByReadFilePath)), outputBaseName)
        intermediateOutputFilePath = os.path.join(getTempDir(formattedMismatchesByReadFilePath),
                                                  outputBaseName.rsplit('.', 1)[0] + "_intermediate.bed")
    else:
        outputBaseName = os.path.basename(originalMismatchesByReadFilePath).rsplit('.',1)[0] + outputSuffix + ".bed"
        outputFilePath = os.path.join(os.path.dirname(os.path.dirname(originalMismatchesByReadFilePath)), outputBaseName)
        intermediateOutputFilePath = os.path.join(getTempDir(originalMismatchesByReadFilePath),
                                                  outputBaseName.rsplit('.', 1)[0] + "_intermediate.bed")
    
    # Create and execute the counter.
    print("Associating mismatches with TFBSs...")
    counter = MismatchesAroundTFBSCounter(formattedMismatchesByReadFilePath, TFBS_FilePath, intermediateOutputFilePath,
                                          encompassingFeatureExtraRadius = midpointRadius, writeIncrementally = ENCOMPASSED_DATA,
                                          suppressOutput = True)
    counter.count()

    # Split up entries where a mismatch aligned to multiple TFBSs
    with open(intermediateOutputFilePath, 'r') as intermediateOutputFile, open(outputFilePath, 'w') as outputFile:
        for line in intermediateOutputFile:
            splitLine = line.strip().split('\t')
            for tfbs, strand in zip(splitLine[7].split(';'), splitLine[8].split(';')):
                outputFile.write('\t'.join(splitLine[:7]+[tfbs]+[strand])+'\n')


def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory=os.path.abspath(os.path.join(__file__,"..","..","data")),
                       title = "Relate Mismatches to TFBS") as dialog:
        dialog.createFileSelector("Mismatches by read:", 0, ("Bed file", ".bed"))
        dialog.createFileSelector("Z-scores file:", 1, ("Tab-separated values file", ".tsv"))
        dialog.createTextField("Z-score cutoff:", 2, 0, defaultText = "4")
        dialog.createTextField("Acceptable read lengths:", 3, 0, defaultText="22-30")
        dialog.createFileSelector("TFBS file:", 4, ("Bed file", ".bed"))
        dialog.createTextField("Output file suffix", 5, 0, defaultText = "_TFBS_related_mismatches")

    # Get the user's input from the dialog.
    selections = dialog.selections

    formattedMismatchesByReadFilePath = filterAndFormatMismatches(selections.getIndividualFilePaths()[0], selections.getIndividualFilePaths()[1],
                                                                  checkForNumber(selections.getTextEntries()[0], True),
                                                                  parseToIterable(selections.getTextEntries()[1], castType = int))
    relateMismatchesToTFBS(formattedMismatchesByReadFilePath, selections.getIndividualFilePaths()[2],
                           selections.getIndividualFilePaths()[0], selections.getTextEntries()[2])


if __name__ == "__main__": main()