# This script takes a file of mismatch positions and a file of genomic features, such as transcription factor binding sites (must both be sorted)
# and determines where those mismatches are relative to the genomic features, within a given radius.

import os
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import ENCOMPASSING_DATA, ENCOMPASSED_DATA
from benbiohelpers.CountThisInThat.SupplementalInformation import SimpleColumnSupInfoHandler
from benbiohelpers.InputParsing.CheckForNumber import checkForNumber
from benbiohelpers.InputParsing.ParseToIterable import parseToIterable
from benbiohelpers.FileSystemHandling.DirectoryHandling import getTempDir

 
class MismatchesAroundFeatureCounter(ThisInThatCounter):

    def setupOutputDataStratifiers(self):
        self.outputDataHandler.addEncompassedFeatureStratifier()
        self.outputDataHandler.addPlaceholderStratifier()
        self.outputDataHandler.addCustomSupplementalInformationHandler(SimpleColumnSupInfoHandler(
            outputName = "Teature_Midpoint", relevantData = ENCOMPASSING_DATA, dataCol = 1, removeDups = False
        ))
        self.outputDataHandler.addCustomSupplementalInformationHandler(SimpleColumnSupInfoHandler(
            outputName = "Feature_Strand", relevantData = ENCOMPASSING_DATA, dataCol = 5, removeDups = False
        ))

    def setupOutputDataWriter(self):
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath, omitFinalStratificationCounts = True)


def relateMismatchesToFeature(formattedMismatchesByReadFilePath: str, featureFilePath: str,
                              outputSuffix = "_feature_related", midpointRadius = 100, verbose = False,
                              enforcedStrand = None):
    """
    Relates mismatches to nearby genomic features, recording the position of the feature midpoint as well as the strand it's on.
    Only mismatches that are nearby a feature, as defined by the midpointRadius argument, are tracked.
    Expects input in the format produced by the formatMismatchesForRelation function.
    NOTE: It is assumed that the midpoint is given directly by the genomic feature bed file. This bed file should NOT give coordinates for the full feature.
    """

    # Create the output file paths.
    outputBaseName = os.path.basename(formattedMismatchesByReadFilePath).rsplit('.',1)[0] + outputSuffix + ".bed"
    if os.path.dirname(formattedMismatchesByReadFilePath).endswith(".tmp"):
        outputFilePath = os.path.join(os.path.dirname(os.path.dirname(formattedMismatchesByReadFilePath)), outputBaseName)
    else:
        outputFilePath = os.path.join(os.path.dirname(formattedMismatchesByReadFilePath), outputBaseName)
    intermediateOutputFilePath = os.path.join(getTempDir(formattedMismatchesByReadFilePath),
                                                outputBaseName.rsplit('.', 1)[0] + "_intermediate.bed")
    
    # Create and execute the counter.
    if verbose: print("Associating mismatches with genomic feature...")
    counter = MismatchesAroundFeatureCounter(formattedMismatchesByReadFilePath, featureFilePath, intermediateOutputFilePath,
                                          encompassingFeatureExtraRadius = midpointRadius, writeIncrementally = ENCOMPASSED_DATA,
                                          suppressOutput = True)
    counter.count()

    # Split up entries where a mismatch aligned to multiple features
    with open(intermediateOutputFilePath, 'r') as intermediateOutputFile, open(outputFilePath, 'w') as outputFile:
        for line in intermediateOutputFile:
            splitLine = line.strip().split('\t')
            if enforcedStrand is None:
                for featureMidpoint, strand in zip(splitLine[7].split(';'), splitLine[8].split(';')):
                    outputFile.write('\t'.join(splitLine[:7]+[featureMidpoint]+[strand])+'\n')
            else:
                for featureMidpoint in splitLine[7].split(';'):
                    outputFile.write('\t'.join(splitLine[:7]+[featureMidpoint]+[enforcedStrand])+'\n')


def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory=os.path.abspath(os.path.join(__file__,"..","..","data")),
                       title = "Relate Mismatches to Feature") as dialog:
        dialog.createFileSelector("Formatted mismatches by read:", 0, ("Bed file", ".bed"))
        dialog.createFileSelector("Feature file:", 1, ("Bed file", ".bed"))
        dialog.createTextField("Output file suffix", 2, 0, defaultText = "_feature_related")
        dialog.createTextField("Midpoint radius:", 3, 0, defaultText = "100")

    # Get the user's input from the dialog.
    selections = dialog.selections

    relateMismatchesToFeature(selections.getIndividualFilePaths()[0], selections.getIndividualFilePaths()[1],
                              selections.getTextEntries()[0], selections.getTextEntries()[1])


if __name__ == "__main__": main()