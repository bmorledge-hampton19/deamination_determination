# Filters a simple text file where each line contains a nucleotide sequence.
# Sequences with nonstandard nucleotides (i.e. any except A, C, G, T, and N) are removed from the output.
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
nucleotides = {'A','C','G','T','N'}


def filterNonstandardNucleotides(inputFilePaths: List[str]):

    outputFilePaths: List[str] = list()

    for inputFilePath in inputFilePaths:

        print(f"Working with {os.path.basename(inputFilePath)}...")

        outputFilePath = inputFilePath.rsplit('.', 1)[0] + "_nonstandard_nuc_filtered.txt"
        outputFilePaths.append(outputFilePath)

        with open(inputFilePath, 'r') as inputFile:
            with open(outputFilePath, 'w') as outputFile:
                for line in inputFile:
                    if nucleotides.issuperset(line.strip()): outputFile.write(line)
                    else: print(line.strip())

    return outputFilePaths


def main():

    with TkinterDialog(workingDirectory = os.path.join(__file__,".."), title = "Filter Non-standard Nucleotides") as dialog:
        dialog.createMultipleFileSelector("Simple text nucleotide files:", 0, "simple_nuc.txt", ("text files", ".txt"))

    filterNonstandardNucleotides(dialog.selections.getFilePathGroups()[0])


if __name__ == "__main__": main()