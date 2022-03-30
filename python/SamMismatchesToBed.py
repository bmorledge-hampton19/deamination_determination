# This script takes a sam file and locates all the mismatches (potential deamination events) in it.
# The results are written to a bed file along with a metadata file containing the total number of
# reads analyzed, and the number of mismatches found.
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog


def samMismatchesToBed(samFilePaths: List[str]):
    
    for samFilePath in samFilePaths:

        # Create output file paths (bed file + metadata)
        outputBedFilePath = samFilePath.rsplit('.',1)[0] + "_mismatches.bed"
        metadataFilePath = outputBedFilePath.rsplit('.',1)[0] + ".metadata"

        # Read through the sam file line by line, looking for mismatches and recording them.
        mismatches = 0
        with open(samFilePath, 'r') as samFile:
            for line in samFile:

                # Skip header lines.
                if line.startswith('@'): continue

                splitLine = line.split()

                # Skip lines that didn't align.
                if splitLine[2] == '*': continue

                # Using the "MD" optinal field to determine if there were any mismatches
                mismatchDesignations = splitLine[17].rsplit(':',1)[1]


def main():
    # Create the Tkinter dialog.
    dialog = TkinterDialog(workingDirectory=os.path.join(__file__,"..",".."))
    dialog.createMultipleFileSelector("Sam Read Files:",0,".sam.gz",("Sam Files",(".sam.gz",".sam")), 
                                      additionalFileEndings = [".sam"])

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections = dialog.selections


if __name__ == "__main__": main()