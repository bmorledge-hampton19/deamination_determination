# Provides a python interface to my preferred weblogo command line arguments.
# Input should be a text file containing one DNA sequence per line.
import os, subprocess
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CustomErrors import InvalidPathError, UserInputError
FIVE_PRIME = "three_prime"
THREE_PRIME = "five_prime"

def createWebLogo(inputFilePaths: List[str], composition):

    outputFilePaths: List[str] = list()

    for inputFilePath in inputFilePaths:

        print(f"Working with {os.path.basename(inputFilePath)}...")

        if "three_prime" in os.path.basename(inputFilePath): strandPolarity = THREE_PRIME
        elif "five_prime" in os.path.basename(inputFilePath): strandPolarity = FIVE_PRIME
        else: raise InvalidPathError(inputFilePath, "Expected strand polarity in basename of file path")

        outputFilePath = inputFilePath.replace("sequence_logo_input","sequence_logo").replace(".txt", ".png")
        outputFilePaths.append(outputFilePath)

        webLogoArgs = ["weblogo", "-F", "png", "--errorbars", "NO", "-C", "red", "A", '"Adenine"', "-C", "blue", "C",
                       '"Cytosine"', "-C", "green", "G", '"Guanine"', "-C", "purple", "T", '"Thymine"', "--aspect-ratio",
                       '3', "--resolution", "1000", "-S", "0.5", "--annotate"]
        if strandPolarity == FIVE_PRIME: webLogoArgs.append("-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7")
        elif strandPolarity == THREE_PRIME: webLogoArgs.append("-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9") 
        
        if composition is None:
            pass
        elif type(composition) is int:
            GCContent = composition
            ATContent = 100 - GCContent
            webLogoArgs += ["--composition",
                            f"\"{{'A':{ATContent/2}, 'C':{GCContent/2}, 'G':{GCContent/2}, 'T':{ATContent/2}}}\""]
        else:
            webLogoArgs += ["--composition", f"\"{composition}\""]

        webLogoArgs += ['<', inputFilePath, '>', outputFilePath]

        subprocess.check_call(' '.join(webLogoArgs), shell = True)

    return outputFilePaths


def main():

    with TkinterDialog(workingDirectory = os.path.join(__file__,".."), title = "Create Web Logo") as dialog:
        dialog.createMultipleFileSelector("Sequence logo input files:", 0, "sequence_logo_input.txt", ("text files", ".txt"),
                                          additionalFileEndings = ["sequence_logo_input_TGG_filtered.txt"])

        with dialog.createDynamicSelector(1, 0) as compositionDynSel:
            compositionDynSel.initDropdownController("Background", ["None", "GC Content", "Species"])
            compositionDynSel.initDisplay("GC Content", "GCContent").createTextField("GC Content:", 0, 0)
            compositionDynSel.initDisplay("Species", "Species").createDropdown("Species",0,0,
                                                                               ["H. sapiens", "S. cerevisiae"])
            
    if compositionDynSel.getControllerVar() == "None":
        composition = None
    elif compositionDynSel.getControllerVar() == "GC Content":
        composition = dialog.selections.getTextEntries("GCContent")[0]
        try:
            composition = int(composition)
        except ValueError:
            raise UserInputError("GC Content must be an integer.")
    elif compositionDynSel.getControllerVar() == "Species":
        composition = dialog.selections.getDropdownSelections("Species")[0]

    createWebLogo(dialog.selections.getFilePathGroups()[0], composition)


if __name__ == "__main__": main()