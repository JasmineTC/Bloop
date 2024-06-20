import ThreeHiggs

def getResourcePath(relativePathToResource: str) -> str:
    """ Gives a safe path to a packaged resource.
    
    Returns
    -------
    Path to the resource file: str.
    """

    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    from importlib.resources import files
    return files(packageName) / relativePathToResource

from ThreeHiggs.UserInput import UserInput
userinput = UserInput()
args = userinput.parse()

## ---- Configure Veff
veffFiles = [getResourcePath(args.loFile)]
if (args.loopOrder >= 1):
    veffFiles.append( getResourcePath(args.nloFile) )
if (args.loopOrder >= 2):
    veffFiles.append( getResourcePath(args.nnloFile) )

hardToSoftFile = getResourcePath(args.hardToSoftFile)
softScaleRGEFile = getResourcePath(args.softScaleRGEFile)
softToUltrasoftFile = getResourcePath(args.softToUltraSoftFile)

from ThreeHiggs.EffectivePotential import EffectivePotential
from ThreeHiggs.parsedmatrix import ParsedMatrix, MatrixDefinitionFiles
effectivePotential = EffectivePotential(['v1', 'v2', 'v3'],
                                        True,
                                        getResourcePath(args.vectorsMassesSquaredFile),
                                        getResourcePath(args.vectorsShortHandsFile),
                                        ParsedMatrix.parseConstantMatrix(getResourcePath(args.scalarPermutationFile)),
                                        [MatrixDefinitionFiles(getResourcePath(args.scalarMassMatrixUpperLeftFile),
                                                               getResourcePath(args.scalarMassMatrixUpperLeftDefinitionsFile)),
                                         MatrixDefinitionFiles(getResourcePath(args.scalarMassMatrixBottomRightFile),
                                                               getResourcePath(args.scalarMassMatrixBottomRightDefinitionsFile))],
                                        getResourcePath(args.scalarRotationFile),
                                        args.loopOrder,
                                        veffFiles)

## Model object setup + load matching relations
from ThreeHiggs.GenericModel import GenericModel
model3HDM = GenericModel(effectivePotential)
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile, softScaleRGEFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)

## Set algorithm to use for Veff minimization
from ThreeHiggs.VeffMinimizer import MinimizationAlgos
model3HDM.effectivePotential.minimizer.setAlgorithm(MinimizationAlgos.eDIRECTGLOBAL)
## Set tolerances used by global and local methods in Veff minimization
## Order is global abs, global rel, local abs, local rel
model3HDM.effectivePotential.minimizer.setTolerances(args.absGlobalTolerance,
                                                     args.relGlobalTolerance, 
                                                     args.absLocalTolerance, 
                                                     args.relLocalTolerance)

## Set benchMarkNumber in minimizer for title and file naming
model3HDM.effectivePotential.minimizer.setBmNumber(args.benchMarkNumber)

with open(args.benchMarkFile) as benchMarkFile:
    from json import load
    inputParams = load(benchMarkFile)[args.benchMarkNumber]

ghdm = inputParams["ghDM"]
model3HDM.effectivePotential.minimizer.setgHDM(ghdm)

from ThreeHiggs.TransitionFinder import TransitionFinder
transitionFinder = TransitionFinder(model=model3HDM)
model3HDM.setInputParams(inputParams)
minimizationResults = transitionFinder.traceFreeEnergyMinimum()

print(f"{minimizationResults=}")

filename = f"{args.resultsDirectory}/BM_{args.benchMarkNumber}"

from pathlib import Path
Path(args.resultsDirectory).mkdir(parents = True, exist_ok = True)

if args.save == True:
    from numpy import savetxt
    savetxt(filename + ".txt", minimizationResults)

if args.plot == True:
    from PlotResult import PlotResult
    PlotResult.PlotData(minimizationResults, args.benchMarkNumber,args.loopOrder, filename)
    
