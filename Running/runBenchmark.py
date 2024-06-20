import numpy as np
from datetime import date

import ThreeHiggs
from ThreeHiggs import GenericModel
from ThreeHiggs import TransitionFinder
from ThreeHiggs import MinimizationAlgos
from ThreeHiggs.parsedmatrix import ParsedMatrix

import Benchmarks.Benchmarks_3HDM

import pickle ##Note 

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

userinput = ThreeHiggs.UserInput()
args = userinput.parse()

## ---- Configure Veff
veffFiles = [getResourcePath(args.loFile)]
if (args.loopOrder >= 1):
    veffFiles.append( getResourcePath(args.nloFile) )
if (args.loopOrder >= 2):
    veffFiles.append( getResourcePath(args.nnloFile) )

## This should be put inside the pickle
hardToSoftFile = getResourcePath(args.hardToSoftFile)
softScaleRGEFile = getResourcePath(args.softScaleRGEFile)
softToUltrasoftFile = getResourcePath(args.softToUltraSoftFile)

from ThreeHiggs.EffectivePotential import EffectivePotential
effectivePotential = EffectivePotential(['v1', 'v2', 'v3'],
                                        True,
                                        getResourcePath(args.vectorsMassesSquaredFile),
                                        getResourcePath(args.vectorsShortHandsFile),
                                        ParsedMatrix.parseConstantMatrix(getResourcePath(args.scalarPermutationFile)),
                                        [ThreeHiggs.MatrixDefinitionFiles(getResourcePath(args.scalarMassMatrixUpperLeftFile),
                                                                          getResourcePath(args.scalarMassMatrixUpperLeftDefinitionsFile)),
                                         ThreeHiggs.MatrixDefinitionFiles(getResourcePath(args.scalarMassMatrixBottomRightFile),
                                                                          getResourcePath(args.scalarMassMatrixBottomRightDefinitionsFile))],
                                        getResourcePath(args.scalarRotationFile),
                                        args.loopOrder,
                                        veffFiles)

## Model object setup + load matching relations
model3HDM = GenericModel(effectivePotential)
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile, softScaleRGEFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)

## Set algorithm to use for Veff minimization
model3HDM.effectivePotential.minimizer.setAlgorithm(MinimizationAlgos.eDIRECTGLOBAL)
## Set tolerances used by global and local methods in Veff minimization
## Order is global abs, global rel, local abs, local rel
model3HDM.effectivePotential.minimizer.setTolerances(args.absGlobalTolerance,
                                                     args.relGlobalTolerance, 
                                                     args.absLocalTolerance, 
                                                     args.relLocalTolerance)

## Set benchMarkNumber in minimizer for title and file naming
model3HDM.effectivePotential.minimizer.setBmNumber(args.benchMarkNumber)

import json
with open(args.benchMarkFile) as benchMarkFile:
    inputParams = json.load(benchMarkFile)[args.benchMarkNumber]

ghdm = inputParams["ghDM"]
model3HDM.effectivePotential.minimizer.setgHDM(ghdm)

transitionFinder = TransitionFinder(model=model3HDM)
model3HDM.setInputParams(inputParams)
minimizationResults = transitionFinder.traceFreeEnergyMinimum()

print(f"{minimizationResults=}")

#filename = f"Results/{date.today()}-BM-{args.benchMarkNumber}-LoopOrder{args.loopOrder}"
#filename = f"Results/Debug/g_01/SS_Off/BM_{args.benchMarkNumber}_gHDM_{ghdm}_SS_Off"

filename = f"{args.resultsDirectory}/BM_{args.benchMarkNumber}"

from pathlib import Path
Path(args.resultsDirectory).mkdir(parents = True, exist_ok = True)

if args.save == True:
    np.savetxt(filename + ".txt", minimizationResults)

if args.plot == True:
    from PlotResult import PlotResult
    PlotResult.PlotData(minimizationResults, args.benchMarkNumber,args.loopOrder, filename)
    
