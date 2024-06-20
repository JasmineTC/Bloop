import numpy as np
from datetime import date

import ThreeHiggs
from ThreeHiggs import GenericModel
from ThreeHiggs import TransitionFinder
from ThreeHiggs import MinimizationAlgos
from ThreeHiggs.parsedmatrix import ParsedMatrix

import Benchmarks.Benchmarks_3HDM

import pickle ##Note 

userinput = ThreeHiggs.UserInput()
args = userinput.parse()

## ---- Configure Veff
veffFiles = [ThreeHiggs.getResourcePath(args.loFile)]
if (args.loopOrder >= 1):
    veffFiles.append( ThreeHiggs.getResourcePath(args.nloFile) )
if (args.loopOrder >= 2):
    veffFiles.append( ThreeHiggs.getResourcePath(args.nnloFile) )

## This should be put inside the pickle
hardToSoftFile = ThreeHiggs.getResourcePath(args.hardToSoftFile)
softScaleRGEFile = ThreeHiggs.getResourcePath(args.softScaleRGEFile)
softToUltrasoftFile = ThreeHiggs.getResourcePath(args.softToUltraSoftFile)

from ThreeHiggs.EffectivePotential import EffectivePotential
effectivePotential = EffectivePotential(['v1', 'v2', 'v3'],
                                        True,
                                        ThreeHiggs.getResourcePath(args.vectorsMassesSquaredFile),
                                        ThreeHiggs.getResourcePath(args.vectorsShortHandsFile),
                                        ParsedMatrix.parseConstantMatrix(ThreeHiggs.getResourcePath(args.scalarPermutationFile)),
                                        [ThreeHiggs.MatrixDefinitionFiles(ThreeHiggs.getResourcePath(args.scalarMassMatrixUpperLeftFile),
                                                                          ThreeHiggs.getResourcePath(args.scalarMassMatrixUpperLeftDefinitionsFile)),
                                         ThreeHiggs.MatrixDefinitionFiles(ThreeHiggs.getResourcePath(args.scalarMassMatrixBottomRightFile),
                                                                          ThreeHiggs.getResourcePath(args.scalarMassMatrixBottomRightDefinitionsFile))],
                                        ThreeHiggs.getResourcePath(args.scalarRotationFile),
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
    
