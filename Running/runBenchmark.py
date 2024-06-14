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

## This should be put inside the pickle
hardToSoftFile = ThreeHiggs.getResourcePath("Data/HardToSoft/softScaleParams_NLO.txt")
softScaleRGEFile = ThreeHiggs.getResourcePath("Data/HardToSoft/softScaleRGE.txt")
softToUltrasoftFile = ThreeHiggs.getResourcePath("Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt")

## Model object setup + load matching relations
model3HDM = GenericModel()
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile, softScaleRGEFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)

hardToSoftFile = ThreeHiggs.getResourcePath("Data/HardToSoft/softScaleParams_NLO.txt")
softScaleRGEFile = ThreeHiggs.getResourcePath("Data/HardToSoft/softScaleRGE.txt")
softToUltrasoftFile = ThreeHiggs.getResourcePath("Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt")

## Model object setup + load matching relations
model3HDM = GenericModel()
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile, softScaleRGEFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)
## ---- Configure Veff
veffFiles = []
veffFiles.append( ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/Veff_LO.txt") )
if (args.loopOrder >= 1):
    veffFiles.append( ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/Veff_NLO.txt") )
if (args.loopOrder >= 2):
    veffFiles.append( ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/Veff_NNLO.txt") )

## EIt's (slightly) faster to generate this config file once and then store it in binary
## Use pickle to store/laod the binary file.
veffConfig = ThreeHiggs.VeffConfig(
    fieldNames = ['v1', 'v2', 'v3'],
    loopOrder = args.loopOrder,
    veffFiles = veffFiles,
    vectorMassFile = ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/vectorMasses.txt"),
    vectorShorthandFile = ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/vectorShorthands.txt"),
    #
    scalarPermutationMatrix = ParsedMatrix.parseConstantMatrix(ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/scalarPermutationMatrix.txt")),
    scalarMassMatrices = [ 
        ThreeHiggs.MatrixDefinitionFiles(ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/scalarMassMatrix_upperLeft.txt"),
                                         ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/scalarMassMatrix_upperLeft_definitions.txt")),
        ThreeHiggs.MatrixDefinitionFiles(ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/scalarMassMatrix_bottomRight.txt"),
                                         ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/scalarMassMatrix_bottomRight_definitions.txt"))
    ],
    scalarRotationMatrixFile = ThreeHiggs.getResourcePath("Data/EffectivePotential_threeFields/scalarRotationMatrix.txt"),
    # We will take abs values of all mass^2
    bAbsoluteMsq = True,
)
model3HDM.effectivePotential.configure(veffConfig)

##TODO pickle this as best as possible
model3HDM.effectivePotential.configure(veffConfig)


## Set algorithm to use for Veff minimization
model3HDM.effectivePotential.minimizer.setAlgorithm(MinimizationAlgos.eDIRECTGLOBAL)
## Set tolerances used by global and local methods in Veff minimization
## Order is global abs, global rel, local abs, local rel
model3HDM.effectivePotential.minimizer.setTolerances(1e-1, 1e-1, 1e-5, 1e-5)
## Set benchMarkNumber in minimizer for title and file naming
model3HDM.effectivePotential.minimizer.setBmNumber(args.benchMarkNumber)


inputParams = Benchmarks.Benchmarks_3HDM.bmList[args.benchMarkNumber]
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
    
