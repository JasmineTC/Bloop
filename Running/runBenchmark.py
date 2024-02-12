import numpy as np
from datetime import date

today = date.today()

import ThreeHiggs

from ThreeHiggs import GenericModel
from ThreeHiggs import TransitionFinder

import Benchmarks.Benchmarks_3HDM
from PlotResult import PlotResult

from ThreeHiggs import MinimizationAlgos

userinput = ThreeHiggs.UserInput()
args = userinput.parse()

hardToSoftFile = ThreeHiggs.getResourcePath("Data/HardToSoft/softScaleParams_NLO.txt")
softScaleRGEFile = ThreeHiggs.getResourcePath("Data/HardToSoft/softScaleRGE.txt")
softToUltrasoftFile = ThreeHiggs.getResourcePath("Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt")


## Model object setup + load matching relations
model3HDM = GenericModel(loopOrder = args.loopOrder)
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile, softScaleRGEFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)

## Set algorithm to use for Veff minimization
model3HDM.effectivePotential.minimizer.setAlgorithm(MinimizationAlgos.eDIRECTGLOBAL)
## Set tolerances used by global and local methods in Veff minimization
## Order is global abs, global rel, local abs, local rel
model3HDM.effectivePotential.minimizer.setTolerances(1e-1, 1e-1, 1e-5, 1e-5)


print("!!!")
print("Currently not matching soft --> ultrasoft, this is WIP. Also: 2-loop masses lack some log terms")
print("!!!")


inputParams = Benchmarks.Benchmarks_3HDM.bmList[args.benchMarkNumber]

transitionFinder = TransitionFinder(model=model3HDM)
model3HDM.setInputParams(inputParams)
minimizationResults = transitionFinder.traceFreeEnergyMinimum()

filename = f"Results/{date.today()}-BM-{args.benchMarkNumber}-LoopOrder{args.loopOrder}"
np.savetxt(filename + ".txt", minimizationResults)
if args.plot == True:
    PlotResult.PlotData(minimizationResults, args.benchMarkNumber, filename)
