import pathlib
import numpy as np
from GenericModel import GenericModel
from TransitionFinder import TransitionFinder

import Benchmarks.Benchmarks_3HDM
from Userinput import Userinput
from PlotResult import PlotResult

userinput = Userinput()
args = userinput.parse()

## hack file path, not a good solution if making this into a proper package

pathToCurrentFile = pathlib.Path(__file__).parent.resolve()

hardToSoftFile = str(pathToCurrentFile) + "/Data/HardToSoft/softScaleParams_NLO.txt"
softToUltrasoftFile = str(pathToCurrentFile) + "/Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt"

## Model object setup + load matching relations
model3HDM = GenericModel(loopOrder = args.loopOrder)
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)

print("!!! \n Currently not matching soft --> ultrasoft, this is WIP. Also: 2-loop masses lack some log terms  \n !!!")

print("Start finite-T stuff")

inputParams = Benchmarks.Benchmarks_3HDM.bmList[args.benchMarkNumber]

transitionFinder = TransitionFinder(model=model3HDM)
model3HDM.setInputParams(inputParams)
minimizationResults = transitionFinder.traceFreeEnergyMinimum()


filename = "Data/Results/bm" + str(args.benchMarkNumber) + ".txt"
np.savetxt(filename, minimizationResults)
if args.plot == True:
    PlotResult.PlotData(minimizationResults, args.benchMarkNumber)
