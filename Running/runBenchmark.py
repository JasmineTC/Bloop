import numpy as np

import ThreeHiggs

from ThreeHiggs import GenericModel
from ThreeHiggs import TransitionFinder

import Benchmarks.Benchmarks_3HDM
from PlotResult import PlotResult

userinput = ThreeHiggs.UserInput()
args = userinput.parse()

hardToSoftFile = ThreeHiggs.getResourcePath("Data/HardToSoft/softScaleParams_NLO.txt")
softToUltrasoftFile = ThreeHiggs.getResourcePath("Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt")


## Model object setup + load matching relations
model3HDM = GenericModel(loopOrder = args.loopOrder)
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)


print("!!!")
print("Currently not matching soft --> ultrasoft, this is WIP. Also: 2-loop masses lack some log terms")
print("!!!")


inputParams = Benchmarks.Benchmarks_3HDM.bmList[args.benchMarkNumber]

transitionFinder = TransitionFinder(model=model3HDM)
model3HDM.setInputParams(inputParams)
minimizationResults = transitionFinder.traceFreeEnergyMinimum()

filename = "Results/BM" + str(args.benchMarkNumber) + ".txt"
np.savetxt(filename, minimizationResults)
if args.plot == True:
    PlotResult.PlotData(minimizationResults, args.benchMarkNumber)
