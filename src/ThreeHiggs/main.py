import numpy as np
import pathlib

from GenericModel import GenericModel
from TransitionFinder import TransitionFinder

import Benchmarks.Benchmarks_3HDM

## hack file path, not a good solution if making this into a proper package
pathToCurrentFile = pathlib.Path(__file__).parent.resolve()

hardToSoftFile = str(pathToCurrentFile) + "/Data/HardToSoft/softScaleParams_NLO.txt"
softToUltrasoftFile = str(pathToCurrentFile) + "/Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt"


## Model object setup + load matching relations
model3HDM = GenericModel()
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)


inputParams = Benchmarks.Benchmarks_3HDM.BM1
#inputParams = Benchmarks.Benchmarks_3HDM.BM_SM_like

transitionFinder = TransitionFinder(model=model3HDM)

## Scanning loops would start here

print("!!!")
print("Currently not matching soft --> ultrasoft, this is WIP. Also: 2-loop masses lack some log terms")
print("!!!")

print("Start finite-T stuff")
model3HDM.setInputParams(inputParams)

transitionFinder.traceFreeEnergyMinimum()





