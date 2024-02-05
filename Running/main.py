import numpy as np
import pathlib

import ThreeHiggs

from ThreeHiggs import GenericModel
from ThreeHiggs import TransitionFinder

import Benchmarks.Benchmarks_3HDM

hardToSoftFile = ThreeHiggs.getResourcePath("Data/HardToSoft/softScaleParams_NLO.txt")
softToUltrasoftFile = ThreeHiggs.getResourcePath("Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt")


## Model object setup + load matching relations
model3HDM = GenericModel()
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)


#inputParams = Benchmarks.Benchmarks_3HDM.BM1
inputParams = Benchmarks.Benchmarks_3HDM.BM_SM_like

transitionFinder = TransitionFinder(model=model3HDM)

## Scanning loops would start here

print("!!!")
print("Currently not matching soft --> ultrasoft, this is WIP. Also: 2-loop masses lack some log terms")
print("!!!")

print("Start finite-T stuff")
model3HDM.setInputParams(inputParams)

transitionFinder.traceFreeEnergyMinimum()