import pathlib
import numpy as np
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

#inputParams = Benchmarks.Benchmarks_3HDM.BM1
#inputParams = Benchmarks.Benchmarks_3HDM.BM_SM_like

transitionFinder = TransitionFinder(model=model3HDM)

print("!!!")
print("Currently not matching soft --> ultrasoft, this is WIP. Also: 2-loop masses lack some log terms")
print("!!!")

print("Start finite-T stuff")


BM_list = [Benchmarks.Benchmarks_3HDM.BM1, Benchmarks.Benchmarks_3HDM.BM2, Benchmarks.Benchmarks_3HDM.BM3,
           Benchmarks.Benchmarks_3HDM.BM4, Benchmarks.Benchmarks_3HDM.BM5, Benchmarks.Benchmarks_3HDM.BM6,
           Benchmarks.Benchmarks_3HDM.BM7, Benchmarks.Benchmarks_3HDM.BM8, Benchmarks.Benchmarks_3HDM.BM9]
for i, inputParams in enumerate(BM_list):

    ## Scanning loops would start here
    model3HDM.setInputParams(inputParams)

    minimizationResults = transitionFinder.traceFreeEnergyMinimum()
    
    fileName =  str(pathToCurrentFile) + "/Data/Results/BM" + str(i+1) + ".txt"
    
    np.savetxt(fileName, minimizationResults)
    
# minimizationResults = transitionFinder.traceFreeEnergyMinimum()

# np.savetxt("results.txt", minimizationResults)