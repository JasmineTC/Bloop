import pathlib
import numpy as np
from GenericModel import GenericModel
from TransitionFinder import TransitionFinder

import Benchmarks.Benchmarks_3HDM
from Userinput import Userinput

userinput = Userinput()
args = userinput.parse()

## hack file path, not a good solution if making this into a proper package

pathToCurrentFile = pathlib.Path(__file__).parent.resolve()

hardToSoftFile = str(pathToCurrentFile) + "/Data/HardToSoft/softScaleParams_NLO.txt"
softToUltrasoftFile = str(pathToCurrentFile) + "/Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt"

##TODO set up loading the file from command line or default arg
##TODO set up mutliprocressing
# ## Model object setup + load matching relations
model3HDM = GenericModel(loopOrder = args.loopOrder)
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)

print("!!! \n Currently not matching soft --> ultrasoft, this is WIP. Also: 2-loop masses lack some log terms  \n !!!")

print("Start finite-T stuff")

# #inputParams = Benchmarks.Benchmarks_3HDM.BM1
# #inputParams = Benchmarks.Benchmarks_3HDM.BM_SM_like

# transitionFinder = TransitionFinder(model=model3HDM)
# model3HDM.setInputParams(inputParams)
# # minimizationResults = transitionFinder.traceFreeEnergyMinimum()
# # np.savetxt("results.txt", minimizationResults)


# pull in the whole list from edited data file
# BM_list = Benchmarks.Benchmarks_3HDM.BM_list

# # we will either supply the string "all" as cmd arg or some integers
# if sys.argv[1] == "all":
#     bmIndex_list = list(range(0, len(BM_list)))
#     BM_list_desired = BM_list
# else:
#     # casting arguments to actually be integers                
#     bmIndex_list = [int(bmIndex) for bmIndex in sys.argv[1:]]
#     # filter the list into a desired list by checking indexes against integer arguments
#     BM_list_desired = [bm for ind, bm in enumerate(BM_list) if ind in bmIndex_list]

# # just changed this to now operate on BM_list_desired
# for i, inputParams in enumerate(BM_list_desired):
#     print ("Now running bench mark %s" %bmIndex_list[i])
#     ## Scanning loops would start here
#     model3HDM.setInputParams(inputParams)

#     minimizationResults = transitionFinder.traceFreeEnergyMinimum()
   
#     fileName =  str(pathToCurrentFile) + "/Data/Results2/BM" + str(bmIndex_list[i]) + ".txt"
    
#     np.savetxt(fileName, minimizationResults)