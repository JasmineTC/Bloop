import pathlib

from ParameterMatching import ParameterMatching
from DimensionalReduction import DimensionalReduction
from Model import Model


## hack file path, not a good solution if making this into a proper package
pathToCurrentFile = pathlib.Path(__file__).parent.resolve()

hardToSoftFile = str(pathToCurrentFile) + "/Data/softScaleParams.txt"
# TODO
#softToUltrasoftFile = str(pathToCurrentFile) + ...


## Setup DR relations
dimensionalReduction = DimensionalReduction()
dimensionalReduction.setupHardToSoftMatching(hardToSoftFile)


inputParams = {
    'T' : 100,
    'Lb' : 1,
    'Lf' : 2,
    'g1' : 0.3,
    'g2' : 0.6,
    'g3' : 1.5,
    'yt3' : 0.98,
    ##
    'lam11' : 1,
    'lam22' : 2,
    'lam33' : 3,
    'lam12' : 12,
    'lam23' : 23,
    'lam31' : 31,
    'lam12p' : 12,
    'lam23p' : 23,
    'lam31p' : 31,
    'lam1Re' : 11,
    'lam1Im' : 11,
    'lam2Re' : 22,
    'lam2Im' : 22,
    'lam3Re' : 33,
    'lam3Im' : 33
}


matchedParams = dimensionalReduction.getSoftScaleParams(inputParams)

print(matchedParams)


model3HDM = Model()

RGScale = 
model3HDM.calculateRenormalizedParameters()

