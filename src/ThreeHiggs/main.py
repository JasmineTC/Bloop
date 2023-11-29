import numpy as np
import pathlib

from ParameterMatching import ParameterMatching
from DimensionalReduction import DimensionalReduction
from GenericModel import GenericModel
from TransitionFinder import TransitionFinder

## hack file path, not a good solution if making this into a proper package
pathToCurrentFile = pathlib.Path(__file__).parent.resolve()

hardToSoftFile = str(pathToCurrentFile) + "/Data/softScaleParams.txt"
# TODO
#softToUltrasoftFile = str(pathToCurrentFile) + ...


## Model object setup + load matching relations
model3HDM = GenericModel()
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile)

## Input at Z pole
inputScale = model3HDM.MZ

## This is Venus' first benchmark point in table 1 of draft
inputParams = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : inputScale,

    ## "Physical" input in the scalar sector
    "mS1" : 67,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 4.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}

transitionFinder = TransitionFinder(model=model3HDM)

## Scanning loops would start here

print("Start finite-T stuff")
model3HDM.setInputParams(inputParams)

transitionFinder.traceFreeEnergyMinimum()





