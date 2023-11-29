import numpy as np
import numpy.typing as npt
from typing import Tuple

from GenericModel import GenericModel


"""Class TransitionFinder -- This handles all logic for tracking the temperature dependence of a model,
identifying phase transitions, determining physical parameters of a transition etc. 
"""
class TransitionFinder:

    model: GenericModel

    def __init__(self, model=None):

        if (model == None):
            model = GenericModel()

        self.model = model



    ## This is a way too big routine 
    def traceFreeEnergyMinimum(self, TRange: npt.ArrayLike = np.arange(50., 200., 1.)) -> Tuple[npt.ArrayLike, npt.ArrayLike]:

        renormalizedParams = self.model.calculateRenormalizedParameters(self.model.inputParams,  self.model.inputParams["RGScale"])

        ## Run the theory to EFT matching scale (usually ~7T). This is where we'd solve beta functions. But for now just match at the input scale
        paramsForMatching = renormalizedParams
        matchingScale = paramsForMatching["RGScale"]

        ## TODO would probs be good to move this part inside DimensionalReduction class

        ## Put T in the dict. I assume that the DR routines below work with T-array input
        paramsForMatching["T"] = TRange

        EulerGamma = 0.5772156649
        ## And T-dependent logs too. Not a particularly nice solution...
        Lb = 2. * np.log(matchingScale / TRange) - 2.*(np.log(4.*np.pi) - EulerGamma)
        paramsForMatching["Lb"] = Lb
        paramsForMatching["Lf"] = Lb + 4.*np.log(2.)

        params3D = self.model.dimensionalReduction.getSoftScaleParams(renormalizedParams)
        
        ## params3D should now be a dict of np.arrays of len(TRange)

        self.model.effectivePotential.setModelParameters(params3D)

        self.model.effectivePotential.findGlobalMinimum()



