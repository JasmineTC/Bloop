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


        """ 
        Now for the temperature loop. I see two options:
            1. Give the TRange array directly to EFT routines and the Veff => DR results dict of arrays of len(TRange).
            This way numpy vectorization would be automatic. HOWEVER, minimization with scipy.optimize.minimize is not 
            directly possible because Veff(x) would be an array of len(TRange), and scipy requires a scalar value.
            So this does not work currently.

            2. Do an ordinary loop over values in TRange. This is necessarily slow in Python but works. Good thing here is that
            we can break the T-loop at any point, eg. when we find a transition of interest.

        Going with option 2. 
        """

        EulerGamma = 0.5772156649

        ## This will contain minimization results in form: 
        ## [ [T, Veff(min), field1, field2, ...], ... ]
        minimizationResults = []
        for T in TRange:

            ## T needs to be in the dict
            paramsForMatching["T"] = T

            ## Put T-dependent logs in the dict too. Not a particularly nice solution...
            Lb = 2. * np.log(matchingScale / T) - 2.*(np.log(4.*np.pi) - EulerGamma)
            paramsForMatching["Lb"] = Lb
            paramsForMatching["Lf"] = Lb + 4.*np.log(2.)

            params3D = self.model.dimensionalReduction.getSoftScaleParams(renormalizedParams)

            self.model.effectivePotential.setModelParameters(params3D)

            minimum, valueVeff = self.model.effectivePotential.findGlobalMinimum()

            minimizationResults.append( [T, valueVeff, *minimum] )

        minimizationResults = np.array(minimizationResults)
        print( minimizationResults )
        np.savetxt("results_test.txt", minimizationResults)


