import numpy as np
import numpy.typing as npt
from typing import Tuple

from GenericModel import GenericModel
from BetaFunctions import BetaFunctions4D


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
    ##Test change increased step size in T from 1 to 2,5
    def traceFreeEnergyMinimum(self, TRange: npt.ArrayLike = np.arange(50., 200., 2.5)) -> Tuple[npt.ArrayLike, npt.ArrayLike]:

        TRange = np.asanyarray(TRange)

        renormalizedParams = self.model.calculateRenormalizedParameters(self.model.inputParams,  self.model.inputParams["RGScale"])

        """RG running. We want to do 4D -> 3D matching at a scale where logs are small; usually a T-dependent scale like 7T.
        To make this work nicely, integrate the beta functions here up to some high enough scale and store the resulting couplings
        in interpolated functions.
        """
        startScale = renormalizedParams["RGScale"]
        endScale = 7.3 * TRange[-1] ## largest T in our range is T[-1] 
        muRange = np.linspace( startScale, endScale, TRange.size*10 )

        ## TODO Is this safe if endScale is smaller than startScale?

        betas = BetaFunctions4D()
        ## Dict of interpolated functions (TODO change)
        interpolatedParams = betas.SolveBetaFunction(renormalizedParams, muRange)

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

            ## Final scale in 3D
            goalRGScale =  T

            matchingScale = 7.055 * T

            ## Fill in a new param dict that contains 4D params at the matching scale
            paramsForMatching = {}
            for key, interpolatedFunction in interpolatedParams.items():
                paramsForMatching[key] = interpolatedFunction(matchingScale)

            ## These need to be in the dict
            paramsForMatching["RGScale"] = matchingScale
            paramsForMatching["T"] = T

            ## Put T-dependent logs in the dict too. Not a particularly nice solution...
            Lb = 2. * np.log(matchingScale / T) - 2.*(np.log(4.*np.pi) - EulerGamma)
            paramsForMatching["Lb"] = Lb
            paramsForMatching["Lf"] = Lb + 4.*np.log(2.)

            params3D = self.model.dimensionalReduction.getEFTParams(paramsForMatching, goalRGScale)

            self.model.effectivePotential.setModelParameters(params3D)

            minimum, valueVeff = self.model.effectivePotential.findGlobalMinimum()

            minimizationResults.append( [T, valueVeff, *minimum] )
            #Test run
            print (T)

        minimizationResults = np.asanyarray(minimizationResults)
        print( minimizationResults )
        np.savetxt("results_test.txt", minimizationResults)


