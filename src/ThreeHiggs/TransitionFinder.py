import numpy as np

from .GenericModel import GenericModel
from .BetaFunctions import BetaFunctions4D

"""Class TransitionFinder -- This handles all logic for tracking the temperature dependence of a model,
identifying phase transitions, determining physical parameters of a transition etc. 
"""
class TransitionFinder:
    def __init__(self, model=None):

        if (model == None):
            model = GenericModel()

        self.model = model

    def traceFreeEnergyMinimum(self, TRange: np.ndarray = np.arange(50., 200., 1.)) -> tuple[np.ndarray, np.ndarray]:
        renormalizedParams = self.model.calculateRenormalizedParameters(self.model.inputParams)
        
        """RG running. We want to do 4D -> 3D matching at a scale where logs are small; usually a T-dependent scale like 7T.
        To make this work nicely, integrate the beta functions here up to some high enough scale and store the resulting couplings
        in interpolated functions.
        """
        startScale = renormalizedParams["RGScale"]
        endScale = 7.3 * TRange[-1] ## largest T in our range is T[-1] 
        muRange = np.linspace( startScale, endScale, TRange.size*10 )

        betas = BetaFunctions4D(muRange, renormalizedParams) ## TODO Are the beta function routines safe if endScale is smaller than startScale?
        
        EulerGamma = 0.5772156649
        EulerGammaPrime = 2.*(np.log(4.*np.pi) - EulerGamma)
        Lfconst = 4.*np.log(2.)
        
        minimizationResults = []
 
        counter = 0
        for T in TRange:
            print (f'Start of temp = {T} loop')
                       
            ## Final scale in 3D
            ## TODO ask Lauri if goalRGscale is ever different from just T
            goalRGScale =  T

            matchingScale = 4.0*np.pi*np.exp(-EulerGamma) * T
            
            paramsForMatching = betas.RunCoupling(matchingScale)
            
            from ThreeHiggs.GenericModel import bIsBounded, bIsPerturbative
            if not bIsBounded(paramsForMatching):
                print ("Model is not bounded from below, exiting")
                return
            
            ## These need to be in the dict
            paramsForMatching["RGScale"] = matchingScale
            paramsForMatching["T"] = T

            ## Put T-dependent logs in the dict too. Not a particularly nice solution...
            Lb = 2. * np.log(matchingScale / T) - EulerGammaPrime
            paramsForMatching["Lb"] = Lb
            paramsForMatching["Lf"] = Lb + Lfconst

            ##This has every coupling needed to compute the EP, computed at the matching scale (I think)
            params3D = self.model.dimensionalReduction.getEFTParams(paramsForMatching, goalRGScale)
            
            self.model.effectivePotential.setModelParameters(params3D)

            minimum, valueVeff = self.model.effectivePotential.findGlobalMinimum(T)
            
            minimizationResults.append( [T, valueVeff, *minimum] )

            if np.all(minimum < 1e-3):
                print (f"Symmetric phase found at temp {T}")
                if counter == 3:
                    break
                counter += 1

        return minimizationResults
