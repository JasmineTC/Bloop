import numpy as np
import math

from .GenericModel import GenericModel
from .BetaFunctions import BetaFunctions4D

def threeDimFieldtoDimensionless(temp: list[float], field: list[float]) -> list[float]:
    return field/np.sqrt(temp)

"""Class TransitionFinder -- This handles all logic for tracking the temperature dependence of a model,
identifying phase transitions, determining physical parameters of a transition etc. 
"""
class TransitionFinder:
    def __init__(self, model=None):

        if (model == None):
            model = GenericModel()

        self.model = model

    def traceFreeEnergyMinimum(self, TRangeStart: float, 
                               TRangeEnd: float, 
                               TRangeStepSize: float) -> tuple[np.ndarray, np.ndarray]:
        renormalizedParams = self.model.calculateRenormalizedParameters(self.model.inputParams)
        
        """RG running. We want to do 4D -> 3D matching at a scale where logs are small; usually a T-dependent scale like 7T.
        To make this work nicely, integrate the beta functions here up to some high enough scale and store the resulting couplings
        in interpolated functions.
        """
        TRange = np.arange(TRangeStart, TRangeEnd, TRangeStepSize )
        startScale = renormalizedParams["RGScale"]
        endScale = 7.3 * TRange[-1] ## largest T in our range is T[-1] 
        muRange = np.linspace( startScale, endScale, TRange.size*10 )

        betas = BetaFunctions4D(muRange, renormalizedParams) ## TODO Are the beta function routines safe if endScale is smaller than startScale?
        
        EulerGamma = 0.5772156649
        EulerGammaPrime = 2.*(math.log(4.*np.pi) - EulerGamma)
        Lfconst = 4.*np.log(2.)
        
        minimizationResults = []
 
        counter = 0
        verbose = False
        for T in TRange:
            if verbose:
                print (f'Start of temp = {T} loop')
                       
            ## Final scale in 3D
            ## TODO ask Lauri if goalRGscale is ever different from just T
            goalRGScale =  T

            matchingScale = 4.0*np.pi*math.exp(-EulerGamma) * T
            
            paramsForMatching = betas.RunCoupling(matchingScale)
            
            from ThreeHiggs.GenericModel import bIsBounded, bIsPerturbative
            if not bIsBounded(paramsForMatching):
                minimizationResults.append( [1, 0, [1], False, False, False] )
                break
                
            
            ## These need to be in the dict
            paramsForMatching["RGScale"] = matchingScale
            paramsForMatching["T"] = T

            ## Put T-dependent logs in the dict too. Not a particularly nice solution...
            Lb = 2. * math.log(matchingScale / T) - EulerGammaPrime
            paramsForMatching["Lb"] = Lb
            paramsForMatching["Lf"] = Lb + Lfconst

            ##This has every coupling needed to compute the EP, computed at the matching scale (I think)
            params3D = self.model.dimensionalReduction.getEFTParams(paramsForMatching, goalRGScale)
            
            self.model.effectivePotential.setModelParameters(params3D)
            initialGuesses = [[0.1,0.1,0.1],
                              [-0.1,0.1,0.1],
                              [1e-3,1e-3,4],
                              [1e-3,1e-3,10],
                              [1e-3,1e-3,25], 
                              [5,5,1e-4],
                              [-5,5,1e-4],
                              [40,40,1e-4], 
                              [-40,40,1e-4],
                              [5,5,5],
                              [-5,5,5], 
                              [40,40,40],
                              [-40,40,40], 
                              [59,59,59], 
                              [-59,59,59]]
            minimumLocation, valueVeff = self.model.effectivePotential.findGlobalMinimum(initialGuesses)
            bReachedUltraSoftScale = self.model.effectivePotential.bReachedUltraSoftScale(minimumLocation, T)


            minimizationResults.append( [T, valueVeff, minimumLocation, bIsPerturbative(paramsForMatching), bReachedUltraSoftScale, 1] )

            if np.all(minimumLocation < 1e-3):
                if verbose:
                    print (f"Symmetric phase found at temp {T}")
                if counter == 3:
                    break
                counter += 1

        return self.convertResultsToDict(minimizationResults)
    
    def convertResultsToDict(self, minimizationResults):
        tempList = [float(result[0]) for result in minimizationResults]
        bReachedUltraSoftScaleList = [result[4] for result in minimizationResults]
        ##Gives the first index where the ultrasoft condition is True or -1 is none is found
        ultraSoftWarning: int = next((i for i, val in enumerate(bReachedUltraSoftScaleList) if val == True), -1) 
        TUltraSoft = tempList[ultraSoftWarning] if ultraSoftWarning >=0 else -1

        return {"T": tempList,
                "valueVeff": [result[1] for result in minimizationResults],
                "minimumLocation": np.transpose([result[2] for result in minimizationResults]).tolist(),
                "bIsPerturbative": all([result[3] for result in minimizationResults]),
                "UltraSoftTemp": TUltraSoft,
                "bBoundFromBelow": all([result[5] for result in minimizationResults]) }
