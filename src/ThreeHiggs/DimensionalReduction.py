import numpy as np
import numpy.typing as npt
from typing import Callable
from copy import deepcopy

from .ParameterMatching import ParameterMatching


## Collects all needed matching relations to go from 4D params to 3D ultrasoft EFT
class DimensionalReduction():

    def __init__(self, hardToSoft, softScaleRGE, softToUltraSoft):
        self.matchToSoft = hardToSoft
        self.softScaleRGE = softScaleRGE
        self.matchToUltrasoft = softToUltraSoft
        self.matchToUltrasoft.matchingRelations = self.__remove3dSuffices(self.matchToUltrasoft.matchingRelations, bRemoveSuffixUS=True)

        print("Setup Hard -> Soft matching relations.")
        print("-- Inputs:")
        print(self.matchToSoft.parameterNames)
        print("-- Outputs:")
        print( list(self.matchToSoft.matchingRelations.keys()) )
        print("")

        print("Setup Soft -> Ultrasoft matching relations.")
        print("-- Inputs:")
        print(self.matchToUltrasoft.parameterNames)
        print("-- Outputs:")
        print( list(self.matchToUltrasoft.matchingRelations.keys()) )
        print("")

    def getEFTParams(self, paramsForMatching: dict[str, float], goalRGScale: float) -> dict[str, float]:
        """This goes from input hard scale parameters to whatever the final EFT is.
        """
        softScaleParams = self.getSoftScaleParams(paramsForMatching, goalRGScale)
        ultrasoftScaleParams = self.matchToUltrasoft(softScaleParams)

        ## HACK this is the RG scale name in Veff
        ultrasoftScaleParams["mu3US"] = goalRGScale

        return ultrasoftScaleParams



    ## NB: T should be in the input dict
    def getSoftScaleParams(self, paramsForMatching: dict[str, float], goalScale: float) -> dict[str, float]:
        """Match hard scale --> soft scale theory
        """
        outParams = self.matchToSoft(paramsForMatching)

        ## RG scale needs to be in the parameter dict
        outParams["RGScale"] = paramsForMatching["RGScale"]
        outParams["goalScale"] = goalScale
        outParams["startScale"] = outParams["RGScale"]

        ## The above gives masses only (usually). So merge it to the initial dict to get all params
        outParams |= self.softScaleRGE(outParams)
        outParams["RGScale"] = goalScale

        return outParams

    # Unused. Remove me.
    def getUltrasoftScaleParams(self, softScaleParams: dict[str, float], goalRGScale: float) -> dict[str, float]:
        return self.matchToUltrasoft(softScaleParams)

    @staticmethod
    def __remove3dSuffices(matchingRelations: dict[str, any], bRemoveSuffixUS=False) -> dict[str, any]:
        """Modifies notation in matching relations so that the matched param dict does not use the "3d" suffix (comes from DRalgo by default)
        """

        newDict = matchingRelations

        # DRalgo also gives gauge couplings as "g13d^2" etc, which is terrible => change to g1sq.
        # Crazy oneliner, creates a new dict where just the key names are different:
        newDict = { key.replace("^2", "sq") : value for key, value in newDict.items() }

        if (bRemoveSuffixUS):
            """ For ultrasoft theory DRalgo appends "US" => remove that too. Gauge couplings again need special treatment."""
            newDict = { key[:-len("US")] if key.endswith("US") else key : value for key, value in newDict.items() }
            newDict = { key.replace("USsq", "sq") if key.endswith("USsq") else key : value for key, value in newDict.items() }

        ## Remove "3d" suffix with even crazier oneliner (suffix meaning that it's removed only from end of the string)
        newDict = { key[:-len("3d")] if key.endswith("3d") else key : value for key, value in newDict.items() }
        ## Gauge couplings are originally of form g3d^2 so account for that too 
        newDict = { key.replace("3dsq", "sq") if key.endswith("3dsq") else key : value for key, value in newDict.items() }
    
        return newDict
