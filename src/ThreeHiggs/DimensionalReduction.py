import numpy as np
import numpy.typing as npt
from typing import Callable
from copy import deepcopy

from .ParameterMatching import ParameterMatching


## Collects all needed matching relations to go from 4D params to 3D ultrasoft EFT
class DimensionalReduction():

    def __init__(self):
        self.matchToSoft = ParameterMatching()
        self.matchToUltrasoft = ParameterMatching()

        ## These default to False and are only toggled on if RGE files are specified when setupping matching relations
        self.bDoSoftScaleRGE = False
        self.bDoUltrasoftScaleRGE = False


    def setupHardToSoftMatching(self, hardToSoftFile: str, softScaleRGEFile: str = None) -> None:
        """softScaleRGEFile specifies where RGEs are loaded from. If left to None, 
        will not perform RG running at soft scale."""

        self.matchToSoft.createMatchingRelations(hardToSoftFile)
        #self.matchToSoft.matchingRelations = self.__remove3dSuffices(self.matchToSoft.matchingRelations)

        print("Setup Hard -> Soft matching relations.")
        print("-- Inputs:")
        print(self.matchToSoft.parameterNames)
        print("-- Outputs:")
        print( list(self.matchToSoft.matchingRelations.keys()) )
        print("")

        if (softScaleRGEFile):
            self.bDoSoftScaleRGE = True
            ## Using ParameterMatching for these because we want to remove 3d suffices etc
            self.softScaleRGE = ParameterMatching()
            self.softScaleRGE.createMatchingRelations(softScaleRGEFile)
            #self.softScaleRGE.matchingRelations = self.__remove3dSuffices(self.softScaleRGE.matchingRelations)
            print("Soft scale RGE")
            print("-- Inputs:")
            print(self.softScaleRGE.parameterNames)
            print("-- Outputs:")
            print( list(self.softScaleRGE.matchingRelations.keys()) )


    def setupSoftToUltrasoftMatching(self, softToUltrasoftFile: str) -> None:
        
        self.matchToUltrasoft.createMatchingRelations(softToUltrasoftFile)
        self.matchToUltrasoft.matchingRelations = self.__remove3dSuffices(self.matchToUltrasoft.matchingRelations, bRemoveSuffixUS=True)

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
        #ultrasoftScaleParams = self.matchToUltrasoft.getMatchedParams(softScaleParams)

        ultrasoftScaleParams = softScaleParams


        ## HACK this is the RG scale name in Veff
        ultrasoftScaleParams["mu3US"] = goalRGScale

        return ultrasoftScaleParams



    ## NB: T should be in the input dict
    def getSoftScaleParams(self, paramsForMatching: dict[str, float], goalRGScale: float) -> dict[str, float]:
        """Match hard scale --> soft scale theory
        """
        outParams = self.matchToSoft.getMatchedParams(paramsForMatching)

        ## RG scale needs to be in the parameter dict
        outParams["RGScale"] = paramsForMatching["RGScale"]

        if (self.bDoSoftScaleRGE):
            outParams = self.solveSoftScaleRGE(outParams, goalRGScale)
        else:
            ## No explicit running, but still need to have the RG scale in the out dict
            outParams["RGscale"] = goalRGScale

        return outParams
    

    def getUltrasoftScaleParams(self, softScaleParams: dict[str, float], goalRGScale: float) -> dict[str, float]:
        """Match soft scale --> ultrasoft scale theory
        """
        outParams = self.matchToUltrasoft.getMatchedParams(softScaleParams)
        ## TODO RG running 

        return outParams

    
    def solveSoftScaleRGE(self, initialParams: dict[str, float], goalScale: float) -> dict[str, float]:
        """Evolves parameters in soft scale EFT to goalScale. Starting scale is read from initialParams["RGScale"].
        This modifies the input dict in place!
        """
        #startScale = initialParams["RGScale"]
        initialParams["goalScale"] = goalScale
        initialParams["startScale"] = initialParams["RGScale"]

        evolvedParams = self.softScaleRGE.getMatchedParams(initialParams)

        ## The above gives masses only (usually). So merge it to the initial dict to get all params
        initialParams.update(evolvedParams)
        initialParams["RGScale"] = goalScale
        return initialParams


    
    @staticmethod
    def __remove3dSuffices(matchingRelations: dict[str, any], bRemoveSuffixUS=False) -> dict[str, any]:
        """Modifies notation in matching relations so that the matched param dict does not use the "3d" suffix (comes from DRalgo by default)
        """

        newDict = deepcopy(matchingRelations)

        # DRalgo also gives gauge couplings as "g13d^2" etc, which is terrible => change to g1sq.
        # Crazy oneliner, creates a new dict where just the key names are different:
        newDict = { key.replace("^2", "sq") : value for key, value in newDict.items() }

        ## Remove "3d" suffix with even crazier oneliner (suffix meaning that it's removed only from end of the string)
        newDict = { key[:-len("3d")] if key.endswith("3d") else key : value for key, value in newDict.items() }
        ## Gauge couplings are originally of form g3d^2 so account for that too 
        newDict = { key.replace("3dsq", "sq") if key.endswith("3dsq") else key : value for key, value in newDict.items() }
    
        if (bRemoveSuffixUS):
            """ For ultrasoft theory DRalgo appends "US" => remove that too. Gauge couplings again need special treatment"""
            newDict = { key[:-len("US")] if key.endswith("US") else key : value for key, value in newDict.items() }
            newDict = { key.replace("USsq", "sq") if key.endswith("USsq") else key : value for key, value in newDict.items() }

        return newDict