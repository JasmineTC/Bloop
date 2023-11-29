import numpy as np
import numpy.typing as npt

from ParameterMatching import ParameterMatching


## Collects all needed matching relations to go from 4D params to 3D ultrasoft EFT
class DimensionalReduction():

    matchToSoft: ParameterMatching
    matchToUltrasoft: ParameterMatching

    def __init__(self):
        self.matchToSoft = ParameterMatching()
        self.matchToUltrasoft = ParameterMatching()



    def setupHardToSoftMatching(self, hardToSoftFile: str) -> None:
        
        self.matchToSoft.createMatchingRelations(hardToSoftFile)
        self.matchToSoft.matchingRelations = self.__remove3dSuffices(self.matchToSoft.matchingRelations)


    def setupSoftToUltrasoftMatching(self, softToUltrasoftFile: str) -> None:
        
        self.matchToUltrasoft.createMatchingRelations(softToUltrasoftFile)


    ## T should be in the dict
    def getSoftScaleParams(self, paramsForMatching: dict[str, float]) -> dict[str, float]:

        return self.matchToSoft.getMatchedParams(paramsForMatching)
    
    ## Modifies notation in matching relations so that the matched param dict does not use the "3d" suffix (comes from DRalgo by default)
    @staticmethod
    def __remove3dSuffices(matchingRelations: dict[str, float]) -> dict[str, float]:

        newDict = matchingRelations
        # DRalgo also gives gauge couplings as "g13d^2" etc, which is terrible => change to g1sq.
        # Crazy oneliner, creates a new dict where just the key names are different
        newDict = { key.replace("^2", "sq") : value for key, value in newDict.items() }

        """ # this didn't work for g13d^2 etc so commented out
        ## Remove "3d" suffix with even crazier oneliner (suffix meaning that it's removed only from end of the string)
        newDict = { key[:-len("3d")] if key.endswith("3d") else key : value for key, value in newDict.items() }
        """
        ## Remove '3d' substrings
        newDict = { key.replace('3d', '') : value for key, value in newDict.items() }

        return newDict