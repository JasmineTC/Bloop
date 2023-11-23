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


    def setupSoftToUltrasoftMatching(self, softToUltrasoftFile: str) -> None:
        
        self.matchToUltrasoft.createMatchingRelations(softToUltrasoftFile)


    def getSoftScaleParams(self, inputParams: dict[str, float]) -> dict[str, float]:

        return self.matchToSoft.getMatchedParams(inputParams)