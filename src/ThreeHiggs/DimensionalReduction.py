## Collects all needed matching relations to go from 4D params to 3D ultrasoft EFT
class DimensionalReduction():

    def __init__(self, hardToSoft, softScaleRGE, softToUltraSoft):
        self.matchToSoft = hardToSoft
        self.softScaleRGE = softScaleRGE
        self.matchToUltrasoft = softToUltraSoft

        verbose = False ## WIP set this from user arg
        if verbose:
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
        ultrasoftScaleParams = self.matchToUltrasoft(softScaleParams, bReturnDict = True)

        ## HACK this is the RG scale name in Veff
        ultrasoftScaleParams["mu3US"] = goalRGScale

        return ultrasoftScaleParams



    ## NB: T should be in the input dict
    def getSoftScaleParams(self, paramsForMatching: dict[str, float], goalScale: float) -> dict[str, float]:
        """Match hard scale --> soft scale theory
        """
        outParams = self.matchToSoft(paramsForMatching, bReturnDict = True)

        ## RG scale needs to be in the parameter dict
        outParams["RGScale"] = paramsForMatching["RGScale"]
        outParams["goalScale"] = goalScale
        outParams["startScale"] = outParams["RGScale"]

        ## The above gives masses only (usually). So merge it to the initial dict to get all params
        outParams |= self.softScaleRGE(outParams, bReturnDict = True)
        outParams["RGScale"] = goalScale

        return outParams

    # Unused. Remove me.
    def getUltrasoftScaleParams(self, softScaleParams: dict[str, float], goalRGScale: float) -> dict[str, float]:
        return self.matchToUltrasoft(softScaleParams)

