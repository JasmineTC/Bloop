import numpy as np
import numpy.typing as npt
from typing import Tuple

from EffectivePotential import EffectivePotential
from DimensionalReduction import DimensionalReduction

## Could be an ABC, but just implementing for 3HDM for now
class GenericModel():

    """ Control variables """
    # Set to true if the input dict contains 3HDM mass splittings (in GeV) instead of masses directly.
    bMassSplittingInput: bool = True

    dimensionalReduction: DimensionalReduction
    effectivePotential: EffectivePotential

    ## This is the "physical" input 
    inputParams: dict[str, float]

    """ Define some constants that are inputs too, but we won't scan over these. 
    """
    ## "Higgs VEV". Consider using Fermi constant instead
    v = 246.22
    ## Higgs mass
    MH = 125.00
    ## Gauge boson masses
    MW = 80.377
    MZ = 91.1876
    ## Top quark mass 
    Mt = 172.76
    ## SU(3) coupling, we neglect it's running. This is the value at Z pole
    g3 = np.sqrt(0.1183 * 4.0 * np.pi)


    def __init__(self):
        self.inputParams = {}
        self.dimensionalReduction = DimensionalReduction()
        ## 3D potential:
        self.effectivePotential = EffectivePotential()

    
    def setInputParams(self, inputParams):
        self.inputParams = inputParams

    ## Convert Venus' mass splitting (deltas) input to mass values. Returns mS2, mSpm1, mSpm2 (pm meaning plus/minus; the charged masses)
    @staticmethod
    def massSplittingsToMasses(mS1: float, delta12: float, delta1c: float, deltac: float) -> Tuple[float, float, float]:

        mS2 = delta12 + mS1
        mSpm1 = delta1c + mS1
        mSpm2 = deltac + mSpm1
        return mS2, mSpm1, mSpm2


    """This goes from whatever "physical" input to parameters in the actions.
    With tree-level matching the renormalization scale does not directly show up in the expressions, but
    needs to be specified for later loop calculations.

    Here the optional arguments are:
        bDarkDemocracyLimit=True if we set various phi1, phi2 specific couplings equal to each other as defined in the draft.
    """ 
    def calculateRenormalizedParameters(self, inputParams, RGScale: float, bDarkDemocracyLimit=True) -> dict[str, float]:

        v = self.v

        mS1 = inputParams["mS1"]

        if (self.bMassSplittingInput):

            mS2, mSpm1, mSpm2 = self.massSplittingsToMasses(mS1, inputParams["delta12"], inputParams["delta1c"], inputParams["deltac"])

        else:
            mS2 = inputParams["mS2"]
            ## Charged masses:
            mSpm1 = inputParams["mSpm1"]
            mSpm2 = inputParams["mSpm2"]

        ## This is the complex phase of lambda_2
        thetaCPV = inputParams["thetaCPV"]
        ghDM = inputParams["ghDM"]

        sinTheta = np.sin(thetaCPV)
        cosTheta = np.cos(thetaCPV)

        res = {}

        ## First copy action params that were given as direct inputs
        res["lam1Re"] = inputParams["lam1Re"]
        res["lam1Im"] = inputParams["lam1Im"]
        res["lam11"] = inputParams["lam11"]
        res["lam22"] = inputParams["lam22"]
        res["lam12"] = inputParams["lam12"]
        res["lam12p"] = inputParams["lam12p"]

        ## Yukawas (not squared! would it be better to store y^2?)
        res["yt3"] = np.sqrt(2.) * self.Mt / v 

        ## Gauge couplings
        res["g1"] = 2.*np.sqrt(self.MZ**2 - self.MW**2) / v     # Hypercharge
        res["g2"] = 2.*self.MW / v                              # SU(2)
        res["g3"] = self.g3                                     # SU(3)

        ## Scalar potential parameters. Some of the RHS of these depend on LHS of others, so need to do in smart order

        res["mu3sq"] = self.MH**2 / 2. 
        res["lam33"] = self.MH**2 / 2. / v**2

        ### mu12^2 is real in our starting basis but not generally so after loop corrections
        mu12sq = 0.5 * (mSpm2**2 - mSpm1**2)
        res["mu12sqRe"] = mu12sq
        res["mu12sqIm"] = 0.0

        # Abs value of lambda2 (which is complex):
        lam2Abs = 1./v**2 * ( mu12sq*cosTheta + 0.25 * np.sqrt( (2.*mu12sq*cosTheta)**2 + (mS2**2 - mS1**2)**2 - (mSpm2**2 - mSpm1**2)**2 ) )

        res["lam2Re"] = lam2Abs * cosTheta
        res["lam2Im"] = lam2Abs * sinTheta

        #### Compute some helper parameters

        ## eq (13)
        LambdaMinus = np.sqrt( mu12sq**2 + v**4*lam2Abs**2 - 2.*v**2*mu12sq*lam2Abs*cosTheta)
        ## eq (12)
        alpha = (-mu12sq + v**2*lam2Abs*cosTheta - LambdaMinus) / (v**2*lam2Abs*sinTheta)

        mu2sq = v**2/2. * ghDM - v**2 / (alpha**2 + 1.) * (2.*alpha* sinTheta + (alpha**2 - 1.)*cosTheta) * lam2Abs - 0.5*(mS2**2 + mS1**2)

        res["mu2sq"] = mu2sq

        res["lam23"] = 1./v**2 * (2.*mu2sq + mSpm2**2 + mSpm1**2)
        res["lam23p"] = 1./v**2 * (mS2**2 + mS1**2 - mSpm2**2 - mSpm1**2)
    

        ## Set remaining params
        if (bDarkDemocracyLimit):
            
            ## eq. (4)
            res["mu1sq"] = res["mu2sq"]
            res["lam3Re"] = res["lam2Re"]
            res["lam3Im"] = res["lam2Im"]
            res["lam31"] = res["lam23"]
            res["lam31p"] = res["lam23p"]

        else:
            print("Input params only implemented in ''dark democracy`` limit!")    
            raise NotImplementedError 


        res["RGScale"] = RGScale

        return res

        

    