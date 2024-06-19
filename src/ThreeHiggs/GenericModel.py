import numpy as np

from .EffectivePotential import EffectivePotential
from .DimensionalReduction import DimensionalReduction

class GenericModel():

    def __init__(self):
        
        self.inputParams = {}
        self.dimensionalReduction = DimensionalReduction()
        ## 3D potential:
        self.effectivePotential = EffectivePotential()

    def setInputParams(self, inputParams):
        self.inputParams = inputParams

    ## 1909.09234 and eq (38) for mass splittings
    def massSplittingsToMasses(self, mS1: float, delta12: float, delta1c: float, deltac: float) -> tuple[float, float, float]:
        mS2 = delta12 + mS1
        mSpm1 = delta1c + mS1
        mSpm2 = deltac + mSpm1
        return mS2, mSpm1, mSpm2
    
    ## Should acutally check vertices but this isn't a feature (yet) in DRalgo
    def checkPerturbativity(param : dict[str, float]) -> None:
        #print ("checkSingleCoupling called")
        if abs(param["lam11"]) > 4*np.pi or abs(param["lam12"]) > 4*np.pi or abs(param["lam12p"]) > 4*np.pi or abs(param["lam1Im"]) > 4*np.pi or abs(param["lam1Re"]) > 4*np.pi or abs(param["lam22"]) > 4*np.pi or abs(param["lam23"]) > 4*np.pi or abs(param["lam23p"]) > 4*np.pi or abs(param["lam2Im"]) > 4*np.pi or abs(param["lam2Re"]) > 4*np.pi or abs(param["lam31"]) > 4*np.pi or abs(param["lam31p"]) > 4*np.pi or abs(param["lam33"]) > 4*np.pi or abs(param["lam3Im"]) > 4*np.pi or abs(param["lam3Re"]) > 4*np.pi or abs(param["g1"]) > 4*np.pi or abs(param["g2"]) > 4*np.pi or abs(param["g3"]) > 4*np.pi:
            print ("Model is at risk of being non-pert")

    def calculateRenormalizedParameters(self, inputParams: dict[str, float]) -> dict[str, float]:
        """Take inputs from the BM file and convert them to parameters in the action.
        With tree-level matching the renormalization scale does not directly show up in the expressions, but
        needs to be specified for later loop calculations."""

        v = 246.22  ## "Higgs VEV". Consider using Fermi constant instead

        mS1 = inputParams["mS1"]

        if inputParams["bMassSplittingInput"]:
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
        res["yt3"] = np.sqrt(2.) * 172.76 / v  ## 172.76 is the top mass
        MW = 80.377 ## W boson mass
        MZ = 91.1876 ## Z boson mass
        ## Gauge couplings
        res["g1"] = 2.*np.sqrt(MZ**2 - MW**2) / v   # U(1)
        res["g2"] = 2.*MW / v                       # SU(2)
        res["g3"] = np.sqrt(0.1183 * 4.0 * np.pi)   # SU(3)

        ## Scalar potential parameters. Some of the RHS of these depend on LHS of others, so need to do in smart order
        MHsq = 125.00**2 ## Higgs mass squared
        res["mu3sq"] = MHsq / 2. 
        res["lam33"] = MHsq / (2.*v**2) 

        ### mu12^2 is real in our starting basis but not generally so after loop corrections
        mu12sq = 0.5 * (mSpm2**2 - mSpm1**2)
        res["mu12sqRe"] = mu12sq
        res["mu12sqIm"] = 0.0

        # Abs value of lambda2 (which is complex):
        lam2Abs = 1./v**2 * ( mu12sq*cosTheta + 0.25 * np.sqrt( (2.*mu12sq*cosTheta)**2 + (mS2**2 - mS1**2)**2 - (mSpm2**2 - mSpm1**2)**2 ) )

        res["lam2Re"] = lam2Abs * cosTheta
        res["lam2Im"] = lam2Abs * sinTheta

        #### Compute some helper parameters

        LambdaMinus = np.sqrt( mu12sq**2 + v**4*lam2Abs**2 - 2.*v**2*mu12sq*lam2Abs*cosTheta) ## Eq 13
        alpha = (-mu12sq + v**2*lam2Abs*cosTheta - LambdaMinus) / ( (v**2*lam2Abs*sinTheta) + 1e-100) ## Eq 12 - Adding 1e-100 to avoid 1/0

        mu2sq = v**2/2. * ghDM - v**2 / (alpha**2 + 1.) * (2.*alpha* sinTheta + (alpha**2 - 1.)*cosTheta) * lam2Abs - 0.5*(mS2**2 + mS1**2)

        res["mu2sq"] = mu2sq

        res["lam23"] = 1./v**2 * (2.*mu2sq + mSpm2**2 + mSpm1**2)
        res["lam23p"] = 1./v**2 * (mS2**2 + mS1**2 - mSpm2**2 - mSpm1**2)
    

        ## Set remaining params
        if inputParams["darkHierarchy"]:
            ## eq. (4)
            res["mu1sq"] = inputParams["darkHierarchy"]*res["mu2sq"]
            res["lam3Re"] = inputParams["darkHierarchy"]*res["lam2Re"]
            res["lam3Im"] = inputParams["darkHierarchy"]*res["lam2Im"]
            res["lam31"] = inputParams["darkHierarchy"]*res["lam23"]
            res["lam31p"] = inputParams["darkHierarchy"]*res["lam23p"]

        else:
            print("Input params only implemented in dark hieracy limit!")    
            raise NotImplementedError 


        res["RGScale"] = inputParams["RGScale"]

        return res
    