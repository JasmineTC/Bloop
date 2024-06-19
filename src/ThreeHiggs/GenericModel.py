import math as m

from .EffectivePotential import EffectivePotential
from .DimensionalReduction import DimensionalReduction

class GenericModel():
    def __init__(self):
        self.inputParams = {}
        self.dimensionalReduction = DimensionalReduction()
        self.effectivePotential = EffectivePotential()  ## 3D potential:

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
        if abs(param["lam11"]) > 4*m.pi or abs(param["lam12"]) > 4*m.pi or abs(param["lam12p"]) > 4*m.pi or abs(param["lam1Im"]) > 4*m.pi or abs(param["lam1Re"]) > 4*m.pi or abs(param["lam22"]) > 4*m.pi or abs(param["lam23"]) > 4*m.pi or abs(param["lam23p"]) > 4*m.pi or abs(param["lam2Im"]) > 4*m.pi or abs(param["lam2Re"]) > 4*m.pi or abs(param["lam31"]) > 4*m.pi or abs(param["lam31p"]) > 4*m.pi or abs(param["lam33"]) > 4*m.pi or abs(param["lam3Im"]) > 4*m.pi or abs(param["lam3Re"]) > 4*m.pi or abs(param["g1"]) > 4*m.pi or abs(param["g2"]) > 4*m.pi or abs(param["g3"]) > 4*m.pi:
            print ("Model is at risk of being non-pert")

    def calculateRenormalizedParameters(self, inputParams: dict[str, float]) -> dict[str, float]:
        """Take inputs from the BM file and convert them to parameters in the action.
        With tree-level matching the renormalization scale does not directly show up in the expressions, but
        needs to be specified for later loop calculations."""
        res = {}
        ## --- SM fermions and gauge bosons ---
        v = 246.22  ## "Higgs VEV". Consider using Fermi constant instead
        res["yt3"] = m.sqrt(2.) * 172.76 / v  ## 172.76 is the top mass (not squared! would it be better to store y^2?)
        MW = 80.377 ## W boson mass
        MZ = 91.1876 ## Z boson mass
        ## Gauge couplings
        res |= {"g1": 2.*m.sqrt(MZ**2 - MW**2)/ v,   # U(1)
                "g2": 2.*MW/ v,                      # SU(2)
                "g3": m.sqrt(0.1183 * 4.0 * m.pi)}   # SU(3)
        
        ## --- BSM scalars ---
        if inputParams["bpreCalculated"]:  
            res |= inputParams["couplings"]
            
        else:
            mS1 = inputParams["mS1"]
            
            if inputParams["bMassSplittingInput"]:
                mS2, mSpm1, mSpm2 = self.massSplittingsToMasses(mS1, inputParams["delta12"], inputParams["delta1c"], inputParams["deltac"])
    
            else:
                mS2 = inputParams["mS2"]
                mSpm1 = inputParams["mSpm1"]
                mSpm2 = inputParams["mSpm2"]
    
            res |= {"lam1Re": inputParams["lam1Re"],
                    "lam1Im": inputParams["lam1Im"],
                    "lam11": inputParams["lam11"],
                    "lam22": inputParams["lam22"],
                    "lam12": inputParams["lam12"],
                    "lam12p": inputParams["lam12p"]}

            ## Scalar potential parameters. Some of the RHS of these depend on LHS of others, so need to do in smart order
            MHsq = 125.00**2 ## Higgs mass squared
            res |= {"mu3sq": MHsq / 2.,
                    "lam33": MHsq / (2.*v**2) }
    
            mu12sq = 0.5 * (mSpm2**2 - mSpm1**2)
            res |= {"mu12sqRe": mu12sq,
                    "mu12sqIm": 0.} ### mu12^2 is real at tree level (basis choice) but generally complex at loop level
            
            sinTheta, cosTheta = m.sin(inputParams["thetaCPV"]), m.cos(inputParams["thetaCPV"])
            lam2Abs = 1./v**2 * ( mu12sq*cosTheta + 0.25 * m.sqrt( (2.*mu12sq*cosTheta)**2 + (mS2**2 - mS1**2)**2 - (mSpm2**2 - mSpm1**2)**2 ) )
            res |= {"lam2Re": lam2Abs * cosTheta,
                    "lam2Im": lam2Abs * sinTheta}
            
    
            #### Compute some helper parameters
            LambdaMinus = m.sqrt( mu12sq**2 + v**4*lam2Abs**2 - 2.*v**2*mu12sq*lam2Abs*cosTheta) ## Eq 13
            alpha = (-mu12sq + v**2*lam2Abs*cosTheta - LambdaMinus) / ( (v**2*lam2Abs*sinTheta) + 1e-100) ## Eq 12 - Adding 1e-100 to avoid 1/0
    
            mu2sq = v**2/2. * inputParams["ghDM"] - v**2 / (alpha**2 + 1.) * (2.*alpha* sinTheta + (alpha**2 - 1.)*cosTheta) * lam2Abs - 0.5*(mS2**2 + mS1**2)
            res |= {"mu2sq": mu2sq,
                    "lam23": 1./v**2 * (2.*mu2sq + mSpm2**2 + mSpm1**2),
                    "lam23p": 1./v**2 * (mS2**2 + mS1**2 - mSpm2**2 - mSpm1**2)}
        
            if inputParams["darkHierarchy"]: ## eq. (4) but generalised to Hieracy (inputParams["darkHierarchy"] = 1 gives democracy)
                res |= {"mu1sq": inputParams["darkHierarchy"]*res["mu2sq"],
                        "lam3Re": inputParams["darkHierarchy"]*res["lam2Re"],
                        "lam3Im": inputParams["darkHierarchy"]*res["lam2Im"],
                        "lam31": inputParams["darkHierarchy"]*res["lam23"],
                        "lam31p": inputParams["darkHierarchy"]*res["lam23p"]}
    
            else:
                print("Input params only implemented in dark hieracy limit!")    
                raise NotImplementedError 
    
            res["RGScale"] = inputParams["RGScale"] ## This needs to be last due to bug in beta function!!!

        return res
    