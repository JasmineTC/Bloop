import numpy as np
from math import sqrt, sin, cos, pi, log, exp

def bIsBounded(param : dict[str, float]) -> bool:
    ## Taking equations 26-31 from the draft that ensure the potential is bounded from below.
    lamx = param["lam12"] + min(0, param["lam12p"] - 2*sqrt(param["lam1Re"]**2 + param["lam1Im"]**2) )
    lamy = param["lam31"] + min(0, param["lam31p"] - 2*sqrt(param["lam3Re"]**2 + param["lam3Im"]**2) )
    lamz = param["lam23"] + min(0, param["lam23p"] - 2*sqrt(param["lam2Re"]**2 + param["lam2Im"]**2) )

    return param["lam11"] > 0 and \
           param["lam22"] > 0 and \
           param["lam33"] > 0 and \
           lamx > -2*sqrt(param["lam11"]*param["lam22"]) and \
           lamy > -2*sqrt(param["lam11"]*param["lam33"]) and \
           lamz > -2*sqrt(param["lam22"]*param["lam33"]) and \
           (sqrt(param["lam33"])*lamx + sqrt(param["lam11"])*lamz + sqrt(param["lam22"])*lamy >= 0 or \
           param["lam33"]*lamx**2 + param["lam11"]*lamz**2 + param["lam22"]*lamy**2 -param["lam11"]*param["lam22"]*param["lam33"] - 2*lamx*lamy*lamz < 0)

def bIsPerturbative(param : dict[str, float]) -> bool:
    ## Should actually check vertices but not a feature in DRalgo at time of writting
    return abs(param["lam11"]) < 4*pi and \
       abs(param["lam12"]) < 4*pi and \
       abs(param["lam12p"]) < 4*pi and \
       abs(param["lam1Im"]) < 4*pi and \
       abs(param["lam1Re"]) < 4*pi and \
       abs(param["lam22"]) < 4*pi and \
       abs(param["lam23"]) < 4*pi and \
       abs(param["lam23p"]) < 4*pi and \
       abs(param["lam2Im"]) < 4*pi and \
       abs(param["lam2Re"]) < 4*pi and \
       abs(param["lam31"]) < 4*pi and \
       abs(param["lam31p"]) < 4*pi and \
       abs(param["lam33"]) < 4*pi and \
       abs(param["lam3Im"]) < 4*pi and \
       abs(param["lam3Re"]) < 4*pi and \
       abs(param["g1"]) < 4*pi and \
       abs(param["g2"]) < 4*pi and \
       abs(param["g3"]) < 4*pi

def massSplittingsToMasses(mS1: float, delta12: float, delta1c: float, deltac: float) -> tuple[float, float, float]:
    ## 1909.09234 and eq (38) for mass splittings
    mS2 = delta12 + mS1
    mSpm1 = delta1c + mS1
    mSpm2 = deltac + mSpm1
    return mS2, mSpm1, mSpm2

def get4DLagranianParams(inputParams: dict[str, float]) -> dict[str, float]:
    """Take inputs from the BM file and convert them to parameters in the action.
    With tree-level matching the renormalization scale does not directly show up in the expressions, but
    needs to be specified for later loop calculations."""
    res = {}
    ## --- SM fermions and gauge bosons ---
    v = 246.22  ## "Higgs VEV". Consider using Fermi constant instead
    res["yt3"] = sqrt(2.) * 172.76 / v  ## 172.76 is the top mass (not squared! would it be better to store y^2?)
    MW = 80.377 ## W boson mass
    MZ = 91.1876 ## Z boson mass
    
    res |= {"g1": 2.*sqrt(MZ**2 - MW**2)/ v,  # U(1)
            "g2": 2.*MW/ v,                   # SU(2)
            "g3": sqrt(0.1183 * 4.0 * pi)}    # SU(3)
    
    ## --- BSM scalars ---
    if inputParams["bPreCalculated"]:  
        res |= inputParams["couplingValues"]
        res["RGScale"] = inputParams["RGScale"]
        
    else:
        mS1 = inputParams["mS1"]
        
        if inputParams["bMassSplittingInput"]:
            mS2, mSpm1, mSpm2 = massSplittingsToMasses(mS1, inputParams["delta12"], inputParams["delta1c"], inputParams["deltac"])

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
        
        sinTheta, cosTheta = sin(inputParams["thetaCPV"]), cos(inputParams["thetaCPV"])
        lam2Abs = 1./v**2 * ( mu12sq*cosTheta + 0.25 * sqrt( (2.*mu12sq*cosTheta)**2 + (mS2**2 - mS1**2)**2 - (mSpm2**2 - mSpm1**2)**2 ) )
        res |= {"lam2Re": lam2Abs * cosTheta,
                "lam2Im": lam2Abs * sinTheta}
        

        #### Compute some helper parameters
        LambdaMinus = sqrt( mu12sq**2 + v**4*lam2Abs**2 - 2.*v**2*mu12sq*lam2Abs*cosTheta) ## Eq 13
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

def threeDimFieldtoDimensionless(temp: list[float], field: list[float]) -> list[float]:
    return field/np.sqrt(temp)

def traceFreeEnergyMinimum(effectivePotential,
                           dimensionalReduction,
                           benchmark,
                           TRangeStart: float, 
                           TRangeEnd: float, 
                           TRangeStepSize: float,
                           verbose = False) -> dict[str: ]:
    """RG running. We want to do 4D -> 3D matching at a scale where logs are small; usually a T-dependent scale ~7T.
    To make this work nicely we integrate the beta functions here up to the largest temp used 
    then interpolate over the beta function."""
    TRange = np.arange(TRangeStart, TRangeEnd, TRangeStepSize)
    
    LagranianParams4D = get4DLagranianParams(benchmark)
    startScale = LagranianParams4D["RGScale"]
    endScale = 7.3 * TRange[-1] 
    muRange = np.linspace(startScale, endScale, TRange.size*10)
    
    from .BetaFunctions import BetaFunctions4D
    betasFunctions = BetaFunctions4D(muRange, LagranianParams4D) 
    
    EulerGammaPrime = 2.*(log(4.*pi) - np.euler_gamma)
    Lfconst = 4.*log(2.)
    
    minimizationResults = {"T": [],
                           "valueVeffReal": [],
                           "valueVeffImag": [],
                           "complex": False,
                           "minimumLocation": [], 
                           "bIsPerturbative": True, 
                           "UltraSoftTemp": None, 
                           "failureReason": None}

    counter = 0
    for T in TRange:
        T = float(T) ## To make compatible with JSON
        if verbose:
            print (f'Start of temp = {T} loop')
        minimizationResults["T"].append(T)
        goalRGScale =  T ## Final scale in 3D -check if goalRGscale is ever different from just T

        matchingScale = 4.*pi*exp(-np.euler_gamma) * T ## Scale that minimises T dependent logs

        paramsForMatching = betasFunctions.runCoupling(matchingScale)
        
        if not bIsBounded(paramsForMatching):
            minimizationResults["failureReason"] = "Unbounded"
            break
        
        ## T dependent Variables needed to compute the 2 loop EP but aren't Lagranian params 
        paramsForMatching["RGScale"] = matchingScale
        paramsForMatching["T"] = T
        Lb = 2. * log(matchingScale / T) - EulerGammaPrime
        paramsForMatching["Lb"] = Lb
        paramsForMatching["Lf"] = Lb + Lfconst

        ##Take the 4D params (at the matching scale) and match them to the 3D params
        params3D = dimensionalReduction.getEFTParams(paramsForMatching, goalRGScale)
        effectivePotential.setModelParameters(params3D)
        
        initialGuesses = [[0.1,0.1,0.1], ## TODO This should go in a config file or something
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
        minimumLocation, minimumValueReal, minimumValueImag, status = effectivePotential.findGlobalMinimum(initialGuesses)

        if status == "NaN": 
            minimizationResults["failureReason"] = "MinimisationFailed"
            break
        
        minimizationResults["valueVeffReal"].append(minimumValueReal)
        minimizationResults["valueVeffImag"].append(minimumValueImag)
        minimizationResults["minimumLocation"].append(minimumLocation)
        if not minimizationResults["complex"]:
            minimizationResults["complex"] = status == "complex"
        
        if not minimizationResults["UltraSoftTemp"]: ## If the ultra soft temp has not yet been set
            if effectivePotential.bReachedUltraSoftScale(minimumLocation, ## Check if ultra soft scale reached
                                                         T, 
                                                         verbose = verbose): 
                minimizationResults["UltraSoftTemp"] = T ## If reached then set that as the ultra soft temp
                                                         ##- this will stop the first if statement from passing
        
        if minimizationResults["bIsPerturbative"]: ##If the potential was perturbative check if it still is
            minimizationResults["bIsPerturbative"] = bIsPerturbative(paramsForMatching) ## If non-pert then value set to false and won't be updated

        if np.all(minimumLocation < 1e-3):
            if verbose:
                print (f"Symmetric phase found at temp {T}")
            if counter == 3:
                break
            counter += 1

    minimizationResults["minimumLocation"] = np.transpose(minimizationResults["minimumLocation"]).tolist()
    return minimizationResults



from unittest import TestCase
class TransitionFinderUnitTests(TestCase):
    def test_bIsPerturbativeTrue(self):
        reference = True
        source = {"lam11": 0,
                  "lam12": 0,
                  "lam12p": 0,
                  "lam1Im": 0,
                  "lam1Re": 0,
                  "lam22": 0,
                  "lam23": 0,
                  "lam23p": 0,
                  "lam2Im": 0,
                  "lam2Re": 0,
                  "lam31": 0,
                  "lam31p": 0,
                  "lam33": 0,
                  "lam3Im": 0,
                  "lam3Re": 0,
                  "g1": 0,
                  "g2": 0,
                  "g3": 0}

        self.assertEqual(reference, bIsPerturbative(source))

    def test_bIsPerturbativeFalse(self):
        reference = False
        source = {"lam11": -999,
                  "lam12": 0,
                  "lam12p": 0,
                  "lam1Im": 0,
                  "lam1Re": 0,
                  "lam22": 0,
                  "lam23": 0,
                  "lam23p": 0,
                  "lam2Im": 0,
                  "lam2Re": 0,
                  "lam31": 0,
                  "lam31p": 0,
                  "lam33": 0,
                  "lam3Im": 0,
                  "lam3Re": 0,
                  "g1": 0,
                  "g2": 0,
                  "g3": 0}

        self.assertEqual(reference, bIsPerturbative(source))

