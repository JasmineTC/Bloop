import numpy as np
from math import sqrt, pi, log, exp

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


def get4DLagranianParams(inputParams: dict[str, float]) -> dict[str, float]:
    """Take inputs from the BM file.
    With tree-level matching the renormalization scale does not directly show up in the expressions, but
    needs to be specified for later loop calculations."""
    langrianParams4D = {}
    ## --- SM fermions and gauge bosons ---
    v = 246.22  ## "Higgs VEV". Consider using Fermi constant instead
    langrianParams4D["yt3"] = sqrt(2.) * 172.76 / v  ## 172.76 is the top mass (not squared! would it be better to store y^2?)
    MW = 80.377 ## W boson mass
    MZ = 91.1876 ## Z boson mass
    
    langrianParams4D |= {"g1": 2.*sqrt(MZ**2 - MW**2)/ v,  # U(1)
            "g2": 2.*MW/ v,                   # SU(2)
            "g3": sqrt(0.1183 * 4.0 * pi)}    # SU(3)
    
    ## --- BSM scalars ---
    langrianParams4D |= inputParams["couplingValues"] | inputParams["massTerms"]
    langrianParams4D["RGScale"] = inputParams["RGScale"]

    return langrianParams4D

def runCoupling(dictt, keyMapping, muEvaulate: float):
    runCoupling = {}
    ##Loop over the dict, for each key compute the spline at muEvaulate, ignore the RGScale in the dict
    for key in keyMapping:
        ##RHS is an array of a single number, so add item() to extract that single number
        runCoupling[key] = dictt[key](muEvaulate).item()
    return runCoupling

def traceFreeEnergyMinimum(effectivePotential,
                           dimensionalReduction,
                           benchmark,
                           TRangeStart: float, 
                           TRangeEnd: float, 
                           TRangeStepSize: float,
                           bVerbose = False) -> dict[str: ]:
    """RG running. We want to do 4D -> 3D matching at a scale where logs are small; usually a T-dependent scale ~7T.
    To make this work nicely we integrate the beta functions here up to the largest temp used 
    then interpolate over the beta function."""
    TRange = np.arange(TRangeStart, TRangeEnd, TRangeStepSize)
    
    LagranianParams4D = get4DLagranianParams(benchmark)
    startScale = LagranianParams4D["RGScale"]
    endScale = 7.3 * TRange[-1] 
    muRange = np.linspace(startScale, endScale, TRange.size*10)
    
    from .BetaFunctions import BetaFunctions4D
    betasFunctions = BetaFunctions4D() 
    BetaSpline4D, keyMapping = betasFunctions.constructSplineDict(muRange, LagranianParams4D)
    
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
    from ThreeHiggs.BmGenerator import bIsBounded
    for T in TRange:
        T = float(T) ## To make compatible with JSON
        if bVerbose:
            print (f'Start of temp = {T} loop')
        minimizationResults["T"].append(T)
        
        goalRGScale =  T ## Final scale in 3D -check if goalRGscale is ever different from just T
        matchingScale = 4.*pi*exp(-np.euler_gamma) * T ## Scale that minimises T dependent logs
        paramsForMatching = runCoupling(BetaSpline4D, keyMapping, matchingScale)
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
        params3D = dimensionalReduction.getUltraSoftParams(paramsForMatching, goalRGScale)
        
        initialGuesses = [[0.1,0.1,0.1], ## TODO This should go in a config file or something
                          [5,5,5],
                          [-5,5,5], 
                          [5,5,5],
                          [-5,5,5], 
                          [40,40,40],
                          [-40,40,40], 
                          [59,59,59], 
                          [-59,59,59]]
        minimumLocation, minimumValueReal, minimumValueImag, status = effectivePotential.findGlobalMinimum(T, params3D, initialGuesses)
        if T == TRangeStart and (minimumLocation[0] > 1 or minimumLocation[1] > 1): ## This is a hack to remove bad benchmark points
            minimizationResults["failureReason"] = "v3NotGlobalMin"
            break
            
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
                                                         params3D): 
                minimizationResults["UltraSoftTemp"] = T ## If reached then set that as the ultra soft temp
                                                         ##- this will stop the first if statement from passing
        
        if minimizationResults["bIsPerturbative"]: ##If the potential was perturbative check if it still is
            minimizationResults["bIsPerturbative"] = bIsPerturbative(paramsForMatching) ## If non-pert then value set to false and won't be updated

        if np.all(minimumLocation < 1e-2):
            if bVerbose:
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

