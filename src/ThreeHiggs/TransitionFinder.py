import numpy as np
from math import sqrt, pi, log, exp

def bIsPerturbative(paramDict4D : dict[str, float], pertSymbols: set) -> bool:
    ## Should actually check vertices but not a feature in DRalgo at time of writting
    for key, value in paramDict4D.items():
        if key in pertSymbols and abs(value) > 4*pi:
            return False
    return True

def runCoupling(betaSpline4D: dict[str, float], keyMapping: list[str], muEvaulate: float):
    runCoupling = {}
    for key in keyMapping:
        runCoupling[key] = betaSpline4D[key](muEvaulate).item()
    return runCoupling

def get4DLagranianParams(inputParams: dict[str, float]) -> dict[str, float]:
    
    langrianParams4D = {}
    ## --- SM fermions and gauge bosons ---
    v = 246.22  ## "Higgs VEV". Consider using Fermi constant instead
    langrianParams4D["yt3"] = sqrt(2.) * 172.76 / v  ## 172.76 is the top mass
    MW = 80.377 ## W boson mass
    MZ = 91.1876 ## Z boson mass
    
    langrianParams4D |= {"g1": 2.*sqrt(MZ**2 - MW**2)/ v,  # U(1)
            "g2": 2.*MW/ v,                   # SU(2)
            "g3": sqrt(0.1183 * 4.0 * pi)}    # SU(3)
    
    ## --- BSM scalars ---
    langrianParams4D |= inputParams["couplingValues"] | inputParams["massTerms"]
    langrianParams4D["RGScale"] = inputParams["RGScale"]

    return langrianParams4D

from dataclasses import dataclass, InitVar,field
from ThreeHiggs.BmGenerator import bIsBounded
@dataclass(frozen=True)
class TraceFreeEnergyMinimum:
    TRange: tuple = (0,)
    
    pertSymbols: frozenset = frozenset({1})
    
    initialGuesses: tuple = (0,)
    
    ## Hack - idk how to type hint this correctly
    effectivePotential: str = "effectivePotentialInstance"
    dimensionalReduction: str = "dimensionalReductionInstance"
    
    bVerbose: bool = False
    
    EulerGammaPrime = 2.*(log(4.*pi) - np.euler_gamma)
    Lfconst = 4.*log(2.)
    
    index2Arg: dict = field(default_factory=dict)
    arg2Index: dict = field(default_factory=dict)
    
    config: InitVar[dict] = None
    
    def __post_init__(self, config: dict):
        if config:
            self.__init__(**config)
            
    def isBad(self, T, 
               minimumLocation, 
               status):
        ## This is a hack to remove bad benchmark points
        if T == self.TRange[0] and (minimumLocation[0] > 1 or minimumLocation[1] > 1):
            return "v3NotGlobalMin"
        if status == "NaN": 
            return "MinimisationFailed"
        return False
    
    
    def updateTDependentConsts(self, T, inputArray):
        matchingScale = 4.*pi*exp(-np.euler_gamma) * T
        Lb = 2. * log(matchingScale / T) - self.EulerGammaPrime
        
        inputArray[self.arg2Index["RGScale"]] = matchingScale
        inputArray[self.arg2Index["T"]] = T
        inputArray[self.arg2Index["Lb"]] = Lb
        inputArray[self.arg2Index["Lf"]] = Lb + self.Lfconst
        return inputArray
    
    def executeMinimisation(self, T, 
                            minimumLocation,
                            betaSpline4D,
                            keyMapping):
        
        paramValuesArray = self.updateTDependentConsts(T, np.zeros(len(self.arg2Index.keys())))
        paramValuesDict = {key: paramValuesArray[value] for key, value in self.arg2Index.items()}
        
        paramsForMatching = paramValuesDict | runCoupling(betaSpline4D, 
                                                          keyMapping, 
                                                          paramValuesDict["RGScale"])
        
        
        
        params3D = self.dimensionalReduction.getUltraSoftParams(paramsForMatching, T)
        a = self.initialGuesses + (minimumLocation, ) 
        return ( *self.effectivePotential.findGlobalMinimum(T, params3D, a), 
                bIsPerturbative(paramsForMatching, self.pertSymbols), 
                bIsBounded(paramsForMatching),
                params3D)
    
    def populateLagranianParams4D(self, inputParams: dict[str, float], argArray: np.array) -> np.array:
        
        langrianParams4D = {}
        ## --- SM fermions and gauge bosons ---
        v = 246.22  ## "Higgs VEV". Consider using Fermi constant instead
        langrianParams4D["yt3"] = sqrt(2.) * 172.76 / v  ## 172.76 is the top mass
        MW = 80.377 ## W boson mass
        MZ = 91.1876 ## Z boson mass
        
        langrianParams4D |= {"g1": 2.*sqrt(MZ**2 - MW**2)/ v,  # U(1)
                "g2": 2.*MW/ v,                   # SU(2)
                "g3": sqrt(0.1183 * 4.0 * pi)}    # SU(3)
        
        ## --- BSM scalars ---
        langrianParams4D |= inputParams["couplingValues"] | inputParams["massTerms"]
        langrianParams4D["RGScale"] = inputParams["RGScale"]
        
        for key, value in langrianParams4D.items():
            argArray[self.arg2Index[key]] = value

        return argArray
            
    def traceFreeEnergyMinimum(self, benchmark:  dict[str: float]) -> dict[str: ]:
        # lagranianParams4D = get4DLagranianParams(benchmark)
        lagranianParams4DArray = self.populateLagranianParams4D(benchmark, np.zeros(len(self.arg2Index)))
        
        print(lagranianParams4D)
        print()
        print(lagranianParams4DArray)
        exit()
        ## RG running. We want to do 4D -> 3D matching at a scale where logs are small; 
        ## usually a T-dependent scale ~7.3T
        muRange = np.linspace(lagranianParams4D["RGScale"], 
                              7.3 * self.TRange[-1],
                              len(self.TRange)*10)
        
        from .BetaFunctions import BetaFunctions4D
        betasFunctions = BetaFunctions4D() 
        betaSpline4D, keyMapping = betasFunctions.constructSplineDict(muRange, lagranianParams4D)
        
        minimizationResults = {"T": [],
                               "valueVeffReal": [],
                               "valueVeffImag": [],
                               "complex": False,
                               "minimumLocation": [], 
                               "bIsPerturbative": True, 
                               "UltraSoftTemp": None, 
                               "failureReason": None}

        counter = 0
        ## Initialise minimumLocation to feed into the minimisation algo so it can
        ## use the location of the previous minimum as a guess for the next
        ## Not ideal as the code has to repeat an initial guess on first T
        minimumLocation = np.array(self.initialGuesses[0])
        for T in self.TRange:
            if self.bVerbose:
                print (f'Start of temp = {T} loop')
                
            minimizationResults["T"].append(T)
            
            minimumLocation, minimumValueReal, minimumValueImag, status, isPert, isBounded, params3D  = self.executeMinimisation(T,
                                                                                                           tuple(minimumLocation.round(5)),                      
                                                                                                           betaSpline4D, 
                                                                                                           keyMapping)
            ##Not ideal name or structure imo
            isBadState = self.isBad(T, minimumLocation, status)
            if isBadState:
                minimizationResults["failureReason"] = isBadState
                break
            
            
            minimizationResults["valueVeffReal"].append(minimumValueReal)
            minimizationResults["valueVeffImag"].append(minimumValueImag)
            minimizationResults["minimumLocation"].append(minimumLocation)
            
            if not minimizationResults["complex"]:
                minimizationResults["complex"] = (status == "complex")
            
            if not minimizationResults["UltraSoftTemp"]:
                if self.effectivePotential.bReachedUltraSoftScale(minimumLocation,
                                                                  T, 
                                                                  params3D): 
                    minimizationResults["UltraSoftTemp"] = T

                
            if np.all( minimizationResults["minimumLocation"][-1] < 0.5):
                if self.bVerbose:
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
        source = {"lam11": 0.7,
                  "lam12": -0.8,
                  "lam12p": 0}
        symbols = {"lam11", "lam12", "lam12p"}

        self.assertEqual(reference, bIsPerturbative(source, symbols))

    def test_bIsPerturbativeFalse(self):
        reference = False
        source = {"lam11": -12.57,
                  "lam12": 0,
                  "lam12p": 0}
        symbols = {"lam11", "lam12", "lam12p"}

        self.assertEqual(reference, bIsPerturbative(source, symbols) )

