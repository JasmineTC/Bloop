import numpy as np
from math import sqrt, pi, log, exp
import scipy

def bIsPerturbative(paramDict4D : dict[str, float], pertSymbols: set) -> bool:
    ## Should actually check vertices but not a feature in DRalgo at time of writting
    for key, value in paramDict4D.items():
        if key in pertSymbols and abs(value) > 4*pi:
            return False
    return True

def constructSplineDictArray(betaFunction4DExpression, muRange, initialConditions, allSymbolsDict) :
    ## -----BUG------
    ## This updates the RGScale with the value of mu
    solutionSoft = scipy.integrate.odeint(lambda initialConditions, mu:  np.array(betaFunction4DExpression.evaluate(initialConditions))/mu,
                                          initialConditions, 
                                          muRange).transpose()

    interpDict = {}
    for key, value in allSymbolsDict.items():
        if key == "RGScale":
            continue
        
        ## Hack to remove all the const entries in the array
        if np.all(solutionSoft[value] == solutionSoft[value][0]):
            continue
        
        interpDict[key] =  scipy.interpolate.CubicSpline(muRange, solutionSoft[value], extrapolate = False)
        
    return interpDict

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
    betaFunction4DExpression: str = "betaFunction4DExpression" 
    
    bVerbose: bool = False
    
    EulerGammaPrime = 2.*(log(4.*pi) - np.euler_gamma)
    Lfconst = 4.*log(2.)
    
    allSymbolsDict: list = field(default_factory=dict)
    
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
        
        inputArray[self.allSymbolsDict["RGScale"]] = matchingScale
        inputArray[self.allSymbolsDict["T"]] = T
        inputArray[self.allSymbolsDict["Lb"]] = Lb
        inputArray[self.allSymbolsDict["Lf"]] = Lb + self.Lfconst
        return inputArray
    
    def updateParams4DRan(self, betaSpline4D: dict, array):
        muEvaulate = array[self.allSymbolsDict["RGScale"]]
        for key, spline in betaSpline4D.items():
            array[self.allSymbolsDict[key]] = float(spline(muEvaulate))
        return array
    
    def executeMinimisation(self, T, 
                            minimumLocation,
                            betaSpline4D): 
        paramValuesArray = self.updateTDependentConsts(T, np.zeros(len(self.allSymbolsDict.keys())))
        
        paramsForMatchingArray = self.updateParams4DRan(betaSpline4D, paramValuesArray)
        paramsForMatchingDict = {key: paramsForMatchingArray[value] for key, value in self.allSymbolsDict.items()}
        
        params3D = self.dimensionalReduction.getUltraSoftParams(paramsForMatchingDict, T)
        return ( *self.effectivePotential.findGlobalMinimum(T, params3D, self.initialGuesses + (minimumLocation, ) ), 
                bIsPerturbative(paramsForMatchingDict, self.pertSymbols), 
                bIsBounded(paramsForMatchingDict),
                params3D)
    
    def populateLagranianParams4D(self, inputParams: dict[str, float]) -> np.array:
        higgsVev = 246.22  #Consider using Fermi constant instead
        ## --- SM fermion and gauge boson masses---
        MW = 80.377 
        MZ = 91.1876 
        MTop = 172.76
        langrianParams4D = {"yt3": sqrt(2.) * MTop/ higgsVev,
                            "g1": 2.*sqrt(MZ**2 - MW**2)/ higgsVev, ## U(1)
                            "g2": 2.*MW/ higgsVev,                  ## SU(2)
                            "g3": sqrt(0.1183 * 4.0 * pi),          ## SU(3)
                            ## BSM stuff from benchmark
                            "RGScale":  inputParams["RGScale"],
                            **inputParams["massTerms"],
                            **inputParams["couplingValues"]}
        
        params4D = np.zeros(len(self.allSymbolsDict))
        for key, value in langrianParams4D.items():
            params4D[self.allSymbolsDict[key]] = value

        return params4D
            
    def traceFreeEnergyMinimum(self, benchmark:  dict[str: float]) -> dict[str: ]:
        lagranianParams4DArray = self.populateLagranianParams4D(benchmark)
               
        ## RG running. We want to do 4D -> 3D matching at a scale where logs are small; 
        ## usually a T-dependent scale 4.*pi*exp(-np.euler_gamma)*T 
        ## TODO FIX for when user RGscale < 7T!!!
        muRange = np.linspace(lagranianParams4DArray[self.allSymbolsDict["RGScale"]], 
                              7.3 * self.TRange[-1],
                              len(self.TRange)*10)
        
        betaSpline4D = constructSplineDictArray(self.betaFunction4DExpression, 
                                                muRange, 
                                                lagranianParams4DArray, 
                                                self.allSymbolsDict)
        
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
                                                                                                           betaSpline4D)
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

