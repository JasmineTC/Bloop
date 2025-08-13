from math import sqrt, pi, log, exp
import numpy as np
import scipy

from dataclasses import dataclass, InitVar,field

# from ThreeHiggs.BmGenerator import bIsBounded
from ThreeHiggs.PDGData import mTop, mW, mZ, higgsVEV

def bIsPerturbative(
    paramValuesArray, 
    pertSymbols, 
    allSymbols
):
    ## Should actually check vertices but not a feature in DRalgo at time of writting
    for pertSymbol in pertSymbols:
        if abs(paramValuesArray[allSymbols.index(pertSymbol)]) > 4*pi:
            return False

    return True

def constructSplineDictArray(
    betaFunction4DExpression, 
    muRange, 
    initialConditions, 
    allSymbols
):
    ## -----BUG------
    ## This updates the RGScale with the value of mu
    ## Weird FP errors can occur here (e.g. including mu in np.real or not)
    solutionSoft = scipy.integrate.solve_ivp(lambda mu, initialConditions: np.real(betaFunction4DExpression.evaluate(initialConditions)/mu),
                                             (muRange[0], muRange[-1]), 
                                             initialConditions, 
                                             t_eval=muRange).y
    interpDict = {}
    for idx, ele in enumerate(allSymbols):
        if ele == "RGScale":
            continue
        
        ## Hack to remove all the const entries in the array
        if np.all(solutionSoft[idx] == solutionSoft[idx][0]):
            continue
        
        interpDict[ele] =  scipy.interpolate.CubicSpline(muRange, solutionSoft[idx])
    return interpDict

@dataclass(frozen=True)
class TrackVEV:
    TRange: tuple = (0,)
    
    pertSymbols: frozenset = frozenset({1})
    
    initialGuesses: tuple = (0,)
    
    ## idk how to type hint this correctly
    effectivePotential: str = "effectivePotentialInstance"
    dimensionalReduction: str = "dimensionalReductionInstance"
    betaFunction4DExpression: str = "betaFunction4DExpression" 
    
    verbose: bool = False
    
    EulerGammaPrime = 2.*(log(4.*pi) - np.euler_gamma)
    Lfconst = 4.*log(2.)
    
    allSymbols: list = field(default_factory = list)
    
    config: InitVar[dict] = None
    
    def __post_init__(self, 
        config
    ):
        if config:
            self.__init__(**config)
    
    def trackVEV(self, 
        benchmark
    ):
        lagranianParams4DArray = self.getLagranianParams4D(benchmark)
               
        ## RG running. We want to do 4D -> 3D matching at a scale where logs are small; 
        ## usually a T-dependent scale 4.*pi*exp(-np.euler_gamma)*T 
        ## TODO FIX for when user RGscale < 7T!!!
        muRange = np.linspace(lagranianParams4DArray[self.allSymbols.index("RGScale")], 
                              7.3 * self.TRange[-1],
                              len(self.TRange)*10)
        
        betaSpline4D = constructSplineDictArray(self.betaFunction4DExpression, 
                                                muRange, 
                                                lagranianParams4DArray, 
                                                self.allSymbols)
        
        minimizationResults = {"T": [],
                               "vevDepthReal": [],
                               "vevDepthImag": [],
                               "vevLocation": [], 
                               "bIsPerturbative": [], 
                               "bIsBounded": []}

        counter = 0
        ## Initialise vevLocation to feed into the minimisation algo so it can
        ## use the location of the previous minimum as a guess for the next
        ## Not ideal as the code has to repeat an initial guess on first T
        vevLocation = np.array(self.initialGuesses[0])
        
        for T in self.TRange:
            if self.verbose:
                print (f'Start of temp = {T} loop')
            
            params4DRan = self.runParams4D(betaSpline4D, T)
            paramsUltraSoft = self.dimensionalReduction.softToUltraSoft.evaluate(
                                  self.dimensionalReduction.softScaleRGE.evaluate(
                                      self.dimensionalReduction.hardToSoft.evaluate(
                                      params4DRan
            )))
            
            vevLocation, vevDepth = self.effectivePotential.findGlobalMinimum(
               T, 
               paramsUltraSoft, 
               self.initialGuesses + [vevLocation]
            )
            
            minimizationResults["T"].append(T)
            minimizationResults["vevDepthReal"].append(vevDepth.real)
            minimizationResults["vevDepthImag"].append(vevDepth.imag)
            minimizationResults["vevLocation"].append(vevLocation)
            #bIsBounded(paramsForMatchingDict)
            minimizationResults["bIsPerturbative"].append(bIsPerturbative(
                params4DRan, 
                self.pertSymbols, 
                self.allSymbols)
            )

            if np.all( vevLocation < 0.5):
                if self.verbose:
                    print (f"Symmetric phase found at temp {T}")
                if counter == 3:
                    break
                counter += 1
            
        minimizationResults["vevLocation"] = np.transpose(minimizationResults["vevLocation"]).tolist()
        
        return minimizationResults
    
    def getLagranianParams4D(self, 
        inputParams
    ):
        ## --- SM fermion and gauge boson masses---
        ## How get g3 from PDG??
        langrianParams4D = {"yt3": sqrt(2.) * mTop/ higgsVEV,
                            "g1": 2.*sqrt(mZ**2 - mW**2)/ higgsVEV, ## U(1)
                            "g2": 2.*mW/ higgsVEV,                  ## SU(2)
                            "g3": sqrt(0.1183 * 4.0 * pi),          ## SU(3)
                            ## BSM stuff from benchmark
                            "RGScale":  inputParams["RGScale"],
                            **inputParams["massTerms"],
                            **inputParams["couplingValues"]}
        
        params4D = np.zeros(len(self.allSymbols))
        for key, value in langrianParams4D.items():
            params4D[self.allSymbols.index(key)] = value

        return params4D
    
    def getTConsts(self, 
        T, 
        inputArray
    ):
        matchingScale = 4.*pi*exp(-np.euler_gamma) * T
        Lb = 2. * log(matchingScale / T) - self.EulerGammaPrime
        
        inputArray[self.allSymbols.index("RGScale")] = matchingScale
        inputArray[self.allSymbols.index("T")] = T
        inputArray[self.allSymbols.index("Lb")] = Lb
        inputArray[self.allSymbols.index("Lf")] = Lb + self.Lfconst
        return inputArray
    
    def runParams4D(
        self, 
        betaSpline4D, 
        T
    ):
        array = self.getTConsts(T, np.zeros(len(self.allSymbols), dtype="float64"))
        
        muEvaulate = array[self.allSymbols.index("RGScale")]
        for key, spline in betaSpline4D.items():
            # Taking real part to avoid complex to real cast warning
            array[self.allSymbols.index(key)] = spline(muEvaulate)
        return array
    
    ############################
    def plotPotential(self, benchmark:  dict[str: float]):
        ## This is just a trimmed version of trace free energy minimum Jasmine uses for plotting
        lagranianParams4DArray = self.populateLagranianParams4D(benchmark)
               
        muRange = np.linspace(lagranianParams4DArray[self.allSymbols.index("RGScale")], 
                              7.3 * self.TRange[-1],
                              len(self.TRange)*10)
        
        betaSpline4D = constructSplineDictArray(self.betaFunction4DExpression, 
                                                muRange, 
                                                lagranianParams4DArray, 
                                                self.allSymbols)
        
        vevLocation = np.array(self.initialGuesses[0])
        
        linestyle = ["-.", "-", "--"]
        v3Max = 0
        yMin = 0
        yMax = 0
        for idx, T in enumerate(self.TRange):
            vevLocation, vevDepthReal, vevDepthImag, status, isPert, isBounded, params3D  = self.executeMinimisation(T,
                                                                                                           tuple(vevLocation.round(5)),                      
                                                                                                           betaSpline4D)
            # print(vevLocation)
            # v3Max = vevLocation[2] if vevLocation[2] > v3Max else v3Max
            # yMinMax = self.effectivePotential.plotPot(T, params3D, linestyle[idx], vevLocation[2], vevDepthReal, v3Max)
            self.effectivePotential.plotPot3D(T, params3D)
            # yMin = yMinMax[0] if yMinMax[0]< yMin else yMin
            # yMax = yMinMax[1] if yMinMax[1]> yMax else yMax

        # import matplotlib.pylab as plt
        # plt.legend(loc = 2)
        # plt.rcParams['text.usetex'] = True
        # plt.xlabel(r"$v_3$ ($\text{GeV}^{\: \frac{1}{2}})$", labelpad=-4)
        # plt.ylabel(r"$\dfrac{\Delta V}{T^3}$", rotation=0, labelpad = +12)
        # plt.hlines(0, 0, v3Max*1.1, colors = 'black')
        # plt.vlines(0, yMin*1.01, yMax*1.01,  colors = 'black')
        # plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        # plt.savefig("Results/StrongPot.png")
        # plt.show()
        return None
        

        




from unittest import TestCase
class TransitionFinderUnitTests(TestCase):
    def test_bIsPerturbativeTrue(self):
        reference = True
        source = [0.7, -0.8, 0]
        pertSymbols = {"lam11", "lam12", "lam12p"}
        allSymbols = ["lam11", "lam12", "lam12p"]

        self.assertEqual(reference, bIsPerturbative(source, pertSymbols, allSymbols) )

    def test_bIsPerturbativeFalse(self):
        reference = False
        source = [-12.57, 0, 0]
        pertSymbols = {"lam11", "lam12", "lam12p"}
        allSymbols = ["lam11", "lam12", "lam12p"]

        self.assertEqual(reference, bIsPerturbative(source, pertSymbols, allSymbols) )

