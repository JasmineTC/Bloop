import numpy as np
import numpy.typing as npt
from typing import Tuple
import scipy.optimize
from dataclasses import dataclass

import pathlib ## for hacking

import Integrals

from ParsedExpression import ParsedSystemOfEquations

""" Helper "struct", just holds mass eigenvalues squared. This is convenient to have for storing both scalar and gauge masses 
in one place, while still having names for each mass. Not very model independent though. (Would a dict be better??) 
"""  
@dataclass 
class MassSquared:
    mWsq: float
    mZsq: float
    


""" Evaluating the potential: 
1. Call setModelParameters() with a dict that sets all parameters in the action. 
This is assumed to be using 3D EFT, so the params are temperature dependent.
2. Call evaluate() with a list that specifies values of background fields. Fields are in 3D units, ie. have dimension GeV^(1/2)
"""
class EffectivePotential:

    ## Defining equations for field-dependent mixing angles
    neutralDiagonalizationEquations: ParsedSystemOfEquations
    chargedDiagonalizationEquations: ParsedSystemOfEquations


    def __init__(self, initialModelParameters: dict = None):
        
        if (initialModelParameters):
            self.setModelParameters(initialModelParameters)

        ## Initialize diagonalization equations etc
        self._initAngleEquations()
        


    def setModelParameters(self, modelParameters: dict) -> None:
        """ This just reads action parameters from a dict and sets them internally for easier/faster(?) access in evaluate 
        """

        self.g1sq = modelParameters["g1sq"]
        self.g2sq = modelParameters["g2sq"]
        ## QCD coupling g3 not needed

        self.mu1sq = modelParameters["mu1sq"]
        self.mu2sq = modelParameters["mu2sq"]
        self.mu3sq = modelParameters["mu3sq"]

        self.lam11 = modelParameters["lam11"]
        self.lam22 = modelParameters["lam22"]
        self.lam33 = modelParameters["lam33"]
        self.lam12 = modelParameters["lam12"]
        self.lam23 = modelParameters["lam23"]
        self.lam31 = modelParameters["lam31"]
        ## Primed params
        self.lam12p = modelParameters["lam12p"]
        self.lam23p = modelParameters["lam23p"]
        self.lam31p = modelParameters["lam31p"]

        self.mu12sqRe = modelParameters["mu12sqRe"]
        self.mu12sqIm = modelParameters["mu12sqIm"]
        self.lam1Re = modelParameters["lam1Re"]
        self.lam1Im = modelParameters["lam1Im"]
        self.lam2Re = modelParameters["lam2Re"]
        self.lam2Im = modelParameters["lam2Im"]
        self.lam3Re = modelParameters["lam3Re"]
        self.lam3Im = modelParameters["lam3Im"]

        self.modelParameter = modelParameters ## Dunno if we want to have this tbh, but store for now


    def calcMassEigenvalues(self, fields: list[float]) -> None:

        fieldsArray = np.asanyarray(fields)
        ## "Full vev^2" in 3HDM:
        vevSq = np.sum(fieldsArray**2) 

        ## Gauges
        mWsq = 0.25 * vevSq * self.g2sq
        mZsq = 0.25 * vevSq * (self.g1sq + self.g2sq)

        ## TODO scalars

        return MassSquared(mWsq=mWsq, mZsq=mZsq)
    

    def calcAngles(self, fields: list[float]):
        ## TODO. do something like:
        #sin1, sin2, sin3, sin4, sin5, sin6 = solveEquationSystem( [eq1, eq2, eq3, eq4, eq5, eq6], fields, initialGuessList)
        # angle1 = ...

        
        return 1, 2, 3, 4, 5, 6


    ## Evaluate Veff(fields, T) using current set of model parameters
    def evaluate(self, fields: list[float], loopOrder: int = 2) -> complex:
        
        VTotal = self.V0(fields)

        if (loopOrder > 0):

            
            angles = self.calcAngles(fields)
            
            massSquared: MassSquared = self.calcMassEigenvalues(fields)

            ## 3D Coleman-Weinberg correction in Landau gauge. Gauge fields get overall (d-1) = 2
            V1 = 4. * Integrals.J3(massSquared.mWsq) + 2* Integrals.J3(massSquared.mZsq)

            ## TODO scalars, vectorize with numpy if using a loop

       

        ## TODO 2-loop

        return VTotal
    

    ## Return value is location, value
    def findLocalMinimum(self, initialGuess: list[float]) -> Tuple[list[float], complex]:
        
        ## I think we need to manually vectorize here if our parameters are arrays (ie. got multiple temperature inputs).
        ## Then self.evaluate would return a numpy array which scipy doesn't know how to work with. 
        ## Here I make an array of lambda functions and minimize those separately


        ## Minimize real part only:
        VeffWrapper = lambda fields: np.real ( self.evaluate(fields) )

        res = scipy.optimize.minimize(VeffWrapper, initialGuess)

        ## res.x = location, res.fun = value, res.success = flag for determining if the algorithm finished successfully
        location = res.x

        ## evaluate once more to get the possible imag parts
        value = self.evaluate(location)
        return location, value
    

    def findGlobalMinimum(self, minimumCandidates: list[list[float]] = None) -> Tuple[list[float], complex]:
        """This calls findLocalMinimum with a bunch of initial guesses and figures out the deepest solution.
        Generally will not work very well if no candidates minima are given. 
        Return value is location, value. value can be complex (but this is probably a sign of failed minimization)
        """
        
        if (minimumCandidates == None or len(minimumCandidates) == 0):
            ## Just search "symmetric" and "broken" in the phi3 direction
            minimumCandidates = [ [0., 0., 1e-4], [0., 0., 20.] ]

        deepest = np.inf
        res = np.nan, np.inf

        ## Should we vectorize??
        for candidate in minimumCandidates:
            location, value = self.findLocalMinimum(candidate)
            if (np.real(value) < deepest):
                res = location, value

        return res


    # just calls self.evaluate
    def __call__(self, temperature: npt.ArrayLike, fields: npt.ArrayLike) -> complex:
        self.evaluate(temperature, fields)


    ## Tree level potential
    def V0(self, fields: npt.ArrayLike) -> float:
        _, _, v3 = fields 

        # this assumes just one background field (in the phi3 doublet)

        res = -0.5*self.mu3sq * v3**2 + 0.25*self.lam33*v3**4
        return res
    

    def _initAngleEquations(self) -> None:

        ## TODO read these paths from a config or something. This is a hack
        pathToCurrentFile = pathlib.Path(__file__).parent.resolve()
        neutralAngleFile = str(pathToCurrentFile) + "/Data/EffectivePotential/neutralDiagonalizationAnglesEquations.txt"
        chargedAngleFile = str(pathToCurrentFile) + "/Data/EffectivePotential/chargedDiagonalizationAnglesEquations.txt"

        self.neutralDiagonalizationEquations = ParsedSystemOfEquations(neutralAngleFile)
        print("-- Parsed neutral scalar angle equations, symbols:")
        print(self.neutralDiagonalizationEquations.functionArguments)

        self.chargedDiagonalizationEquations = ParsedSystemOfEquations(chargedAngleFile)
        print("-- Parsed charged scalar angle equations, symbols:")
        print(self.chargedDiagonalizationEquations.functionArguments)


    ## This is for putting things in correct order for angle equation lambdas. But very ugly and WIP 
    def _wrapAngleEquations(self, unknownNames: list[str], fields: list[float], params: dict[str, float]) -> None:

        ## Hack this for now. This is probably not very efficient
        _, _, v3 = fields
        newDict = params.copy() 
        newDict["v3"] = v3
        
        def wrapper(s1, s2, s3, s4, s5, s6):

            tempList = [s1, s2, s3, s4, s5, s6]
            ## .....

