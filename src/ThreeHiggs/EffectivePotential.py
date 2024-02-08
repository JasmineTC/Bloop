import numpy as np
import numpy.typing as npt
from typing import Tuple
import scipy.optimize
from dataclasses import dataclass
import matplotlib.pylab as plt

from scipy.optimize import least_squares
from scipy.interpolate import CubicSpline

import pathlib ## for hacking

from .Integrals import J3

from .ParsedExpression import ParsedExpressionSystem, SystemOfEquations, ParsedExpression
from .CommonUtils import combineDicts

from .VeffMinimizer import VeffMinimizer

class MixingAngleEquations(SystemOfEquations):
    """This is just a SystemOfEquations but we add a common solve() routine that can be configured for solving sines of angles, 
    and then returning the angles."""

    def __init__(self, fileName: str, unknownVariables: list[str]):
        super().__init__(fileName, unknownVariables)


    def solve(self, knownVariables: dict[str, float]) -> list[float]:

        ## Get list of "extra" variables in correct order
        otherArgs = tuple( self.getOtherVariablesFromDict(knownVariables) )

        ## Initial guesses for the sines
        initialGuess = 0., 0., 0., 0., 0., 0.

        ## function signature needs to be f(x, *args) if x are the unknowns.
        ## This is now a very stupid implementation, too much list <-> tuple <-> numpy array etc conversions. TODO need to clean this up
        def evaluateWrapper(x: npt.ArrayLike, *args):
            return self.evaluateSystem( x.tolist() + list(args) )
        
        ##Least squares will minimise a function given an initial guess and some bounds (needed so abs(S_i) <= 1)
        ##ftol controls the difference in function evalution, xtol how small the change in arguement can be, gtol how much the gradient changes
        ##the default tol of 1e-8 was not sufficient for outputs to agree on different machines
        ##Changed the jacobian calculation to use 3 point method, twice as many computations but should be more accurate
        #res = (least_squares(evaluateWrapper, initialGuess, bounds=(-1,1), args=(otherArgs), ftol = 1e-9, xtol=1e-9, gtol=1e-9, jac='3-point')).x
        res = (least_squares(evaluateWrapper, initialGuess, bounds=(-1,1), args=(otherArgs))).x
        ## The solver can throw warnings if it tries to evaluate the eqs outside sine's value range. So TODO change this to some algorithm that can limit the search range 

        ## We solved sines, not angles, but this is OK for now at least. So just return the sines.
    
        ## return list, not numpy array for now
        return res.tolist()
    
    def solveAsDict(self, knownVariables: dict[str, float]) -> dict[str, float]:
        """Same as solve() but returns a dict with named output.
        """

        solns = self.solve(knownVariables)
        return dict(zip( self.unknownVariables, solns ))
        




## everything we need to evaluate the potential. this is very WIP
class VeffParams:
    """ Usage after initialization: 
        1. call setActionParams(dict)
        2. call evaluateAll
    """

    """ Order of computations:
        1. Fix action params (and possibly temperature). This needs to be done via setActionParams() before evaluateAll(fields) 
        2. Fix background fields 
        3. Solve diagonalization conditions
        4. Evaluate post-diagonalization shorthand symbols
        5. Evaluate masses
        6. Evaluate pre-Veff shorthand symbols (so stuff that does into rotation matrices)
    """

    fields: dict[str, float]
    actionParams: dict[str, float]
    angles: dict[str, float]
    masses: dict[str, float]

    ## TODO export these from mathematica and read in here
    fieldNames: list[str] = ["v3"]

    """Equations that need to be solved numerically to guarantee diagonal mass matrix (at tree level).
    We allow there to be many independent systems of equations.
    """
    diagonalizationConditions: list[MixingAngleEquations]
    
    ## Field-dependent masses
    vectorMassesSquared: ParsedExpressionSystem
    scalarMassesSquared: ParsedExpressionSystem

    ## Shorthand symbols
    postDiagonalizationShorthands: ParsedExpressionSystem
    preVeffShorthands: ParsedExpressionSystem

    def __init__(self):
        self.initMassExpressions()
        self.initDiagonalizationConditions()
        self.initShorthandSymbols()


    def setActionParams(self, inputParams: dict[str, float]) -> None:
        self.actionParams = inputParams


    def evaluateAll(self, fields: list[float], bNeedsDiagonalization=True) -> dict[str, float]:
        """This should return a dict that fixes all symbols needed for Veff 2-loop evaluation.
        """
        ## To make this work nicely we need to put the field in same dict as our other inputs, so hack it here (will need optimization)
        _, _, v3 = fields
        #v3 = fields
        knownParamsDict = self.actionParams.copy()
        knownParamsDict["v3"] = v3
        
        if (bNeedsDiagonalization):
            ## Diagonalization conditions
            diagCondDict = self.solveDiagonalizationConditions(knownParamsDict)
            knownParamsDict = combineDicts(knownParamsDict, diagCondDict)

            ## Post-diagonalization symbols
            shorthands = self.postDiagonalizationShorthands.evaluateSystemWithDict(knownParamsDict, bReturnDict=True)
            knownParamsDict = combineDicts(knownParamsDict, shorthands)

            ## Masses
            masses = self.evaluateMasses(knownParamsDict)
            knownParamsDict = combineDicts(knownParamsDict, masses)

            ## pre-Veff symbols
            shorthands = self.preVeffShorthands.evaluateSystemWithDict(knownParamsDict, bReturnDict=True)
            knownParamsDict = combineDicts(knownParamsDict, shorthands)

        return knownParamsDict
    

    def evaluateScalarMasses(self, knownParams: dict[str, float]) -> dict[str, float]:
        scalarMassesSquaredDict =  self.scalarMassesSquared.evaluateSystemWithDict(knownParams, bReturnDict=True)
        ## Take abs. TODO handle better
        for key in scalarMassesSquaredDict:
            scalarMassesSquaredDict[key] = abs(scalarMassesSquaredDict[key])
        return scalarMassesSquaredDict

    def evaluateVectorMasses(self, knownParams: dict[str, float]) -> dict[str, float]:
        vectorMassesSquaredDict = self.vectorMassesSquared.evaluateSystemWithDict(knownParams, bReturnDict=True)
        ## Take abs. TODO handle better
        for key in vectorMassesSquaredDict:
            vectorMassesSquaredDict[key] = abs(vectorMassesSquaredDict[key])
        return vectorMassesSquaredDict

    def evaluateMasses(self, knownParams: dict[str, float]) -> dict[str, float]:
        return combineDicts( self.evaluateScalarMasses(knownParams), self.evaluateVectorMasses(knownParams) )


    def solveDiagonalizationConditions(self, knownParams: dict[str, float]) -> dict[str, float]:
        assert self.diagonalizationConditions != None
        
        solutionsDict = {}
        for cond in self.diagonalizationConditions:
             
            solutionsDict = combineDicts( solutionsDict, cond.solveAsDict(knownParams) )
        
        return solutionsDict


    def initMassExpressions(self):

        ## TODO read these paths from a config or something. This is a hack
        pathToCurrentFile = pathlib.Path(__file__).parent.resolve()
        scalarMassesFile = pathToCurrentFile / "Data/EffectivePotential/scalarMasses.txt"
        vectorMassesFile = pathToCurrentFile / "Data/EffectivePotential/vectorMasses.txt"
        
        self.vectorMassesSquared = ParsedExpressionSystem(vectorMassesFile)
        self.scalarMassesSquared = ParsedExpressionSystem(scalarMassesFile)


    def initDiagonalizationConditions(self):

        self.diagonalizationConditions = []

        ## TODO read these paths from a config or something. This is a hack
        pathToCurrentFile = pathlib.Path(__file__).parent.resolve()
        neutralAngleFile = pathToCurrentFile / "Data/EffectivePotential/neutralDiagonalizationAnglesEquations.txt"
        chargedAngleFile = pathToCurrentFile / "Data/EffectivePotential/chargedDiagonalizationAnglesEquations.txt"

        ## Currently we have to explicitly list the unknown variables by name:

        self.diagonalizationConditions.append( MixingAngleEquations(neutralAngleFile, ['S1Ne', 'S2Ne', 'S3Ne', 'S4Ne', 'S5Ne', 'S6Ne']) )

        self.diagonalizationConditions.append( MixingAngleEquations(chargedAngleFile, ['S1Ch', 'S2Ch', 'S3Ch', 'S4Ch', 'S5Ch', 'S6Ch']) )


    def initShorthandSymbols(self):

        ## TODO read these paths from a config or something. This is a hack
        pathToCurrentFile = pathlib.Path(__file__).parent.resolve()
        postDiagonalizationShorthandsFile = pathToCurrentFile / "Data/EffectivePotential/postDiagonalizationShorthands.txt"
        preVeffShorthandsFile = pathToCurrentFile / "Data/EffectivePotential/preVeffShorthands.txt"

        self.postDiagonalizationShorthands = ParsedExpressionSystem(postDiagonalizationShorthandsFile)
        self.preVeffShorthands = ParsedExpressionSystem(preVeffShorthandsFile)


""" Evaluating the potential: 
1. Call setModelParameters() with a dict that sets all parameters in the action. 
This is assumed to be using 3D EFT, so the params are temperature dependent.
2. Call evaluate() with a list that specifies values of background fields. Fields are in 3D units, ie. have dimension GeV^(1/2)
"""
class EffectivePotential:

    
    loopOrder: int
    ## One expression for each loop order
    expressions: ParsedExpressionSystem

    minimizer: VeffMinimizer

    def __init__(self, loopOrder, initialModelParameters: dict = None):
        """loopOrder specifies how many perturbative orders we take. 
        This CANNOT be changed at runtime, you will have to make a new object instead.
        Order count starts from 0. Tree level is 0, 1-loop is 1 etc.
        """

        self.loopOrder = loopOrder

        self.initExpressions()

        self.bNeedsDiagonalization = (loopOrder > 0)

        self.params = VeffParams()
        
        if (initialModelParameters):
            self.setModelParameters(initialModelParameters)

        self.minimizer = VeffMinimizer(3) # currently the numVariables is not used by minimizer

        


    def setModelParameters(self, modelParameters: dict) -> None:
        """ This just reads action parameters from a dict and sets them internally for easier/faster(?) access in evaluate 
        """

        self.params.setActionParams(modelParameters)

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
        """

    def initExpressions(self):

        self.expressions = []
        
        ## TODO read these paths from a config or something. This is a hack
        pathToCurrentFile = pathlib.Path(__file__).parent.resolve()

        veffFiles = []
        veffFiles.append( pathToCurrentFile / "Data/EffectivePotential/Veff_LO.txt")
        if (self.loopOrder >= 1):
            veffFiles.append( pathToCurrentFile / "Data/EffectivePotential/Veff_NLO.txt")
        if (self.loopOrder >= 2):
            veffFiles.append( pathToCurrentFile / "Data/EffectivePotential/Veff_NNLO.txt")


        ## Hack: combine these into a one file so that ParsedExpressionSystem understand it

        tempFileName = pathToCurrentFile / "tempFile1424522343.txt"

        with open(tempFileName, 'w') as tempf:
            for filename in veffFiles:
                with open(filename, 'r') as f:
                    content = f.read()
                    tempf.write(content)
                    tempf.write("\n")

        self.expressions = ParsedExpressionSystem(tempFileName)

    def evaluate(self, fields: list[float]) -> complex:
        """Evaluate Veff at specified field values. Uses the currently set model parameters.
        """
        
        self.params.fields = fields
    
        ## This has masses, angles, all shorthand symbols etc. Everything we need to evaluate loop corrections
        paramDict = self.params.evaluateAll(fields, bNeedsDiagonalization=self.bNeedsDiagonalization)
        
        #print(f"{paramDict=}")
        #input()

        ## summing works because the result is a list [V0, V1, ...]
        res = sum( self.expressions.evaluateSystemWithDict(paramDict) )
        return res



    ## Return value is location, value
    def findLocalMinimum(self, initialGuess: list[float]) -> Tuple[list[float], complex]:
        
        ## I think we need to manually vectorize here if our parameters are arrays (ie. got multiple temperature inputs).
        ## Then self.evaluate would return a numpy array which scipy doesn't know how to work with. 
        ## Here I make an array of lambda functions and minimize those separately

        ## Minimize real part only:
        VeffWrapper = lambda fields: np.real ( self.evaluate(fields) )

        ##Added bounds to minimize to reduce the IR senstivity coming from low mass modes
        bounds = ((1e-6, 1e-6), (1e-6, 1e-6), (1e-6, 1e3))
        #res = scipy.optimize.minimize(VeffWrapper, initialGuess, tol = 1e-8, bounds=bnds)
        location, value = self.minimizer.minimize(VeffWrapper, initialGuess, bounds)
        #print (f"for an initial guess of {initialGuess} the local minimum found is {location}")


        if np.any(np.isnan(location)):
            location  = [np.nan] * len[initialGuess]
            location = np.asarray(location)
            value = np.nan
        else:
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
            minimumCandidates = [ [0., 0., 1e-4], [0., 0., 30.] ]

        deepest = np.inf
        res = [np.nan]*3, np.inf

        ## Should we vectorize??
        for candidate in minimumCandidates:
            location, value = self.findLocalMinimum(candidate)
            if (np.real(value) < deepest):
                res = location, value
                deepest = np.real(value)

        return res


    # just calls self.evaluate
    def __call__(self, temperature: npt.ArrayLike, fields: npt.ArrayLike) -> complex:
        self.evaluate(temperature, fields)
