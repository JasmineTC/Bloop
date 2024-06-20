import numpy as np
import numpy.typing as npt
from typing import Tuple
import tempfile # hacking

from scipy import linalg
from dataclasses import dataclass

from .parsedmatrix import ParsedMatrix, MatrixDefinitionFiles
from .ParsedExpression import ParsedExpressionSystem, SystemOfEquations
from .CommonUtils import combineDicts, diagonalizeSymmetric

from .VeffMinimizer import VeffMinimizer

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

    def __init__(self, 
                 fieldNames, 
                 bAbsoluteMsq, 
                 vectorMassFile, 
                 vectorShorthandFile, 
                 scalarPermutationMatrix, 
                 scalarMassMatrices, 
                 scalarRotationMatrixFile):
        self.fieldNames = fieldNames
        self.bAbsoluteMsq = bAbsoluteMsq

        self.vectorMassesSquared = ParsedExpressionSystem(vectorMassFile)
        self.vectorShorthands = ParsedExpressionSystem(vectorShorthandFile)
        self.scalarPermutationMatrix = scalarPermutationMatrix
        ## can have many matrices if we've block-diagonalized already
        ## ASSUME: the blocks are given in order: upper left to lower right. 
        ##TODO improve this
        self.scalarMassMatrices = [ParsedMatrix(matrix.matrixFile, matrix.expressionsFile) for matrix in scalarMassMatrices]
        self.scalarRotationMatrix = ParsedMatrix(scalarRotationMatrixFile)

    def setActionParams(self, inputParams: dict[str, float]) -> None:
        self.actionParams = inputParams


    def evaluateAll(self, fields: list[float], bNeedsDiagonalization=True) -> dict[str, float]:
        """This should return a dict that fixes all symbols needed for Veff 2-loop evaluation.
        """
        # Gradually build a dict containing (key, val) for all needed symbols
        # TODO this will need to be optimized
        knownParamsDict = self.actionParams.copy()

        ## Background fields
        for i, value in enumerate(fields):
            knownParamsDict[self.fieldNames[i]] = value

        ## Vectors
        vectorShorthands = self.vectorShorthands.evaluateSystemWithDict(knownParamsDict, bReturnDict=True)
        knownParamsDict = combineDicts(knownParamsDict, vectorShorthands)
        vectorMasses = self.vectorMassesSquared.evaluateSystemWithDict(knownParamsDict, bReturnDict=True)

        if (self.bAbsoluteMsq):
            for key, val in vectorMasses.items():
                vectorMasses[key] = np.abs(val)

        knownParamsDict = combineDicts(knownParamsDict, vectorMasses)

        ## Scalars
        diagDict = self.diagonalizeScalars(knownParamsDict)
        
        knownParamsDict = combineDicts(knownParamsDict, diagDict)

        return knownParamsDict

    def diagonalizeScalars(self, params: dict[str, float]) -> dict[str, float]:
        """Finds a rotation matrix that diagonalizes the scalar mass matrix
        and returns a dict with diagonalization-specific params
        """
        # TODO optimize, comment this with some matrix eqs

        outDict = {}

        # Diagonalize blocks separately
        blockRot = []
        blockM = []
        for m in self.scalarMassMatrices:
            
            numericalM = m.evaluateWithDict(params)
            eigenValue, vects = diagonalizeSymmetric( numericalM, bCheckFinite=False )
            ## NOTE: vects has the eigenvectors on columns => D = V^T . M . V is diagonal

            ## Quick check that the numerical mass matrix is diagonal after being rotated by vects
            diagonalBlock = np.transpose(vects) @ numericalM @ vects
            for i in range(6):
                for j in range(6):
                    ## If on the diagonal compute abs % diff with eigenvalues and rotated matrix, if large then something went wrong
                    if i == j and abs((eigenValue[i]-diagonalBlock[i,i])/eigenValue[i]) > 1e-5:
                        print (f"Large difference in eigenValues at index {i},{j}")
                        print (f'The rotated mass matrix is {diagonalBlock}')
                    if i != j and diagonalBlock[i,j] > 1e-8:
                        print (f"Off diagonal element {i}{j} is larger than 1e-8, may not be diagonal")
                        print (f'The rotated mass matrix is {diagonalBlock}')
                        
            blockRot.append(vects)
            blockM.append(numericalM)

        # Note: @ is short for numpy matrix multiplication
        
        ## This diagonalizes the block-diagonal mass matrix
        blockDiagRot = linalg.block_diag(*blockRot)
        ## Diagonalized mass matrix
        diag = np.transpose(blockDiagRot) @ linalg.block_diag(*blockM) @ blockDiagRot

        ## Rotation that diagonalizes the original, unpermuted mass matrix.
        ## So we undo the permutation, and we need to transpose to match what we gave DRalgo
        rot = np.transpose(blockDiagRot) @ self.scalarPermutationMatrix

        ## OK we have the matrices that DRalgo used. But we now need to assign a correct value to each
        ## matrix element symbol in the Veff expressions. This is currently very hacky 
        
        outDict = self.scalarRotationMatrix.matchSymbols(rot)

        ## TODO improve this. currently I just hardcode scalar mass names
        massNames = ["MSsq01", "MSsq02", "MSsq03", "MSsq04", "MSsq05", "MSsq06", "MSsq07", "MSsq08", "MSsq09", "MSsq10", "MSsq11", "MSsq12"]

        if (self.bAbsoluteMsq):
            for i, msq in enumerate(np.diagonal(diag)):
                outDict[massNames[i]] = np.abs(msq)
        else:
            for i, msq in enumerate(np.diagonal(diag)):
                outDict[massNames[i]] = msq

        return outDict
    

""" Evaluating the potential: 
1. Call setModelParameters() with a dict that sets all parameters in the action. 
This is assumed to be using 3D EFT, so the params are temperature dependent.
2. Call evaluate() with a list that specifies values of background fields. Fields are in 3D units, ie. have dimension GeV^(1/2)
"""
class EffectivePotential:

    def __init__(self,
                 fieldNames, 
                 bAbsoluteMsq, 
                 vectorMassFile, 
                 vectorShorthandFile, 
                 scalarPermutationMatrix, 
                 scalerMassMatrices, 
                 scalarRotationMatrixFile,
                 loopOrder,
                 veffFiles):
        ## How many background fields do we depend on
        self.fieldNames = fieldNames
        self.nbrFields = len(self.fieldNames)

        self.params = VeffParams(fieldNames, 
                                 bAbsoluteMsq, 
                                 vectorMassFile, 
                                 vectorShorthandFile, 
                                 scalarPermutationMatrix, 
                                 scalerMassMatrices, 
                                 scalarRotationMatrixFile)
        
        self.loopOrder = loopOrder

        ## HACK: combine these into one file so that ParsedExpressionSystem understand it
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tempf:
            for filename in veffFiles:
                with open(filename, 'r') as f:
                    content = f.read()
                    tempf.write(content)
                    tempf.write("\n")

            ## close here because we need to re-open for parsing
            tempf.close()
            self.expressions = ParsedExpressionSystem(tempf.name)

        self.bNeedsDiagonalization = (self.loopOrder > 0)
        self.minimizer = VeffMinimizer(self.nbrFields) # currently the numVariables is not used by minimizer


    def setModelParameters(self, modelParameters: dict) -> None:
        """ This just reads action parameters from a dict and sets them internally for easier/faster(?) access in evaluate 
        """
        self.params.setActionParams(modelParameters)


    def initExpressions(self, filesToParse: list[str]) -> None:

        self.expressions = []


    def evaluate(self, fields: list[float]) -> complex:
        """Evaluate Veff at specified field values. Uses the currently set model parameters.
        """
        
        self.params.fields = fields
    
        ## This has masses, angles, all shorthand symbols etc. Everything we need to evaluate loop corrections
        paramDict = self.params.evaluateAll(fields, bNeedsDiagonalization=self.bNeedsDiagonalization)

        ## summing works because the result is a list [V0, V1, ...]
        res = sum( self.expressions.evaluateSystemWithDict(paramDict) )

        return res

    ## Return value is location, value
    def findLocalMinimum(self, T, initialGuess: list[float]) -> Tuple[list[float], complex]:
        
        ## I think we need to manually vectorize here if our parameters are arrays (ie. got multiple temperature inputs).
        ## Then self.evaluate would return a numpy array which scipy doesn't know how to work with. 
        ## Here I make an array of lambda functions and minimize those separately

        ## Minimize real part only:
        VeffWrapper = lambda fields: np.real ( self.evaluate(fields) )

        ##Added bounds to minimize to reduce the IR senstivity coming from low mass modes
        bounds = ((1e-6, 1e-6), (1e-6, 1e-6), (1e-6, 1e3))

        location, value = self.minimizer.minimize(T, VeffWrapper, initialGuess, bounds)


        if np.any(np.isnan(location)):
            location  = [np.nan] * len[initialGuess]
            location = np.asarray(location)
            value = np.nan
        else:
            ## evaluate once more to get the possible imag parts
            value = self.evaluate(location)

        return location, value
    

    def findGlobalMinimum(self, T, minimumCandidates: list[list[float]] = None) -> Tuple[list[float], complex]:
        """This calls findLocalMinimum with a bunch of initial guesses and figures out the deepest solution.
        Generally will not work very well if no candidates minima are given. 
        Return value is location, value. value can be complex (but this is probably a sign of failed minimization)
        """
        
        if (minimumCandidates == None or len(minimumCandidates) == 0):
            ## Initial guesses for minimiser
            minimumCandidates = [ [0., 0., 1e-4]]

        deepest = np.inf
        res = [np.nan]*3, np.inf

        ## Should we vectorize??
        for candidate in minimumCandidates:
            location, value = self.findLocalMinimum(T,candidate)
            if (np.real(value) < deepest):
                res = location, value
                deepest = np.real(value)

        return res


    # just calls self.evaluate
    def __call__(self, temperature: npt.ArrayLike, fields: npt.ArrayLike) -> complex:
        self.evaluate(temperature, fields)
