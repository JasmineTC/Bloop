import numpy as np
import tempfile # hacking

from scipy import linalg
from dataclasses import dataclass

from .parsedmatrix import ParsedMatrix, MatrixDefinitionFiles
from .ParsedExpression import ParsedExpressionSystem

from .VeffMinimizer import VeffMinimizer

def diagonalizeSymmetric(matrix: np.ndarray, method: str = "np") -> tuple[np.ndarray, np.ndarray]:
    """Diagonalizes a symmetric matrix. 
    Returns eigvalues, eigvectors in a matrix form
    For a general model the vector masses will also need this so don't factorise too hard!
    """
    if method == "np":
        return np.linalg.eigh(matrix)
    elif method == "mp":
        import mpmath as mp
        mp.mp.dps = 30
        eigenValue, eigenVector = mp.eigsy(mp.matrix(matrix), eigvals_only = False)
        return (np.array(eigenValue.tolist(),dtype=np.float64), np.array(eigenVector.tolist(),dtype=np.float64))
    elif method == "scipy":
        import scipy
        return scipy.linalg.eigh(matrix, check_finite = False)
    else:
        print(f"{method} is not assigned to a method in diagonalizeSymmetric, exiting program.")
        exit(-1)

@dataclass(slots=True)
class VeffConfig:
    # Names of the background-field variables that the Veff depends on
    fieldNames: list[str]
    # 0 = tree etc
    loopOrder: int
    # list of file names from where we parse expressions for the Veff. These get added together when evaluating Veff(fields)
    veffFiles: list[str]
    # Vector boson mass squares
    vectorMassFile: str
    # Shorthand symbols for vectors. Gets evaluated before masses
    vectorShorthandFile: str

    # Constant matrix transformation that brings the mass matrix into block diagonal form
    scalarPermutationMatrix: np.ndarray
    # Specify scalar mass matrices to diagonalize, can have many. The full mass matrix should be block diagonal in these (after the permutation transform)
    scalarMassMatrices: list[MatrixDefinitionFiles]
    # Full rotation from unpermutated basis to diagonal basis
    scalarRotationMatrixFile: str

    # Take absolute value of mass squares?
    ## TODO Didn't we set this to true in the config?
    bAbsoluteMsq: bool = False
    

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
        knownParamsDict |= self.vectorShorthands.evaluateSystemWithDict(knownParamsDict, bReturnDict=True)
        vectorMasses = self.vectorMassesSquared.evaluateSystemWithDict(knownParamsDict, bReturnDict=True)

        if (self.bAbsoluteMsq):
            for key, val in vectorMasses.items():
                vectorMasses[key] = np.abs(val)
        else:
            for key, val in vectorMasses.items():
                vectorMasses[key] = complex(val)

        knownParamsDict |= vectorMasses

        ## Scalars       
        knownParamsDict |= self.diagonalizeScalars(knownParamsDict)

        return knownParamsDict

    def diagonalizeScalars(self, params: dict[str, float]) -> dict[str, float]:
        """Finds a rotation matrix that diagonalizes the scalar mass matrix
        and returns a dict with diagonalization-specific params
        """
        # Diagonalize blocks separately
        subRotationMatrix = []
        subMassMatrix = []
        
        for matrix in self.scalarMassMatrices:
            numericalM = matrix.evaluateWithDict(params)
            eigenValue, vects = diagonalizeSymmetric( numericalM, "scipy")
            ## NOTE: vects has the eigenvectors on columns => D = V^T . M . V is diagonal
            verbose = False
            if verbose: ## 'Quick' check that the numerical mass matrix is within tol after being rotated by vects
                diagonalBlock = np.transpose(vects) @ numericalM @ vects
                offDiagonalIndex = np.where(~np.eye(diagonalBlock.shape[0],dtype=bool))
                if np.any(diagonalBlock[offDiagonalIndex]) > 1e-8:
                    print (f"Detected off diagonal element larger than 1e-8 tol,  'diagonal' mass matrix is: {diagonalBlock}")
   
                            
            subRotationMatrix.append(vects)
            subMassMatrix.append(numericalM)

        fullRotationMatrix = linalg.block_diag(*subRotationMatrix)

        fullMassMatrixDiag = np.transpose(fullRotationMatrix) @ linalg.block_diag(*subMassMatrix) @ fullRotationMatrix

        """ At the level of DRalgo we permuted the mass matrix to make it block diagonal, 
        we need to undo that permutation before we give the rotation matrix to the effectivate potential or something. 
        I am not 100% on this"""
        drAlgoRot = np.transpose(fullRotationMatrix) @ self.scalarPermutationMatrix

        ## OK we have the matrices that DRalgo used. But we now need to assign a correct value to each
        ## matrix element symbol in the Veff expressions. This is currently very hacky 
        outDict = self.scalarRotationMatrix.matchSymbols(drAlgoRot)

        massNames = ["MSsq01", "MSsq02", "MSsq03", "MSsq04", "MSsq05", "MSsq06", "MSsq07", "MSsq08", "MSsq09", "MSsq10", "MSsq11", "MSsq12"]

        if self.bAbsoluteMsq:
            for i, msq in enumerate(np.diagonal(fullMassMatrixDiag)):
                outDict[massNames[i]] = np.abs(msq)
        else:
            for i, msq in enumerate(np.diagonal(fullMassMatrixDiag)):
                outDict[massNames[i]] = complex(msq)

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
                 veffFiles,
                 minimizationAlgo):
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
        self.minimizationAlgo = minimizationAlgo

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


    def evaluatePotential(self, fields: list[float]) -> complex:
        """Evaluate Veff at specified field values. Uses the currently set model parameters.
        """
        
        self.params.fields = fields
    
        ## This has masses, angles, all shorthand symbols etc. Everything we need to evaluate loop corrections
        paramDict = self.params.evaluateAll(fields, bNeedsDiagonalization=self.bNeedsDiagonalization)

        return sum( self.expressions.evaluateSystemWithDict(paramDict) ) ## Sum because the result is a list of tree, 1loop etc 

    ## Return value is location, value
    def findLocalMinimum(self, T, initialGuess: list[float]) -> tuple[list[float], complex]:

        VeffWrapper = lambda fields: np.real ( self.evaluatePotential(fields) ) ## Minimize real part only:

        location, value = self.minimizer.minimize(VeffWrapper, initialGuess, self.minimizationAlgo)

        if np.any(np.isnan(location)):
            location = np.full(len[initialGuess], np.nan )
            value = np.nan
        else:
            value = self.evaluatePotential(location) ## evaluate once more to get the possible imag parts

        return location, value
    

    def findGlobalMinimum(self, T, minimumCandidates: list[list[float]] = None) -> tuple[list[float], complex]:
        """This calls findLocalMinimum with a bunch of initial guesses and figures out the deepest solution.
        Generally will not work very well if no candidates minima are given. 
        Return value is location, value. value can be complex (but this is probably a sign of failed minimization)
        """
        
        if not minimumCandidates:
            minimumCandidates = [ [1e-4, 1e-4, 1e-4]]

        deepest = np.inf
        res = [np.nan]*3, np.inf

        for candidate in minimumCandidates:
            location, value = self.findLocalMinimum(T,candidate)
            if (np.real(value) < deepest):
                res = location, value
                deepest = np.real(value)

        return res


    # just calls self.evaluate
    def __call__(self, temperature: np.ndarray, fields: np.ndarray) -> complex:
        self.evaluate(temperature, fields)
