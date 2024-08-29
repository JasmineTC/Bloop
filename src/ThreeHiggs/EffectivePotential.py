import numpy as np
from scipy import linalg
from dataclasses import dataclass

from .VeffMinimizer import VeffMinimizer

def diagonalizeSymmetric(matrix: np.ndarray, method: str = "np") -> tuple[np.ndarray, np.ndarray]:
    """Diagonalizes a symmetric matrix. 
    Returns eigvalues, eigvectors in a matrix form
    For a general model the vector masses will also need this so don't factorise too hard!
    """
    if method == "np" or method == "numpy":
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
                 vectorMassesSquared, 
                 vectorShortHands, 
                 scalarPermutationMatrix, 
                 scalarMassMatrices, 
                 scalarRotationMatrix,
                 diagonalizationAlgo):
        
        self.fieldNames = fieldNames
        self.bAbsoluteMsq = bAbsoluteMsq
        self.diagonalizationAlgo = diagonalizationAlgo

        self.vectorMassesSquared = vectorMassesSquared
        self.vectorShortHands = vectorShortHands

        self.scalarPermutationMatrix = scalarPermutationMatrix
        ## can have many matrices if we've block-diagonalized already
        ## ASSUME: the blocks are given in order: upper left to lower right. 
        ##TODO improve this
        self.scalarMassMatrices = [matrix for matrix in scalarMassMatrices]

        self.scalarRotationMatrix = scalarRotationMatrix

    def setActionParams(self, inputParams: dict[str, float]) -> None:
        self.actionParams = inputParams


    def evaluateAll(self, fields: list[float], bNeedsDiagonalization=True, verbose = False) -> dict[str, float]:
        """This should return a dict that fixes all symbols needed for Veff 2-loop evaluation.
        """
        # Gradually build a dict containing (key, val) for all needed symbols
        # TODO this will need to be optimized
        knownParamsDict = self.actionParams.copy()

        ## Background fields
        for i, value in enumerate(fields):
            knownParamsDict[self.fieldNames[i]] = value

        ## Vectors
        knownParamsDict |= self.vectorShortHands(knownParamsDict, bReturnDict=True)
        vectorMasses = self.vectorMassesSquared(knownParamsDict, bReturnDict=True)

        if (self.bAbsoluteMsq):
            for key, val in vectorMasses.items():
                vectorMasses[key] = np.abs(val)
        else:
            for key, val in vectorMasses.items():
                vectorMasses[key] = complex(val)

        knownParamsDict |= vectorMasses

        ## Scalars       
        knownParamsDict |= self.diagonalizeScalars(knownParamsDict, verbose)

        return knownParamsDict

    def diagonalizeScalars(self, params: dict[str, float], verbose = False) -> dict[str, float]:
        """Finds a rotation matrix that diagonalizes the scalar mass matrix
        and returns a dict with diagonalization-specific params
        """
        # Diagonalize blocks separately
        subRotationMatrix = []
        subMassMatrix = []

        for matrix in self.scalarMassMatrices:
            numericalM = matrix(params)
            eigenValue, vects = diagonalizeSymmetric( numericalM, self.diagonalizationAlgo)

            ## NOTE: vects has the eigenvectors on columns => D = V^T . M . V is diagonal
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

        outDict = self.scalarRotationMatrix(drAlgoRot)

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
                 vectorMassesSquared, 
                 vectorShorthands, 
                 scalarPermutationMatrix, 
                 scalerMassMatrices, 
                 scalarRotationMatrix,
                 loopOrder,
                 veff,
                 minimizationAlgo,
                 diagonalizationAlgo,
                 absGlobalTolerance,
                 relGlobalTolerance, 
                 absLocalTolerance, 
                 relLocalTolerance,
                 v1Bounds,
                 v2Bounds,
                 v3Bounds):
        ## How many background fields do we depend on
        self.fieldNames = fieldNames
        self.nbrFields = len(self.fieldNames)

        self.params = VeffParams(fieldNames, 
                                 bAbsoluteMsq, 
                                 vectorMassesSquared, 
                                 vectorShorthands, 
                                 scalarPermutationMatrix, 
                                 scalerMassMatrices, 
                                 scalarRotationMatrix,
                                 diagonalizationAlgo)
        
        self.loopOrder = loopOrder
        self.minimizationAlgo = minimizationAlgo
        self.expressions = veff
        self.bNeedsDiagonalization = (self.loopOrder > 0)
        self.minimizer = VeffMinimizer(self.nbrFields,                 
                                       absGlobalTolerance,
                                       relGlobalTolerance, 
                                       absLocalTolerance, 
                                       relLocalTolerance,
                                       v1Bounds,
                                       v2Bounds,
                                       v3Bounds)

    def setModelParameters(self, modelParameters: dict) -> None:
        """ This just reads action parameters from a dict and sets them internally for easier/faster(?) access in evaluate 
        """
        self.params.setActionParams(modelParameters)


    def initExpressions(self, filesToParse: list[str]) -> None:

        self.expressions = []


    def evaluatePotential(self, fields: list[float], verbose = False) -> complex:
        """Evaluate Veff at specified field values. Uses the currently set model parameters.
        """
        
        self.params.fields = fields
    
        ## This has masses, angles, all shorthand symbols etc. Everything we need to evaluate loop corrections
        paramDict = self.params.evaluateAll(fields, 
                                            bNeedsDiagonalization=self.bNeedsDiagonalization, 
                                            verbose = verbose)

        return sum(self.expressions(paramDict)) ## Sum because the result is a list of tree, 1loop etc 

    def findLocalMinimum(self, initialGuess: list[float], algo, verbose = False) -> tuple[list[float], complex]:
        
        VeffWrapper = lambda fields: np.real ( self.evaluatePotential(fields, 
                                                                      verbose = verbose) ) ## Minimize real part only:

        location, value = self.minimizer.minimize(VeffWrapper, initialGuess, algo)

        value = self.evaluatePotential(location, verbose = verbose) ## evaluate once more to get the possible imag parts
        if abs(value.imag)/abs(value.real) > 1e-8: ## If the imaginary part of the potential is large then minimisation has found unphysical
            value = np.inf ##Set the value to inf so it won't be ever be smaller than best result
        return location, value.real ## Taking the real part so compactiable with json

    def findGlobalMinimum(self, minimumCandidates: list[list[float]] = None, verbose = False) -> tuple[list[float], complex]:
        """This calls findLocalMinimum with a bunch of initial guesses and figures out the deepest solution.
        Return value is location, value. value can be complex (but this is probably a sign of failed minimization)
        """
        
        if not minimumCandidates:
            minimumCandidates = [ [1e-4, 1e-4, 1e-4] ]

        bestResult = ((np.full(3, np.nan)), np.inf)
        
        if self.minimizationAlgo == "combo":
            result = self.findLocalMinimum(minimumCandidates[0], "directGlobal", verbose = verbose)
            if result[1] < bestResult[1]:
                bestResult = result
            
            for candidate in minimumCandidates:
                result = self.findLocalMinimum(candidate, "BOBYQA", verbose = verbose)
                if result[1] < bestResult[1]:
                    bestResult = result
                    
        else:
            for candidate in minimumCandidates:
                result = self.findLocalMinimum(candidate, self.minimizationAlgo, verbose = verbose)
                if result[1] < bestResult[1]:
                    bestResult = result
                    
        if any(np.isnan(bestResult[0])) or np.isinf(bestResult[1]):
            return ((np.full(3, None)), None)
        else:
            return bestResult
    
    
    def getUltraSoftScale(self, paramDict, T: float) -> float:
        '''Returns
        -------
        float
            Given a field input and temperature compute the scale of ultra soft physics which is g^2 T/ 16 pi
            g is taken to be the largest coupling in the theory to give the tightest constraint'''

        couplings = np.array([paramDict['lam1Re'], paramDict['lam1Im'],
                              paramDict['lam11'], paramDict['lam22'],
                              paramDict['lam12'], paramDict['lam12p'],
                              paramDict['g1'], paramDict['g2'],
                              paramDict['g3'], paramDict['lam33'],
                              paramDict['lam2Re'], paramDict['lam2Im'],
                              paramDict['lam23'], paramDict['lam23p'],
                              paramDict['lam3Re'], paramDict['lam3Im'],
                              paramDict['lam31'], paramDict['lam31p']]) 
        return np.max(np.abs(couplings))**2 * T / (16*np.pi)

    def bReachedUltraSoftScale(self, fields: list[complex], T: float, verbose = False) -> bool:
        '''
        Check if we can trust the results by comparing the masses we find at the minimum to the ultra soft scale i.e.
        Are all physical masses > g^2 T/16pi, we use the largest coupling in the theory to do the comparrsion 
        --Note we expect some goldstone bosons from the symmetry breaking so we check the number of light modes = goldstone modes
        ----Get someone to check the logic of this
        2) Return true if # of light modes is less than the # of goldstone modes'''
        ## If all field values are close to zero then we are probably in the symmetric phase with no goldstone bosons
        ## TODO: the symmetric phase cannot be described by pert theory so shouldn't be trusted, so should prob just return False right away
        ## But get someone to check this
        goldStone = 0 if np.all(np.abs(fields) < 0.1) else 3
        paramDict = self.params.evaluateAll(fields, 
                                            bNeedsDiagonalization=self.bNeedsDiagonalization,
                                            verbose = verbose)

        ultraSoftScale = self.getUltraSoftScale(paramDict, T)
    
        ## Convert mass into real type to do comparisons 
        massList = np.real([paramDict["MSsq01"], paramDict["MSsq02"],
                            paramDict["MSsq03"], paramDict["MSsq04"],
                            paramDict["MSsq05"], paramDict["MSsq06"],
                            paramDict["MSsq07"], paramDict["MSsq08"],
                            paramDict["MSsq09"], paramDict["MSsq10"],
                            paramDict["MSsq11"], paramDict["MSsq12"]])
    
        return len([lowMass for lowMass in massList if lowMass < ultraSoftScale]) > goldStone

    # just calls self.evaluate
    def __call__(self, temperature: np.ndarray, fields: np.ndarray) -> complex:
        self.evaluate(temperature, fields)

from unittest import TestCase
class EffectivePotentialUnitTests(TestCase):
    def test_diagonalizeSymmetricNumpy(self):
        reference = [[-1, 1], [[-0.7071067811865475, 0.7071067811865475], 
                               [0.7071067811865475, 0.7071067811865475]]]

        source = [[0, 1], [1, 0]]

        self.assertEqual(reference, 
                         list(map(lambda x: x.tolist(), 
                                  diagonalizeSymmetric(source, "numpy"))))

    def test_diagonalizeSymmetricMp(self):
        reference = [[[-1], [1]], [[0.7071067811865476, 0.7071067811865476], 
                                   [-0.7071067811865476, 0.7071067811865476]]]

        source = [[0, 1], [1, 0]]

        self.assertEqual(reference, 
                         list(map(lambda x: x.tolist(), 
                                  diagonalizeSymmetric(source, "mp"))))

    def test_diagonalizeSymmetricScipy(self):
        reference = [[-1, 1], [[-0.7071067811865475, 0.7071067811865475], 
                               [0.7071067811865475, 0.7071067811865475]]]

        source = [[0, 1], [1, 0]]

        self.assertEqual(reference, 
                         list(map(lambda x: x.tolist(), 
                                  diagonalizeSymmetric(source, "scipy"))))

