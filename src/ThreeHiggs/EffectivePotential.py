import numpy as np
from scipy import linalg

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
        return linalg.eigh(matrix, check_finite = False)
    else:
        print(f"{method} is not assigned to a method in diagonalizeSymmetric, exiting program.")
        exit(-1)

def evaluateAll(fields: list[float], 
                T:float, 
                params3D, 
                fieldNames, 
                scalarPermutationMatrix,
                scalarMassMatrices,
                scalarRotationMatrix,
                diagAlgo,
                vectorShortHands,
                vectorMassesSquared,
                bAbsoluteMsq,
                bNeedsDiagonalization=True, 
                bVerbose = False) -> dict[str, float]:
    """This should return a dict that fixes all symbols needed for Veff 2-loop evaluation."""
    knownParamsDict = params3D.copy()

    ## Background fields
    for i, value in enumerate(fields):
        knownParamsDict[fieldNames[i]] = value

    ## Vectors
    knownParamsDict |= vectorShortHands(knownParamsDict, bReturnDict=True)
    vectorMasses = vectorMassesSquared(knownParamsDict, bReturnDict=True)

    for key, val in vectorMasses.items():
        vectorMasses[key] = np.abs(val) if bAbsoluteMsq else complex(val)

    knownParamsDict |= vectorMasses

    ## Scalars       
    knownParamsDict |= diagonalizeScalars(knownParamsDict, 
                                          T, 
                                          diagAlgo, 
                                          scalarPermutationMatrix,
                                          scalarMassMatrices,
                                          scalarRotationMatrix,
                                          bAbsoluteMsq,
                                          bVerbose)

    return knownParamsDict

def diagonalizeScalars(params: dict[str, float], 
                       T: float, 
                       diagAlgo, 
                       scalarPermutationMatrix,
                       scalarMassMatrices,
                       scalarRotationMatrix,
                       bAbsoluteMsq,
                       bVerbose = False) -> dict[str, float]:
    """Finds a rotation matrix that diagonalizes the scalar mass matrix
    and returns a dict with diagonalization-specific params"""
    # Diagonalize blocks separatey
    subRotationMatrix = []
    subEigenValues = []

    for matrix in scalarMassMatrices:
        numericalM = np.asarray(matrix(params))/T**2
        eigenValue, vects = diagonalizeSymmetric(numericalM, diagAlgo)
        eigenValue *=T**2
        ## NOTE: vects has the eigenvectors on columns => D = V^T . M . V, such that D is diagonal
        if bVerbose: ## 'Quick' check that the numerical mass matrix is within tol after being rotated by vects
            diagonalBlock = np.transpose(vects) @ numericalM @ vects
            offDiagonalIndex = np.where(~np.eye(diagonalBlock.shape[0],dtype=bool))
            if np.any(diagonalBlock[offDiagonalIndex] > 1e-8):
                print (f"Detected off diagonal element larger than 1e-8 tol,  'diagonal' mass matrix is: {diagonalBlock}")

        subEigenValues.append(eigenValue)                    
        subRotationMatrix.append(vects)
    fullRotationMatrix = linalg.block_diag(*subRotationMatrix)

    """ At the level of DRalgo we permuted the mass matrix to make it block diagonal, 
    we need to undo that permutation before we give the rotation matrix to the effectivate potential or something. 
    I am not 100% on this"""
    drAlgoRot = np.transpose(fullRotationMatrix) @ scalarPermutationMatrix

    ## OK we have the matrices that DRalgo used. But we now need to assign a correct value to each
    ## matrix element symbol in the Veff expressions. This is currently very hacky 
    outDict = scalarRotationMatrix(drAlgoRot)

    ##TODO this could be automated better if mass names were MSsq{i}, i.e. remove the 0 at the begining.
    ##But should probably be handled by a file given from mathematica
    massNames = ["MSsq01", "MSsq02", "MSsq03", "MSsq04", "MSsq05", "MSsq06", "MSsq07", "MSsq08", "MSsq09", "MSsq10", "MSsq11", "MSsq12"]
    from itertools import chain
    for i, msq in enumerate(tuple(chain(*subEigenValues))):
        outDict[massNames[i]] = abs(msq) if bAbsoluteMsq else complex(msq)

    return outDict 

import nlopt
def callNlopt(method: nlopt, 
              numVariables: int, 
              function: callable, 
              v1Bounds: tuple[float], 
              v2Bounds: tuple[float], 
              v3Bounds: tuple[float], 
              AbsTol: float, 
              relTol: float, 
              initialGuess: list[float]):
    
    opt = nlopt.opt(method, numVariables)
    functionWrapper = lambda fields, grad: function(fields) 
    opt.set_min_objective(functionWrapper)
    opt.set_lower_bounds((v1Bounds[0], v2Bounds[0], v3Bounds[0]))
    opt.set_upper_bounds((v1Bounds[1], v2Bounds[1], v3Bounds[1]))
    opt.set_xtol_abs(AbsTol)
    opt.set_xtol_rel(relTol)
    return opt.optimize(initialGuess),  opt.last_optimum_value()

def minimize(function: callable, 
             initialGuess: np.ndarray, 
             minimizationAlgo: str,
             numVariables: int, 
             globalAbs: float,
             globalRel: float,
             localAbs: float,
             localRel: float,
             v1Bounds: tuple[float],
             v2Bounds: tuple[float],
             v3Bounds: tuple[float]) -> tuple[np.ndarray, float]:
    """Even though we don't use the gradient, nlopt still tries to pass a grad arguemet to the function, so the function needs to be 
    wrapped to give it room for the grad arguement"""

    if minimizationAlgo == "scipy":
            import scipy.optimize
            bounds = ((v1Bounds[0], v1Bounds[1]), (v2Bounds[0], v2Bounds[1]), (v3Bounds[0], v3Bounds[1]))
            minimizationResult = scipy.optimize.minimize(function, initialGuess, bounds=bounds, tol = 1e-6)
            return minimizationResult.x, minimizationResult.fun
               
    elif minimizationAlgo == "directGlobal":
            location, _ = callNlopt(nlopt.GN_DIRECT_NOSCAL, 
                                        numVariables, 
                                        function, 
                                        v1Bounds, 
                                        v2Bounds, 
                                        v3Bounds, 
                                        globalAbs, 
                                        globalRel, 
                                        initialGuess)
                        
            return callNlopt(nlopt.LN_BOBYQA, 
                                            numVariables, 
                                            function, 
                                            v1Bounds, 
                                            v2Bounds, 
                                            v3Bounds, 
                                            localAbs, 
                                            localRel, 
                                            location)
    elif minimizationAlgo == "BOBYQA":            
            return callNlopt(nlopt.LN_BOBYQA, 
                                            numVariables, 
                                            function, 
                                            v1Bounds, 
                                            v2Bounds, 
                                            v3Bounds, 
                                            localAbs, 
                                            localRel, 
                                            initialGuess)
    
    else:
        print(f"ERROR: {minimizationAlgo} does not match any of our minimzationAlgos, attempting to exit")
        exit(-1)


    
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
                 vectorShortHands, 
                 scalarPermutationMatrix, 
                 scalarMassMatrices, 
                 scalarRotationMatrix,
                 loopOrder,
                 veff,
                 minimizationAlgo,
                 diagAlgo,
                 absGlobalTolerance,
                 relGlobalTolerance, 
                 absLocalTolerance, 
                 relLocalTolerance,
                 v1Bounds,
                 v2Bounds,
                 v3Bounds):
        self.fieldNames = fieldNames
        self.nbrFields = len(self.fieldNames)

        self.bAbsoluteMsq = bAbsoluteMsq
        self.diagAlgo = diagAlgo

        self.vectorMassesSquared = vectorMassesSquared
        self.vectorShortHands = vectorShortHands

        self.scalarPermutationMatrix = scalarPermutationMatrix
        ## can have many matrices if we've block-diagonalized already
        ## ASSUME: the blocks are given in order: upper left to lower right. 
        ##TODO improve this
        self.scalarMassMatrices = [matrix for matrix in scalarMassMatrices]

        self.scalarRotationMatrix = scalarRotationMatrix
        
        self.loopOrder = loopOrder
        self.minimizationAlgo = minimizationAlgo
        self.expressions = veff
        self.bNeedsDiagonalization = (self.loopOrder > 0)
        self.absGlobalTolerance = absGlobalTolerance
        self.relGlobalTolerance = relGlobalTolerance
        self.absLocalTolerance = absLocalTolerance
        self.relLocalTolerance = relLocalTolerance
        self.v1Bounds = v1Bounds
        self.v2Bounds = v2Bounds
        self.v3Bounds = v3Bounds

    def initExpressions(self, filesToParse: list[str]) -> None:
        self.expressions = []

    def evaluatePotential(self, fields: list[float], T:float, params3D, bVerbose = False) -> complex:
        ## This has masses, angles, all shorthand symbols etc. Everything we need to evaluate loop corrections
        ## Sum because the result is a list of tree, 1loop etc 
        return sum(self.expressions(evaluateAll(fields,
                                                T,
                                                params3D,
                                                self.fieldNames,
                                                self.scalarPermutationMatrix, 
                                                self.scalarMassMatrices, 
                                                self.scalarRotationMatrix,
                                                self.diagAlgo,
                                                self.vectorShortHands,
                                                self.vectorMassesSquared,
                                                self.bAbsoluteMsq,
                                                bNeedsDiagonalization=self.bNeedsDiagonalization, 
                                                bVerbose = bVerbose)))

    def findLocalMinimum(self, initialGuess: list[float],T:float, params3D, algo, bVerbose = False) -> tuple[list[float], complex]:
        ## Minimize real part only:
        VeffWrapper = lambda fields: np.real ( self.evaluatePotential(fields,
                                                                      T,
                                                                      params3D,
                                                                      bVerbose = bVerbose) )

        return minimize(VeffWrapper, 
                        initialGuess, 
                        algo,
                        self.nbrFields,
                        self.absGlobalTolerance,
                        self.relGlobalTolerance,
                        self.absLocalTolerance,
                        self.relLocalTolerance,
                        self.v1Bounds,
                        self.v2Bounds,
                        self.v3Bounds)

    def findGlobalMinimum(self,T:float, 
                          params3D,
                          minimumCandidates: list[list[float]] = None,
                          bVerbose = False) -> tuple[list[float], float, float, str]:
        bestResult = ((np.full(3, np.nan)), np.inf)
        
        if self.minimizationAlgo == "combo":
            result = self.findLocalMinimum(minimumCandidates[0], T, params3D, "directGlobal", bVerbose = bVerbose)
            if result[1] < bestResult[1]:
                bestResult = result
            
            for candidate in minimumCandidates:
                result = self.findLocalMinimum(candidate, T, params3D, "BOBYQA", bVerbose = bVerbose)
                if result[1] < bestResult[1]:
                    bestResult = result
                    
        else:
            for candidate in minimumCandidates:
                result = self.findLocalMinimum(candidate,T, params3D, self.minimizationAlgo, bVerbose = bVerbose)
                if result[1] < bestResult[1]:
                    bestResult = result
        
        if any(np.isnan(bestResult[0])) or np.isinf(bestResult[1]):
            return (np.full(3, None)), None, None, "NaN"
        
        potentialAtMin = self.evaluatePotential(bestResult[0], T, params3D, bVerbose = bVerbose) ## Compute the potential at minimum to check if its complex
        if abs(potentialAtMin.imag)/abs(potentialAtMin.real) > 1e-8: 
            return bestResult[0], potentialAtMin.real, potentialAtMin.imag, "complex" ## Flag minimum with imag > tol
        return bestResult[0], potentialAtMin.real, None, None
    
    
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

    def bReachedUltraSoftScale(self, fields: list[complex], T: float, params3D, bVerbose = False) -> bool:
        '''Check if we can trust the results by comparing the masses we find at the minimum to the ultra soft scale i.e.
        Are all physical masses > g^2 T/16pi, we use the largest coupling in the theory to do the comparrsion 
        --Note we expect some goldstone bosons from the symmetry breaking so we check the number of light modes = goldstone modes
        ----Get someone to check the logic of this
        2) Return true if # of light modes is less than the # of goldstone modes'''
        goldStone = 0 if np.all(np.abs(fields) < 0.1) else 3
        paramDict = evaluateAll(fields,
                                T,
                                params3D,
                                self.fieldNames,
                                self.scalarPermutationMatrix, 
                                self.scalarMassMatrices, 
                                self.scalarRotationMatrix,
                                self.diagAlgo,
                                self.vectorShortHands,
                                self.vectorMassesSquared,
                                self.bAbsoluteMsq,
                                bNeedsDiagonalization=self.bNeedsDiagonalization,
                                bVerbose = bVerbose)

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

