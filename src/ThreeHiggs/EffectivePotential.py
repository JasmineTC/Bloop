import nlopt
import numpy as np

from scipy import linalg
from itertools import chain
from dataclasses import dataclass, InitVar

from .Veff import Veff
from .Veff import Veff_params


def compFieldDepParams(
        fields: list[float], 
        T:float, 
        params3D, 
        fieldNames, 
        scalarPermutationMatrix,
        scalarMassMatrices,
        scalarRotationMatrix,
        vectorShortHands,
        vectorMassesSquared,
        bNumba,
        bVerbose
    ) -> dict[str, float]:
    ## Background fields
    for i, value in enumerate(fields):
        params3D[fieldNames[i]] = value

    ## Vectors
    params3D |= vectorShortHands.evaluate(params3D, bReturnDict=True)
    
    params3D |= vectorMassesSquared.evaluate(params3D, bReturnDict=True)

    ## Scalars       
    params3D |= diagonalizeScalars(
        params3D, 
        T,
        scalarPermutationMatrix,
        scalarMassMatrices,
        scalarRotationMatrix,
        bNumba,
        bVerbose
    )
    return params3D


def diagonalizeScalars(
        params: dict[str, float], 
        T: float,  
        scalarPermutationMatrix,
        scalarMassMatrices,
        scalarRotationMatrix,
        bNumba,
        bVerbose
    ) -> dict[str, float]:
    """Finds a rotation matrix that diagonalizes the scalar mass matrix
    and returns a dict with diagonalization-specific params"""
    subMassMatrix = np.array( [matrix.evaluate(params) for matrix in scalarMassMatrices ]).real / T**2
    if bNumba:
        from ThreeHiggs.diagonalizeNumba import diagonalizeNumba
        subEigenValues, subRotationMatrix = diagonalizeNumba(subMassMatrix, len(subMassMatrix), len(subMassMatrix[0][0]), T)
    else:
        subRotationMatrix = []
        subEigenValues = []
        for matrix in subMassMatrix:
            eigenValue, vects = np.linalg.eigh(matrix)
            eigenValue *=T**2
            
            subEigenValues.append(eigenValue)                    
            subRotationMatrix.append(vects)
            
            ## NOTE: vects has the eigenvectors on columns => D = V^T . M . V, such that D is diagonal
            ## 'Quick' check that the numerical mass matrix is within tol after being rotated by vects
            if not bVerbose:
                continue
            
            diagonalBlock = np.transpose(vects) @ matrix @ vects
            offDiagonalIndex = np.where(~np.eye(diagonalBlock.shape[0],dtype=bool))
            
            if not np.any(diagonalBlock[offDiagonalIndex] > 1e-8):
                continue
            print (f"Detected off diagonal element larger than 1e-8 tol,  'diagonal' mass matrix is: {diagonalBlock}")

    """ At the level of DRalgo we permuted the mass matrix to make it block diagonal, 
    so we need to undo the permutatation"""
    outDict = scalarRotationMatrix.evaluate(scalarPermutationMatrix @ linalg.block_diag(*subRotationMatrix))

    ##TODO this could be automated better if mass names were MSsq{i}, i.e. remove the 0 at the begining.
    ##But should probably be handled by a file given from mathematica (such a list is already made in mathematica)
    massNames = ["MSsq01", "MSsq02", "MSsq03", "MSsq04", "MSsq05", "MSsq06", "MSsq07", "MSsq08", "MSsq09", "MSsq10", "MSsq11", "MSsq12"]

    return outDict | {name: float(msq) for name, msq in zip(massNames, chain(*subEigenValues))}


@dataclass(frozen=True)
class cNlopt:
    nbrVars: int = 0
    varLowerBounds: tuple[float] = (0,) 
    varUpperBounds: tuple[float] = (0,) 
    absLocalTol: float = 0
    relLocalTol: float = 0
    absGlobalTol: float = 0
    relGlobalTol: float = 0
    config: InitVar[dict] = None
    
    ##Regular init method doesn't work with frozen data classes,
    ##Need to manually init by passing the class a dict i.e. class(config = dict)
    def __post_init__(self, config: dict):
        if config:
            self.__init__(**config)

    def nloptGlobal(self, func: callable, 
                 initialGuess: list[float]):
        opt = nlopt.opt(nlopt.GN_DIRECT_NOSCAL, self.nbrVars)
       	opt.set_min_objective(func)
       	opt.set_lower_bounds(self.varLowerBounds)
       	opt.set_upper_bounds(self.varUpperBounds)
        opt.set_xtol_abs(self.absGlobalTol)
       	opt.set_xtol_rel(self.relGlobalTol)
       	return self.nloptLocal(func, opt.optimize(initialGuess)) 
       
    def nloptLocal(self, func: callable, 
                 initialGuess: list[float]):
        opt = nlopt.opt(nlopt.LN_BOBYQA, self.nbrVars)
       	opt.set_min_objective(func)
       	opt.set_lower_bounds(self.varLowerBounds)
       	opt.set_upper_bounds(self.varUpperBounds)
       	opt.set_xtol_abs(self.absLocalTol) 
       	opt.set_xtol_rel(self.relLocalTol)
       	return opt.optimize(initialGuess),  opt.last_optimum_value()


""" Evaluating the potential: 
1. Call setModelParameters() with a dict that sets all parameters in the action. 
This is assumed to be using 3D EFT, so the params are temperature dependent.
2. Call evaluate() with a list that specifies values of background fields. Fields are in 3D units, ie. have dimension GeV^(1/2)
"""
class EffectivePotential:
    def __init__(
            self,
            fieldNames, 
            loopOrder,
            bNumba,
            bVerbose,
            nloptInst,
            vectorMassesSquared, 
            vectorShortHands, 
            scalarPermutationMatrix, 
            scalarMassMatrices, 
            scalarRotationMatrix,
            veff
        ):
        self.fieldNames = fieldNames
        self.loopOrder  = loopOrder
        self.bNumba     = bNumba
        self.bVerbose   = bVerbose
        self.nloptInst  = nloptInst
        
        self.vectorMassesSquared = vectorMassesSquared
        self.vectorShortHands = vectorShortHands

        self.scalarPermutationMatrix = scalarPermutationMatrix
        ## can have many matrices if we've block-diagonalized already
        ## ASSUME: the blocks are given in order: upper left to lower right. 
        ##TODO improve this
        self.scalarMassMatrices = [matrix for matrix in scalarMassMatrices]
        self.scalarRotationMatrix = scalarRotationMatrix
        self.expressions = veff

    def evaluatePotential(self, fields: list[float], T:float, params3D) -> complex:
        # This has masses, angles, all shorthand symbols etc. Everything we 
        # need to evaluate loop corrections Sum because the result is a list of 
        # tree, 1loop etc
        params = compFieldDepParams(
            fields,
            T,
            params3D,
            self.fieldNames,
            self.scalarPermutationMatrix, 
            self.scalarMassMatrices, 
            self.scalarRotationMatrix,
            self.vectorShortHands,
            self.vectorMassesSquared,
            self.bNumba,
            self.bVerbose
        )
        
        ps = Veff_params(params)
        return sum(Veff(*ps))
        
        # return sum(self.expressions.evaluate(params))


    def findGlobalMinimum(self,T:float, 
                          params3D,
                          minimumCandidates: list[list[float]] = None) -> tuple[list[float], float, float, str]:
        
        """For physics reasons we only minimise the real part,
        for nlopt reasons we need to give a redunant grad arg"""
        VeffWrapper = lambda fields, grad : np.real ( self.evaluatePotential(fields,
                                                                      T,
                                                                      params3D) )
        
        bestResult = self.nloptInst.nloptGlobal(VeffWrapper, minimumCandidates[0])
        
        for candidate in minimumCandidates:
            result = self.nloptInst.nloptLocal(VeffWrapper, candidate)
            if result[1] < bestResult[1]:
                bestResult = result
                    
        potentialAtMin = self.evaluatePotential(bestResult[0], T, params3D) ## Compute the potential at minimum to check if its complex
        if abs(potentialAtMin.imag)/abs(potentialAtMin.real) > 1e-8: 
            return bestResult[0], potentialAtMin.real, potentialAtMin.imag, "complex" ## Flag minimum with imag > tol
        return bestResult[0], potentialAtMin.real, None, None
    
    
    def getUltraSoftScale(self, paramDict, T: float) -> float:
        '''Given a field input and temperature compute the scale of ultra soft physics 
        which is g^2 T/ 16 pi, g is taken to be the largest coupling in the 
        theory to give the tightest constraint'''

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

    def bReachedUltraSoftScale(self, fields: list[complex], T: float, params3D) -> bool:
        '''Check if we can trust the results by comparing the masses we find at the minimum to the ultra soft scale i.e.
        Are all physical masses > g^2 T/16pi, we use the largest coupling in the theory to do the comparrsion 
        --Note we expect some goldstone bosons from the symmetry breaking so we check the number of light modes = goldstone modes
        ----Get someone to check the logic of this
        2) Return true if # of light modes is less than the # of goldstone modes'''
        goldStone = 0 if np.all(np.abs(fields) < 0.1) else 3
        paramDict = compFieldDepParams(fields,
                                T,
                                params3D,
                                self.fieldNames,
                                self.scalarPermutationMatrix, 
                                self.scalarMassMatrices, 
                                self.scalarRotationMatrix,
                                self.vectorShortHands,
                                self.vectorMassesSquared,
                                self.bNumba,
                                self.bVerbose)

        ## Convert mass into real type to do comparisons 
        massList = np.real([paramDict["MSsq01"], paramDict["MSsq02"],
                            paramDict["MSsq03"], paramDict["MSsq04"],
                            paramDict["MSsq05"], paramDict["MSsq06"],
                            paramDict["MSsq07"], paramDict["MSsq08"],
                            paramDict["MSsq09"], paramDict["MSsq10"],
                            paramDict["MSsq11"], paramDict["MSsq12"]])
    
        return len([lowMass for lowMass in massList if lowMass < self.getUltraSoftScale(paramDict, T)]) > goldStone

from unittest import TestCase
class EffectivePotentialUnitTests(TestCase):
    def test_diagonalizeNumba(self):
        reference = [[[-4., 4.], [-4.649110640673516, 20.64911064067352]], 
                     [[[-0.7071067811865475, 0.7071067811865475],                                                                                                                                        
                       [0.7071067811865475, 0.7071067811865475]],                                                                                                                                        
                      [[-0.9870874576374967, -0.1601822430069672],                                                                                                                                       
                       [-0.1601822430069672, 0.9870874576374967]]]  ]
        
        source = np.array( [ [[0, 1], [1, 0]], 
                             [[-1, 5], [-1, 5.0]] ] )
      
        from ThreeHiggs.diagonalizeNumba import diagonalizeNumba
        self.assertEqual(reference, 
                         list(map(lambda x: x.tolist(), 
                                  diagonalizeNumba(source, 2, 2, 2))))   
    
