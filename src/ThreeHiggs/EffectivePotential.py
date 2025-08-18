import numpy as np
from scipy import linalg
from numba import njit
from itertools import chain
import nlopt
from dataclasses import dataclass, InitVar

@njit
def diagonalizeNumba(matrices, matrixNumber, matrixSize, T):
    subEigenValues = np.empty( (matrixNumber, matrixSize) ) 
    subRotationMatrix = np.empty( (matrixNumber, matrixSize, matrixSize) )
    for idx, matrix in enumerate(matrices):
         subEigenValues[idx], subRotationMatrix[idx] = np.linalg.eigh(matrix)
    return subEigenValues*T**2, subRotationMatrix

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
    def __init__(self,
        fieldNames, 
        loopOrder,
        verbose,
        nloptInst,
        vectorMassesSquared, 
        vectorShortHands, 
        scalarPermutationMatrix,
        scalarMassMatrices, 
        scalarRotationMatrix,
        veffArray,
        allSymbols
    ):
        self.fieldNames = fieldNames
        
        self.loopOrder = loopOrder
        
        self.verbose = verbose
        
        self.nloptInst = nloptInst

        self.vectorMassesSquared = vectorMassesSquared
        self.vectorShortHands = vectorShortHands

        self.scalarPermutationMatrix = np.asarray(scalarPermutationMatrix, dtype = int)
        ## can have many matrices if we've block-diagonalized already
        ## ASSUME: the blocks are given in order: upper left to lower right. 
        ##TODO improve this
        self.scalarMassMatrices = [matrix for matrix in scalarMassMatrices]
        self.scalarRotationMatrix = scalarRotationMatrix
        self.expressionsArray = veffArray
        self.allSymbols = allSymbols

    def findGlobalMinimum(self,
          T, 
          params3D,
          minimumCandidates
    ):
        
        """For physics reasons we only minimise the real part,
        for nlopt reasons we need to give a redunant grad arg"""
        VeffWrapper = lambda fields, grad : np.real ( self.evaluatePotential(fields,
                                                                             T,
                                                                             params3D) )
        
        bestResult = self.nloptInst.nloptGlobal(VeffWrapper, 
                                                minimumCandidates[0])

        for candidate in minimumCandidates:
            result = self.nloptInst.nloptLocal(VeffWrapper,
                                               candidate)
            if result[1] < bestResult[1]:
                bestResult = result
                
        ## Potential computed again in case its complex
        return bestResult[0], self.evaluatePotential(bestResult[0],
                                                T, 
                                                params3D) 
    
    def evaluatePotential(self, 
          fields, 
          T, 
          params3D
    ):
        array = self.expressionsArray.dictToArray(self.computeMasses(
            fields,
            T,
            params3D
        ))

        return sum(self.expressionsArray.evaluateUnordered(array))
    
    
    def computeMasses(self,
        fields, 
        T, 
        params3D
    ) :
        for i, value in enumerate(fields):
            params3D[self.allSymbols.index(self.fieldNames[i])] = value

        params3D = self.vectorShortHands.evaluate(params3D)
        params3D = self.vectorMassesSquared.evaluate(params3D)
        ## Scalar masses doesn't take array (yet) so convert to dict
        params3D = {key: value for (key, value) in zip(self.allSymbols, params3D) }

        return self.diagonalizeScalars(
            params3D, 
            T
        )

    def diagonalizeScalars(self,
        params3D, 
        T
    ):
        """Finds a rotation matrix that diagonalizes the scalar mass matrix
        and returns a dict with diagonalization-specific params"""
        subMassMatrix = np.array( [matrix.evaluate(params3D) for matrix in self.scalarMassMatrices ]).real / T**2

        subEigenValues, subRotationMatrix = diagonalizeNumba(subMassMatrix, subMassMatrix.shape[0], subMassMatrix.shape[1], T)

        """ At the level of DRalgo we permuted the mass matrix to make it block diagonal, 
        so we need to undo the permutatation"""
        params3D |= self.scalarRotationMatrix.evaluate(self.scalarPermutationMatrix @ linalg.block_diag(*subRotationMatrix))
        ##TODO load names from mathematica
        massNames = ["MSsq01", "MSsq02", "MSsq03", "MSsq04", "MSsq05", "MSsq06", "MSsq07", "MSsq08", "MSsq09", "MSsq10", "MSsq11", "MSsq12"]
        return params3D | {name: float(msq) for name, msq in zip(massNames, chain(*subEigenValues))}
    
    ##Jasmine plotting tools
    def plotPot(self, T, params3D, linestyle, v3Min, potMin, v3Max):
        VeffWrapper = lambda fields, grad : np.real ( self.evaluatePotential(fields, T, params3D) )
        
        v3Range = np.linspace(1e-4, v3Max*1.1, 80)
        potArray = np.zeros(len(v3Range))
        for idx, v3 in enumerate(v3Range):
            potArray[idx] = VeffWrapper([1e-4, 1e-4, v3], 1)
        potArray = potArray/T**3
        import matplotlib.pylab as plt
        # plt.rcParams.update({'font.size': 12.5})
        plt.plot(v3Range, potArray-potArray[0], label = f"{T=:.1f} GeV", linestyle = linestyle)
        plt.scatter(v3Min, potMin/T**3 - potArray[0])
        return min(potArray-potArray[0]), max(potArray-potArray[0])
    
    def plotPot3D(self, T, params3D):
        import matplotlib.pyplot as plt
        import numpy as np

        from matplotlib import cm
        from matplotlib.ticker import LinearLocator
        
        VeffWrapper = lambda fields, grad : np.real ( self.evaluatePotential(fields, T, params3D) )
        n = 100
        v3Range = np.linspace(-1.5, 13, n)
        v2Range = np.linspace(-1.5, 5, n)
        v2Mesh, v3Mesh = np.meshgrid(v2Range, v3Range)
        potList = []
        for v3 in v3Range:
            for v2 in v2Range:
                potList.append(VeffWrapper([v2, v2, v3], None))
        minPot = min(potList)
        potArray = np.array(potList).reshape(n, n)
        potArray = (potArray-minPot)/T**3
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        zMax = 0.005
        ax.set_zlim(None, zMax)
        
        ax.plot_surface(v2Mesh, v3Mesh, potArray, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False, vmax = zMax)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.axes.set_xlabel(r"$v_1$, $v_2$ ($\text{GeV}^{\: \frac{1}{2}})$")
        ax.axes.set_ylabel(r"$v_3$ ($\text{GeV}^{\: \frac{1}{2}})$")
        ax.zaxis.set_rotate_label(False)
        ax.axes.set_zlabel(r"$\dfrac{\Delta V}{T^3}$", labelpad = 10, rotation = 0)
        # A StrMethodFormatter is used automatically
        ax.zaxis.set_major_formatter('{x:.03f}')

        plt.show()
        return None
