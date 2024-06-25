from enum import Enum
import nlopt
import numpy as np

class MinimizationAlgos(Enum):
    eScipy = 1 
    eBOBYQA = 2 
    eDIRECTGLOBAL = 3 
    
class VeffMinimizer:

    def __init__(self, numVariables: int):
        self.numVariables = numVariables

    def setAlgorithm(self, algo: MinimizationAlgos) -> None:
        self.__algo = algo
        
    def setTolerances(self, globalAbs : float, globalRel : float, localAbs : float, localRel : float) -> None:
        self.globalAbs = globalAbs
        self.globalRel = globalRel
        
        self.localAbs = localAbs
        self.localRel = localRel
        
        
    def minimize(self, T, function: callable, initialGuess: np.ndarray, bounds) -> tuple[np.ndarray, float]:
        """Give bounds in format ((min1, max1), (min2, max2)) etc, one pair for each variable.
        Returns: 
        location, Veff(location)
        Note on notation for NLopt Algorithms, G/L refer to global/local and N/D refer to no gradient/gradient based, we want N methods
        Even though we don't use the gradient, nlopt still tries to pass a grad arguemet to the function, so the function needs to be 
        wrapped a second time to give it room for the useless grad arguement"""

        match(self.__algo):
            case MinimizationAlgos.eScipy:
                import scipy.optimize
                minimizationResult = scipy.optimize.minimize(function, initialGuess, bounds=bounds, tol = 1e-6)
                location, value = minimizationResult.x, minimizationResult.fun
                   
            case MinimizationAlgos.eDIRECTGLOBAL:
                ##The idea of this case is to use a global minimiser to get the ballpark of the global minimum
                ##then use that as initial guess for a local solver
                opt = nlopt.opt(nlopt.GN_DIRECT_NOSCAL, self.numVariables)
                functionWrapper = lambda fields, grad: function(fields) 
                opt.set_min_objective(functionWrapper)
                opt.set_lower_bounds((1e-6, 1e-6, 1e-6))
                opt.set_upper_bounds((1e-6, 1e-6, 100))
                opt.set_xtol_abs(self.globalAbs)
                opt.set_xtol_rel(self.globalRel)
                location = opt.optimize(initialGuess)
                
                opt2 = nlopt.opt(nlopt.LN_BOBYQA, self.numVariables)
                opt2.set_min_objective(functionWrapper)
                opt2.set_lower_bounds((1e-6, 1e-6, 1e-6))
                opt2.set_upper_bounds((1e-6, 1e-6, 100))
                opt2.set_xtol_abs(self.localAbs)
                opt2.set_xtol_rel(self.localRel)
                
                location, value = opt2.optimize(location),  opt2.last_optimum_value()
                
            case MinimizationAlgos.eBOBYQA:
                opt = nlopt.opt(nlopt.LN_BOBYQA, self.numVariables)
                functionWrapper = lambda fields, grad: function(fields) 
                opt.set_min_objective(functionWrapper)
                opt.set_lower_bounds((1e-6, 1e-6, 1e-6))
                opt.set_upper_bounds((1e-6, 1e-6, 100))
                opt.set_xtol_abs(self.localAbs)
                opt.set_xtol_rel(self.localRel)
                
                location, value = opt.optimize(initialGuess),  opt.last_optimum_value()
               
        return location, value
