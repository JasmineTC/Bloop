from enum import Enum
import nlopt
import numpy as np
#from numpy import *
import scipy.optimize
from typing import Callable, Tuple

##TODO, once global min is found, use as initial guess to find local min
##This is because global min methods focus on scanning parameter space and aren't as accurate as local min routines

class MinimizationAlgos(Enum):
    ##Enums work by setting the LHS, some important name you want to refer to later and keep fixed, to some unique number you don't care about
    eScipy = 1 ## direct scipy.optimize ~5.5mins to run bm1 at 2loop (does awfully)
    eBOBYQA = 2 ##An optimzation routine in nlopt ~ 1.5mins to run bm1 at 2loop
    eDIRECT = 3
    eNlopt = 999 ## some nlopt algo, TODO specify which one


class VeffMinimizer:

    numVariables: int # Not used currently
    __algo: MinimizationAlgos
    minimizer: Callable

    def __init__(self, numVariables: int):
        self.numVariables = numVariables
        self.__algo = MinimizationAlgos.eNlopt


    def setAlgorithm(self, algo: MinimizationAlgos) -> None:
        self.__algo = algo
    def setNumVariables(self, numVariables: int) -> None:
        self.numVariables = numVariables    
    def minimize(self, function: Callable, initialGuess: np.ndarray, bounds) -> Tuple[np.ndarray, float]:
    #def minimize(self, function: Callable, initialGuess: ndarray, bounds) -> Tuple[ndarray, float]:
        """Give bounds in format ((min1, max1), (min2, max2)) etc, one pair for each variable.
        Returns: 
        location, Veff(location)
        Note on notation for NLopt Algorithms, G/L refer to global/local and N/D refer to no gradient/gradient based, we want N methods
        Even though we don't use the gradient, nlopt still tries to pass a grad arguemet to the function, so the function needs to be 
        wrapped a second time to give it room for the useless grad arguement  
        """
        
        match(self.__algo):
            case MinimizationAlgos.eNlopt:
                opt = nlopt.opt(self.__algo, self.numVariables)


            case MinimizationAlgos.eScipy:
                minimizationResult = scipy.optimize.minimize(function, initialGuess, bounds=bounds, tol=1e-6)
                #print (f"number of interations = {minimizationResult.nit}")
                print (minimizationResult)
                
                location, value = minimizationResult.x, minimizationResult.fun
                
                
            case MinimizationAlgos.eDIRECT:
                ##The idea of this case is to use a global minimiser to get the ballpark of the global minimum
                ##then use that as initial guess for a local solver
                opt = nlopt.opt(nlopt.GN_DIRECT_NOSCAL, 3)
                ##Set function to minimise
                functionWrapper = lambda fields, grad: function(fields) 
                opt.set_min_objective(functionWrapper)
                ##Set lower bound variables on the minimisation varables to an array of length number of variables filled with 0
                #opt.set_lower_bounds(np.full(3, 1e-6))
                opt.set_lower_bounds((1e-6, 1e-6, 1e-6))
                ##Set upper bound variables on the minimisation varables  to an array of length number of variables filled with 100
                opt.set_upper_bounds((1e-6, 1e-6, 100))
                ##Set abs and rel tol on background field value
                opt.set_xtol_abs(0.1)
                opt.set_xtol_rel(0.1)
                location = opt.optimize(initialGuess)
                print (f"For an initial guess of {initialGuess} the global minimum is found to be {location}")
                
                opt2 = nlopt.opt(nlopt.LN_BOBYQA, 3)
                ##Set function to minimise
                opt2.set_min_objective(functionWrapper)
                ##Set lower bound variables on the minimisation varables to an array of length number of variables filled with 0
                opt2.set_lower_bounds(np.full(3, 1e-6))
                ##Set upper bound variables on the minimisation varables  to an array of length number of variables filled with 100
                opt2.set_upper_bounds((1e-6, 1e-6, 100))
                ##Set abs and rel tol on background field value
                opt2.set_xtol_abs(0.01)
                opt2.set_xtol_rel(0.001)
                
                location, value = opt2.optimize(location),  opt2.last_optimum_value()
                print(f"Given the above initial guess the local minimum is found to be {location}")
                
                
            case MinimizationAlgos.eBOBYQA:
                ##TODO use the bounds given to the function, rather than hard coding
                opt = nlopt.opt(nlopt.LN_BOBYQA, 3)
                ##Set function to minimise
                functionWrapper = lambda fields, grad: function(fields) 
                opt.set_min_objective(functionWrapper)
                ##Set lower bound variables on the minimisation varables to an array of length number of variables filled with 0
                opt.set_lower_bounds(np.full(3, 1e-6))
                ##Set upper bound variables on the minimisation varables  to an array of length number of variables filled with 100
                opt.set_upper_bounds((1e-6, 1e-6, 100))
                ##Set abs and rel tol on background field value
                opt.set_xtol_abs(0.0001)
                opt.set_xtol_rel(0.001)
                
                location, value = opt.optimize(initialGuess),  opt.last_optimum_value()
                print (location[2])
                
        return location, value
