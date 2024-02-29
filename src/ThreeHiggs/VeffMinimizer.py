from enum import Enum
import nlopt
import numpy as np
import scipy.optimize
from typing import Callable, Tuple

class MinimizationAlgos(Enum):
    ##Enums work by setting the LHS, some important name you want to refer to later and keep fixed, to some unique number you don't care about
    eScipy = 1 ## direct scipy.optimize 
    ##~5.5mins to run bm1 at 2loop (does awfully)
    eBOBYQA = 2 ##An optimzation routine in nlopt 
    ##~ 1.5mins to run bm1 at 2loop - senstive to initial conditions
    eDIRECTGLOBAL = 3 ##An optimzation routine in nlopt 
    ##~ 5mins at 2loop - (seems to be) indepedent of initital guess as it uses global then local search
    eNelderMead = 4 
    ##~1m20s - dependent on initial condition
    eSbplx = 5
    ##Seems to be bad don't bother
class VeffMinimizer:

    numVariables: int # Not used currently
    __algo: MinimizationAlgos
    minimizer: Callable

    def __init__(self, numVariables: int):
        self.numVariables = numVariables

    

    def setAlgorithm(self, algo: MinimizationAlgos) -> None:
        self.__algo = algo
        
    def setTemp(self, temp: float) -> None:
        self.temp = temp
        print ("set temp called")
        
    def setgHDM(self, ghdm: float) -> None:
        self.ghdm = ghdm
        
    def setNumVariables(self, numVariables: int) -> None:
        self.numVariables = numVariables
        
    def setTolerances(self, globalAbs : float, globalRel : float, localAbs : float, localRel : float) -> None:
        self.globalAbs = globalAbs
        self.globalRel = globalRel
        
        self.localAbs = localAbs
        self.localRel = localRel
        
    def setBmNumber(self, bmNumber : int) -> None:
        self.bmNumber = bmNumber
        
    def minimize(self, function: Callable, initialGuess: np.ndarray, bounds) -> Tuple[np.ndarray, float]:
        """Give bounds in format ((min1, max1), (min2, max2)) etc, one pair for each variable.
        Returns: 
        location, Veff(location)
        Note on notation for NLopt Algorithms, G/L refer to global/local and N/D refer to no gradient/gradient based, we want N methods
        Even though we don't use the gradient, nlopt still tries to pass a grad arguemet to the function, so the function needs to be 
        wrapped a second time to give it room for the useless grad arguement  
        """
        ##TODO take bounds as user input, rather than hard coding

        match(self.__algo):

            case MinimizationAlgos.eScipy:
                minimizationResult = scipy.optimize.minimize(function, initialGuess, bounds=bounds, tol = 1e-6)
                #print (f"number of interations = {minimizationResult.nit}")
                print (minimizationResult)
                
                location, value = minimizationResult.x, minimizationResult.fun
                
                
            case MinimizationAlgos.eDIRECTGLOBAL:
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
                opt.set_xtol_abs(self.globalAbs)
                opt.set_xtol_rel(self.globalRel)
                location = opt.optimize(initialGuess)
                #print (f"For an initial guess of {initialGuess} the global minimum is found to be {location}")
                
                opt2 = nlopt.opt(nlopt.LN_BOBYQA, 3)
                ##Set function to minimise
                opt2.set_min_objective(functionWrapper)
                ##Set lower bound variables on the minimisation varables to an array of length number of variables filled with 0
                opt2.set_lower_bounds(np.full(3, 1e-6))
                ##Set upper bound variables on the minimisation varables  to an array of length number of variables filled with 100
                opt2.set_upper_bounds((1e-6, 1e-6, 100))
                ##Set abs and rel tol on background field value
                opt2.set_xtol_abs(self.localRel)
                opt2.set_xtol_rel(self.localRel)
                
                location, value = opt2.optimize(location),  opt2.last_optimum_value()
                ##For testing how well minimiser does
                # points = 150
                # xMax = max(2, location[2]) 
                # xList = np.linspace(1e-4, xMax*1.4, points)
                # yList = np.zeros(points)
                # for i, value in enumerate(xList):
                #     yList[i] = function( (0, 0, value) )
                # yList = yList - yList[0]
                # plt.plot(xList, yList, '.')
                # plt.xlabel('v3')
                # plt.ylabel('V')
                # plt.vlines(location[2], min(yList), max(yList))
                # plt.title(f'Benchmark {self.bmNumber} - gHDM {self.ghdm}, T {T} 1loop')
                # plt.savefig(f"Results/Debug/g_01/1loop/BM_{self.bmNumber}_gHDM_{self.ghdm}_T_{T}_1loop.png")
                # plt.close()
                
            case MinimizationAlgos.eNelderMead:
                opt = nlopt.opt(nlopt.LN_NELDERMEAD, 3)
                ##Set function to minimise
                functionWrapper = lambda fields, grad: function(fields) 
                opt.set_min_objective(functionWrapper)
                ##Set lower bound variables on the minimisation varables to an array of length number of variables filled with 0
                opt.set_lower_bounds(np.full(3, 1e-6))
                ##Set upper bound variables on the minimisation varables  to an array of length number of variables filled with 100
                opt.set_upper_bounds((1e-6, 1e-6, 100))
                ##Set abs and rel tol on background field value
                opt.set_xtol_abs(self.localAbs)
                opt.set_xtol_rel(self.localRel)
                
                location, value = opt.optimize(initialGuess),  opt.last_optimum_value()
                print (f"For an initial guess of {initialGuess} the local minimum is found to be {location}")
            case MinimizationAlgos.eSbplx:
                    opt = nlopt.opt(nlopt.LN_SBPLX, 3)
                    ##Set function to minimise
                    functionWrapper = lambda fields, grad: function(fields) 
                    opt.set_min_objective(functionWrapper)
                    ##Set lower bound variables on the minimisation varables to an array of length number of variables filled with 0
                    opt.set_lower_bounds(np.full(3, 1e-6))
                    ##Set upper bound variables on the minimisation varables  to an array of length number of variables filled with 100
                    opt.set_upper_bounds((1e-6, 1e-6, 100))
                    ##Set abs and rel tol on background field value
                    opt.set_xtol_abs(self.localAbs)
                    opt.set_xtol_rel(self.localRel)
                    
                    location, value = opt.optimize(initialGuess),  opt.last_optimum_value()
                    print (f"For an initial guess of {initialGuess} the local minimum is found to be {location}")   
                

                
        return location, value
