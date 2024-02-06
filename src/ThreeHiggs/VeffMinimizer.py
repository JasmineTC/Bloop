import numpy as np
import scipy.optimize
from enum import Enum
from typing import Callable, Tuple

class MinimizationAlgos(Enum):
    eScipy = 1 ## direct scipy.optimize
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

    
    def minimize(self, function: Callable, initialGuess: np.ndarray, bounds) -> Tuple[np.ndarray, float]:
        """Give bounds in format ((min1, max1), (min2, max2)) etc, one pair for each variable.
        Returns: 
        location, Veff(location)
        """

        match(self.__algo):
            case MinimizationAlgos.eNlopt:

                ## TODO
                exit("NLOPT minimization is TODO")

            case MinimizationAlgos.eScipy:
                minimizationResult = scipy.optimize.minimize(function, initialGuess, bounds=bounds)
                location, value = minimizationResult.x, minimizationResult.fun
                
        
        return location, value



    