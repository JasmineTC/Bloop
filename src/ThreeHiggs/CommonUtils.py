import importlib.resources
import numpy as np
from typing import Tuple
from scipy import linalg
#import unicodedata

def replaceGreekSymbols(string: str) -> str:

    ## Unicode magic, this is definitely not ideal
    lowerCaseMu = u"\u03BC"
    lowerCaseLambda = u"\u03BB"

    newString = string 

    newString = newString.replace(lowerCaseLambda, "lam")
    newString = newString.replace(lowerCaseMu, "mu")

    
    """ TODO use unicodedata package here to do magic.
    """

    # NOTE: Manual replacements are definitely not a general solution. Consider problematic case: expression that contains both unicode lambda and separate symbol "lam" 
    # So tbh I'd like to keep the symbols are they are. But parse_mathematica from sympy does not seem to manage greek symbols at all!! 

    return newString


def dictToOrderedList(inDict: dict[any, any], keyOrder: list[any]):
    """ Picks values from a dict and puts them in a list with specified order.
    Currently this performs no checks on whether keys are found or not
    """
    outList = [None] * len(keyOrder)
    for i in range(len(outList)):
        outList[i] = inDict[ keyOrder[i] ]

    return outList


def combineDicts(dict1: dict[any, any], dict2: dict[any, any]) -> dict[any, any]:
    """Combine dicts by unpacking both"""
    ## Combine dicts by unpacking both
    return {**dict1, **dict2}


def diagonalizeSymmetric(matrix: np.ndarray, bCheckFinite: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Diagonalizes a symmetric matrix. Setting bCheckFinite = False may improve performance.
    Returns eigvalues, eigvectors in a matrix form
    """
    return linalg.eigh(matrix, check_finite=bCheckFinite)

