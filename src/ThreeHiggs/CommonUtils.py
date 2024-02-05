import importlib.resources
#import unicodedata

def getSafePathToResource(relativePathToResource: str) -> str:
    """ Gives a safe path to a packaged resource.
    
    Returns
    -------
    Path to the resource file: str.
    """

    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    return importlib.resources.files(packageName) / relativePathToResource


def replaceGreekSymbols(string: str) -> str:

    ## Unicode magic, this is definitely not ideal
    lowerCaseMu = u"\u03BC"
    lowerCaseLambda = u"\u03BB"

    newString = string 

    newString = newString.replace(lowerCaseLambda, "lam")
    newString = newString.replace(lowerCaseMu, "mu")

    
    """ TODO use unicodedata package here to do magic.
    """

    # NOTE: Manual replacements are definitely not a general solution. Consider: expression that contains both unicode lambda and separate symbol "lam" 
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
    """"""
    ## Combine dicts by unpacking both
    return {**dict1, **dict2}