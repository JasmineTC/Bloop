import importlib.resources
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

def combineDicts(dict1: dict[any, any], dict2: dict[any, any]) -> dict[any, any]:
    """Combine dicts by unpacking both"""
    ## Combine dicts by unpacking both
    return {**dict1, **dict2}

