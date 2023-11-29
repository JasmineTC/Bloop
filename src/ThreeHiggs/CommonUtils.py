import importlib.resources


def getPackagedDataPath(relativeModulePath: str, fileName: str) -> str:
    """ Common routine for accessing packaged data files within a package, using modern importlib practices.
        Usage: if the file is <packageName>/Data/Something/example.txt, 
        call this as getPackagedDataPath("<packageName>.Data.Something", "example.txt").

        Returns
        -------
        Path to the resource file: str.
    """
    return str( importlib.resources.files(relativeModulePath).joinpath(fileName) )


def replaceGreekSymbols(string: str) -> str:

    ## Unicode magic, this is definitely not ideal

    lowerCaseMu = u"\u03BC"
    lowerCaseLambda = u"\u03BB"

    newString = string 

    newString = newString.replace(lowerCaseLambda, "lam")
    newString = newString.replace(lowerCaseMu, "mu")

    return newString
