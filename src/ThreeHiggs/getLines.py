def getLines(relativePathToResource):
    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    from importlib.resources import files
    path = files(packageName) / relativePathToResource

    return open(path, encoding = "utf-8").readlines()

