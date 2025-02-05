from json import load
## Adding a mode is not ideal but idk what else to do
def getLines(relativePathToResource, mode = "default"):
    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    from importlib.resources import files
    path = files(packageName) / relativePathToResource
    
    if mode == "json":
        return load(open(path, "r"))
    return open(path, "r" , encoding = "utf-8").readlines()