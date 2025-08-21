from json import load
from importlib.resources import files

## Adding a mode is not ideal but idk what else to do
def getLines(relativePathToResource, mode="default"):
    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "Bloop"

    path = files(packageName) / relativePathToResource

    if mode == "json":
        with open(path, "r") as fp:
            return load(fp)
    with open(path, "r", encoding="utf-8") as fp:
        return fp.readlines()
