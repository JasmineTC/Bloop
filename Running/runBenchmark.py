import ThreeHiggs

def parseConstantMatrix(lines):
    symbols = [[symbol for symbol in line.rstrip()
                                         .rstrip('}')
                                         .lstrip('{')
                                         .split(',')] for line in lines]

    from sympy import Matrix
    sympyMatrix = Matrix(symbols)

    from numpy import array, float64
    return array(sympyMatrix.tolist()).astype(float64)

def getResourcePath(relativePathToResource: str) -> str:
    """ Gives a safe path to a packaged resource.
    
    Returns
    -------
    Path to the resource file: str.
    """

    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    from importlib.resources import files
    return files(packageName) / relativePathToResource

from ThreeHiggs.UserInput import UserInput
userinput = UserInput()
args = userinput.parse()

## ---- Configure Veff
veffFiles = [getResourcePath(args.loFile)]
if (args.loopOrder >= 1):
    veffFiles.append( getResourcePath(args.nloFile) )
if (args.loopOrder >= 2):
    veffFiles.append( getResourcePath(args.nnloFile) )

hardToSoftFile = getResourcePath(args.hardToSoftFile)
softScaleRGEFile = getResourcePath(args.softScaleRGEFile)
softToUltrasoftFile = getResourcePath(args.softToUltraSoftFile)

from ThreeHiggs.EffectivePotential import EffectivePotential
from ThreeHiggs.parsedmatrix import ParsedMatrix, MatrixDefinitionFiles
from ThreeHiggs.VeffMinimizer import MinimizationAlgos
effectivePotential = EffectivePotential(['v1', 'v2', 'v3'],
                                        True,
                                        getResourcePath(args.vectorsMassesSquaredFile),
                                        getResourcePath(args.vectorsShortHandsFile),
                                        parseConstantMatrix(open(getResourcePath(args.scalarPermutationFile), 
                                                            encoding = "utf-8").readlines()),
                                        [[getResourcePath(args.scalarMassMatrixUpperLeftFile),
                                          getResourcePath(args.scalarMassMatrixUpperLeftDefinitionsFile)],
                                         [getResourcePath(args.scalarMassMatrixBottomRightFile),
                                          getResourcePath(args.scalarMassMatrixBottomRightDefinitionsFile)]],
                                        getResourcePath(args.scalarRotationFile),
                                        args.loopOrder,
                                        veffFiles,
                                        args.minimizationAlgo, ## Set algorithm to use for Veff minimization
                                        args.DiagAlgo) ## Set algorithm for scalar mass diag to use

## Model object setup + load matching relations
from ThreeHiggs.GenericModel import GenericModel
model3HDM = GenericModel(effectivePotential)
model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile, softScaleRGEFile)
model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)

model3HDM.effectivePotential.minimizer.setTolerances(args.absGlobalTolerance,
                                                     args.relGlobalTolerance, 
                                                     args.absLocalTolerance, 
                                                     args.relLocalTolerance)

with open(args.benchMarkFile) as benchMarkFile:
    from json import load
    benchMarks = load(benchMarkFile)

    for index, benchMark in enumerate(benchMarks):
        if args.benchMarkNumber:
            if index != args.benchMarkNumber:
                continue
        
        from ThreeHiggs.TransitionFinder import TransitionFinder
        transitionFinder = TransitionFinder(model=model3HDM)
        model3HDM.setInputParams(benchMark)
        minimizationResults = transitionFinder.traceFreeEnergyMinimum(args.TRangeStart,
                                                                      args.TRangeEnd,
                                                                      args.TRangeStepSize)
        
        filename = f"{args.resultsDirectory}/BM_{args.benchMarkNumber}"
        
        from pathlib import Path
        Path(args.resultsDirectory).mkdir(parents = True, exist_ok = True)
        
        if args.save == True:
            from numpy import savetxt
            savetxt(filename + ".txt", minimizationResults)
        
        if args.plot == True:
            from PlotResult import PlotResult
            PlotResult.PlotData(minimizationResults, args.benchMarkNumber,args.loopOrder, filename)
    
