import ThreeHiggs

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

hardToSoftFile = getResourcePath(args.hardToSoftFile)
softScaleRGEFile = getResourcePath(args.softScaleRGEFile)
softToUltrasoftFile = getResourcePath(args.softToUltraSoftFile)

from ThreeHiggs.MathematicaParsers import parseExpressionSystem
from ThreeHiggs.ParsedExpression import ParsedExpressionSystem
vectorMassesSquared = ParsedExpressionSystem(parseExpressionSystem(open(getResourcePath(args.vectorMassesSquaredFile), 
                                                                        encoding = "utf-8").readlines()))

vectorShortHands = ParsedExpressionSystem(parseExpressionSystem(open(getResourcePath(args.vectorShortHandsFile), 
                                                                     encoding = "utf-8").readlines()))

from ThreeHiggs.MathematicaParsers import parseConstantMatrix
scalarPermutationMatrix = parseConstantMatrix(open(getResourcePath(args.scalarPermutationFile), 
                                                   encoding = "utf-8").readlines())["matrix"]


from ThreeHiggs.MathematicaParsers import parseMassMatrix
from ThreeHiggs.ParsedExpression import MassMatrix
scalarMassMatrixUpperLeft = MassMatrix(parseMassMatrix(open(getResourcePath(args.scalarMassMatrixUpperLeftFile), encoding = "utf-8").readlines())["matrix"],
                                       ParsedExpressionSystem(parseExpressionSystem(open(getResourcePath(args.scalarMassMatrixUpperLeftDefinitionsFile),
                                                                                         encoding = "utf-8").readlines())))

scalarMassMatrixBottomRight = MassMatrix(parseMassMatrix(open(getResourcePath(args.scalarMassMatrixBottomRightFile), encoding = "utf-8").readlines())["matrix"],
                                       ParsedExpressionSystem(parseExpressionSystem(open(getResourcePath(args.scalarMassMatrixBottomRightDefinitionsFile),
                                                                                         encoding = "utf-8").readlines())))
scalarMassMatrices = [scalarMassMatrixUpperLeft, scalarMassMatrixBottomRight]

from ThreeHiggs.MathematicaParsers import parseRotationMatrix
from ThreeHiggs.ParsedExpression import RotationMatrix
scalarRotationMatrix = RotationMatrix(parseRotationMatrix(open(getResourcePath(args.scalarRotationFile), 
                                                               encoding = "utf-8").readlines())["matrix"])

veffLines = open(getResourcePath(args.loFile), encoding = "utf-8").readlines()
if (args.loopOrder >= 1):
    veffLines += open(getResourcePath(args.nloFile), encoding = "utf-8").readlines()
if (args.loopOrder >= 2):
    veffLines += open(getResourcePath(args.nnloFile), encoding = "utf-8").readlines()

veff = ParsedExpressionSystem(parseExpressionSystem(veffLines))

from ThreeHiggs.EffectivePotential import EffectivePotential
effectivePotential = EffectivePotential(['v1', 'v2', 'v3'],
                                        True,
                                        vectorMassesSquared,
                                        vectorShortHands,
                                        scalarPermutationMatrix,
                                        scalarMassMatrices,
                                        scalarRotationMatrix,
                                        args.loopOrder,
                                        veff,
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
    
