import ThreeHiggs

def getResourcePath(relativePathToResource):
    """ Gives a safe path to a packaged resource.
    
    Returns
    -------
    Path to the resource file: str.
    """


def getLines(relativePathToResource):
    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    from importlib.resources import files
    path = files(packageName) / relativePathToResource

    return open(path, encoding = "utf-8").readlines()

from ThreeHiggs.UserInput import UserInput
userinput = UserInput()
args = userinput.parse()

## ---- Configure Veff

hardToSoftFile = getLines(args.hardToSoftFile)
softScaleRGEFile = getLines(args.softScaleRGEFile)
softToUltrasoftFile = getLines(args.softToUltraSoftFile)

from ThreeHiggs.MathematicaParsers import parseExpressionSystem
from ThreeHiggs.ParsedExpression import ParsedExpressionSystem
vectorMassesSquared = ParsedExpressionSystem(parseExpressionSystem(getLines(args.vectorMassesSquaredFile)))
vectorShortHands = ParsedExpressionSystem(parseExpressionSystem(getLines(args.vectorShortHandsFile)))

from ThreeHiggs.MathematicaParsers import parseConstantMatrix
scalarPermutationMatrix = parseConstantMatrix(getLines(args.scalarPermutationFile))["matrix"]

from ThreeHiggs.MathematicaParsers import parseMassMatrix
from ThreeHiggs.ParsedExpression import MassMatrix
scalarMassMatrixUpperLeft = MassMatrix(parseMassMatrix(getLines(args.scalarMassMatrixUpperLeftFile))["matrix"],
                                       ParsedExpressionSystem(parseExpressionSystem(getLines(args.scalarMassMatrixUpperLeftDefinitionsFile))))

scalarMassMatrixBottomRight = MassMatrix(parseMassMatrix(getLines(args.scalarMassMatrixBottomRightFile))["matrix"],
                                         ParsedExpressionSystem(parseExpressionSystem(getLines(args.scalarMassMatrixBottomRightDefinitionsFile))))

scalarMassMatrices = [scalarMassMatrixUpperLeft, scalarMassMatrixBottomRight]

from ThreeHiggs.MathematicaParsers import parseRotationMatrix
from ThreeHiggs.ParsedExpression import RotationMatrix
scalarRotationMatrix = RotationMatrix(parseRotationMatrix(getLines(args.scalarRotationFile)))

veffLines = getLines(args.loFile)
if (args.loopOrder >= 1):
    veffLines += getLines(args.nloFile)
if (args.loopOrder >= 2):
    veffLines += getLines(args.nnloFile)

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
    
