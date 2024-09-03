def doMinimization(indexAndBenchMark):
    _, benchMark = indexAndBenchMark
    if not args.firstBenchmark <= benchMark['bmNumber'] <= args.lastBenchmark:
        return

    from ThreeHiggs.TransitionFinder import traceFreeEnergyMinimum
    minimizationResult = traceFreeEnergyMinimum(effectivePotential, 
                                                dimensionalReduction, 
                                                benchMark,
                                                args.TRangeStart,
                                                args.TRangeEnd,
                                                args.TRangeStepSize,
                                                verbose = args.verbose)
  
    filename = f"{args.resultsDirectory}/BM_{benchMark['bmNumber']}"
    
    from pathlib import Path
    Path(args.resultsDirectory).mkdir(parents = True, exist_ok = True)
    from json import dumps
    if args.bSave:
        open(f"{filename}.json", "w").write(dumps(minimizationResult, indent = 4))
      
    if args.bPlot:
        from PlotResult import PlotResult
        PlotResult.PlotData(minimizationResult, benchMark['bmNumber'], args.loopOrder, filename)

    if args.bProcessMin:
        from ThreeHiggs.ProcessMinimization import interpretData
        open(f"{filename}_interp.json", "w").write(dumps(interpretData(minimizationResult,
                                                                        benchMark["bmNumber"],
                                                                        benchMark["bmInput"]),
                                                         indent = 4))

def getLines(relativePathToResource):
    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    from importlib.resources import files
    path = files(packageName) / relativePathToResource

    return open(path, encoding = "utf-8").readlines()

from ThreeHiggs.UserInput import UserInput
userinput = UserInput()
args = userinput.parse()

from ThreeHiggs.UserInput import Stages
if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    with open(args.parsedExpressionsFile, "w") as parsedExpressionsFile:
        veffLines = getLines(args.loFile)
        if (args.loopOrder >= 1):
            veffLines += getLines(args.nloFile)
        if (args.loopOrder >= 2):
            veffLines += getLines(args.nnloFile)

        from ThreeHiggs.MathematicaParsers import (parseExpressionSystem,
                                                   parseConstantMatrix,
                                                   parseMassMatrix,
                                                   parseRotationMatrix)


        from json import dump
        dump({"vectorMassesSquared": parseExpressionSystem(getLines(args.vectorMassesSquaredFile)),
              "vectorShortHands": parseExpressionSystem(getLines(args.vectorShortHandsFile)),
              "scalarPermutationMatrix": parseConstantMatrix(getLines(args.scalarPermutationFile)),
              "scalarMassMatrixUpperLeft": parseMassMatrix(getLines(args.scalarMassMatrixUpperLeftDefinitionsFile),
                                                           getLines(args.scalarMassMatrixUpperLeftFile)),
              "scalarMassMatrixBottomRight": parseMassMatrix(getLines(args.scalarMassMatrixBottomRightDefinitionsFile),
                                                             getLines(args.scalarMassMatrixBottomRightFile)),
              "scalarRotationMatrix": parseRotationMatrix(getLines(args.scalarRotationFile)),
              "veff": parseExpressionSystem(veffLines),
              "hardToSoft": parseExpressionSystem(getLines(args.hardToSoftFile)),
              "softScaleRGE": parseExpressionSystem(getLines(args.softScaleRGEFile)),
              "softToUltraSoft": parseExpressionSystem(getLines(args.softToUltraSoftFile))},
             parsedExpressionsFile,
             indent = 4)

if args.firstStage <= Stages.minimization <= args.lastStage:
    with open(args.parsedExpressionsFile, "r") as parsedExpressionsFile:
        from json import load
        parsedExpressions = load(parsedExpressionsFile)

        from ThreeHiggs.EffectivePotential import EffectivePotential
        from ThreeHiggs.ParsedExpression import (ParsedExpressionSystem,
                                                 MassMatrix,
                                                 RotationMatrix)

        # TODO: This should not need to be treated as global in 3HDM.
        effectivePotential = EffectivePotential(['v1', 'v2', 'v3'],
                                                args.bAbsMass,
                                                ParsedExpressionSystem(parsedExpressions["vectorMassesSquared"]),
                                                ParsedExpressionSystem(parsedExpressions["vectorShortHands"]),
                                                parsedExpressions["scalarPermutationMatrix"]["matrix"],
                                                [MassMatrix(parsedExpressions["scalarMassMatrixUpperLeft"]), 
                                                 MassMatrix(parsedExpressions["scalarMassMatrixBottomRight"])],
                                                RotationMatrix(parsedExpressions["scalarRotationMatrix"]),
                                                args.loopOrder,
                                                ParsedExpressionSystem(parsedExpressions["veff"]),
                                                args.minimizationAlgo, ## Set algorithm to use for Veff minimization
                                                args.DiagAlgo, ## Set algorithm for scalar mass diag to use
                                                args.absGlobalTolerance,
                                                args.relGlobalTolerance, 
                                                args.absLocalTolerance, 
                                                args.relLocalTolerance,
                                                args.v1Bounds,
                                                args.v2Bounds,
                                                args.v3Bounds) 

    from ThreeHiggs.DimensionalReduction import DimensionalReduction
    dimensionalReduction = DimensionalReduction(ParsedExpressionSystem(parsedExpressions["hardToSoft"]),
                                                ParsedExpressionSystem(parsedExpressions["softScaleRGE"]),
                                                ParsedExpressionSystem(parsedExpressions["softToUltraSoft"]),
                                                verbose = args.verbose)
    
    with open(args.benchMarkFile) as benchMarkFile:
        if args.bPool:
            from multiprocessing import Pool
            with Pool(args.cores) as pool:
                from ijson import items
                pool.map(doMinimization, enumerate(items(benchMarkFile, "item", use_float = True)))
        else:
            from json import load
            benchMarkFile = load(benchMarkFile)
            for indexAndBenchMark in enumerate(benchMarkFile):
                doMinimization( indexAndBenchMark )
    

