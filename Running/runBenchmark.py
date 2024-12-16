def doMinimization(parameters):
    benchmark = parameters["benchmark"] if "benchmark" in parameters else None
    effectivePotential = parameters["effectivePotential"] if "effectivePotential" in parameters else None
    dimensionalReduction = parameters["dimensionalReduction"] if "dimensionalReduction" in parameters else None

    if args.bVerbose:
        print(f"Starting benchmark: {benchmark['bmNumber']}")

    if not args.firstBenchmark <= benchmark['bmNumber'] <= args.lastBenchmark:
        if args.bVerbose:
            print(f"Benchmark {benchmark['bmNumber']} has been rejected as outside benchmark range.")

        return

    from ThreeHiggs.TransitionFinder import traceFreeEnergyMinimum
    minimizationResult = traceFreeEnergyMinimum(effectivePotential, 
                                                dimensionalReduction, 
                                                benchmark,
                                                args.TRangeStart,
                                                args.TRangeEnd,
                                                args.TRangeStepSize,
                                                bVerbose = args.bVerbose)
  
    filename = f"{args.resultsDirectory}/BM_{benchmark['bmNumber']}"
    
    from pathlib import Path
    Path(args.resultsDirectory).mkdir(parents = True, exist_ok = True)
    from json import dumps
    if args.bSave:
        if args.bVerbose:
            print(f"Saving {benchmark['bmNumber']} to {filename}")
        open(f"{filename}.json", "w").write(dumps(minimizationResult, indent = 4))
      
    if args.bPlot:
        if args.bVerbose:
            print(f"Plotting {benchmark['bmNumber']}")

        from ThreeHiggs.PlotResult import plotData
        plotData(minimizationResult, benchmark['bmNumber'], args.loopOrder, filename)

    if args.bProcessMin:
        if args.bVerbose:
            print(f"Processing {benchmark['bmNumber']} to {filename+'_interp'}")
        from ThreeHiggs.ProcessMinimization import interpretData
        open(f"{filename}_interp.json", "w").write(dumps(interpretData(minimizationResult,
                                                                        benchmark["bmNumber"],
                                                                        benchmark["bmInput"]),
                                                         indent = 4))

def getLines(relativePathToResource):
    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    from importlib.resources import files
    path = files(packageName) / relativePathToResource

    return open(path, encoding = "utf-8").readlines()

def convertMathematica(args):
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
         open(args.parsedExpressionsFile, "w"),
         indent = 4)

def minimization(args):
    from json import load
    parsedExpressions = load(open(args.parsedExpressionsFile, "r"))

    from ThreeHiggs.EffectivePotential import EffectivePotential
    from ThreeHiggs.ParsedExpression import (ParsedExpressionSystem,
                                             MassMatrix,
                                             RotationMatrix)

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
                                            {"minAlgo" : args.minimizationAlgo, 
                                             "absGlobalTol" : args.absGlobalTolerance,\
                                             "relGlobalTol" :args.relGlobalTolerance, 
                                             "absLocalTol" : args.absLocalTolerance, 
                                             "relLocalTol" : args.relLocalTolerance,
                                             "v1Bounds" : args.v1Bounds,
                                             "v2Bounds" : args.v2Bounds,
                                             "v3Bounds" : args.v3Bounds},
                                            args.bNumba) 

    from ThreeHiggs.DimensionalReduction import DimensionalReduction
    dimensionalReduction = DimensionalReduction(ParsedExpressionSystem(parsedExpressions["hardToSoft"]),
                                                ParsedExpressionSystem(parsedExpressions["softScaleRGE"]),
                                                ParsedExpressionSystem(parsedExpressions["softToUltraSoft"]),
                                                bVerbose = args.bVerbose)
    
    from ThreeHiggs.BmGenerator import notRandomBmParam
    benchmarks = notRandomBmParam() if args.benchmarkMode == "handPicked" else []

    if args.bPool:
        from pathos.multiprocessing import Pool
        with Pool(args.cores) as pool:
            pool.map(doMinimization, ({"benchmark": item,
                                       "effectivePotential": effectivePotential,
                                       "dimensionalReduction": dimensionalReduction} for item in benchmarks))

    else:
        for parameters in ({"benchmark": item,
                            "effectivePotential": effectivePotential,
                            "dimensionalReduction": dimensionalReduction} for item in benchmarks):
            doMinimization(parameters)

from ThreeHiggs.UserInput import UserInput
userinput = UserInput()
args = userinput.parse()

from ThreeHiggs.UserInput import Stages
if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    convertMathematica(args)

if args.firstStage <= Stages.minimization <= args.lastStage:
   minimization(args)

