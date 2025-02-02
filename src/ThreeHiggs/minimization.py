def doMinimization(parameters):
    benchmark = parameters["benchmark"] if "benchmark" in parameters else None
    effectivePotential = parameters["effectivePotential"] if "effectivePotential" in parameters else None
    dimensionalReduction = parameters["dimensionalReduction"] if "dimensionalReduction" in parameters else None
    args = parameters["args"] if "args" in parameters else None

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

def minimization(args):
    from json import load
    parsedExpressions = load(open(args.parsedExpressionsFile, "r"))

    from ThreeHiggs.ParsedExpression import (ParsedExpressionSystem,
                                             MassMatrix,
                                             RotationMatrix)
    from ThreeHiggs.EffectivePotential import EffectivePotential, cNlopt
    nloptInst = cNlopt(config = {"nbrVars": 3, ##TODO this should be len(fieldNames) 
                                 "absGlobalTol" : args.absGlobalTolerance,
                                 "relGlobalTol" :args.relGlobalTolerance, 
                                 "absLocalTol" : args.absLocalTolerance, 
                                 "relLocalTol" : args.relLocalTolerance,
                                 "varLowerBounds" : args.varLowerBounds,
                                 "varUpperBounds" : args.varUpperBounds})
    effectivePotential = EffectivePotential(['v1', 'v2', 'v3'],
                                            args.loopOrder,
                                            args.bNumba,
                                            args.bVerbose,
                                            nloptInst,
                                            ParsedExpressionSystem(parsedExpressions["vectorMassesSquared"]),
                                            ParsedExpressionSystem(parsedExpressions["vectorShortHands"]),
                                            parsedExpressions["scalarPermutationMatrix"]["matrix"],
                                            [MassMatrix(parsedExpressions["scalarMassMatrixUpperLeft"]), 
                                             MassMatrix(parsedExpressions["scalarMassMatrixBottomRight"])],
                                            RotationMatrix(parsedExpressions["scalarRotationMatrix"]),
                                            ParsedExpressionSystem(parsedExpressions["veff"])) 

    from ThreeHiggs.DimensionalReduction import DimensionalReduction
    dimensionalReduction = DimensionalReduction(config = {"hardToSoft": ParsedExpressionSystem(parsedExpressions["hardToSoft"]),
                                                          "softScaleRGE": ParsedExpressionSystem(parsedExpressions["softScaleRGE"]),
                                                          "softToUltraSoft": ParsedExpressionSystem(parsedExpressions["softToUltraSoft"])})
    
    with open(args.benchmarkFile) as benchmarkFile:
        if args.bPool:
            from pathos.multiprocessing import Pool
            with Pool(args.cores) as pool:
                from ijson import items
                pool.map(doMinimization, ({"benchmark": item,
                                           "effectivePotential": effectivePotential,
                                           "dimensionalReduction": dimensionalReduction,
                                           "args": args} for item in items(benchmarkFile, "item", use_float = True)))

        else:
            for parameters in ({"benchmark": item,
                                "effectivePotential": effectivePotential,
                                "dimensionalReduction": dimensionalReduction,
                                "args": args} for item in load(benchmarkFile)):
                doMinimization(parameters)

