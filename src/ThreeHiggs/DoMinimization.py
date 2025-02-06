import json
from typing import Generator
import decimal

## This avoids floating point error in T gotten by np.arange or linspace
## However one must be careful as 1 = decimal.Decimal(1.000000000000001) 
def _drange(start: float, end: float, jump: str) -> Generator:
    start =  decimal.Decimal(start) 
    while start <= end:
        yield float(start)
        start += decimal.Decimal(jump)

def _doMinimization(parameters):
    ## This should be doable with **unpacking but difficult with pool (starmap?)
    benchmark = parameters["benchmark"] 
    effectivePotential = parameters["effectivePotential"] 
    dimensionalReduction = parameters["dimensionalReduction"] 
    pertSymbols = parameters["pertSymbols"] 
    argumentMap = parameters["argumentMap"]
    args = parameters["args"]
    if args.bVerbose:
        print(f"Starting benchmark: {benchmark['bmNumber']}")

    if not args.firstBenchmark <= benchmark['bmNumber'] <= args.lastBenchmark:
        if args.bVerbose:
            print(f"Benchmark {benchmark['bmNumber']} has been rejected as outside benchmark range.")

        return
    
    from ThreeHiggs.TransitionFinder import TraceFreeEnergyMinimum
    traceFreeEnergyMinimumInst = TraceFreeEnergyMinimum(config = {
        "effectivePotential":effectivePotential, 
        "dimensionalReduction": dimensionalReduction, 
        "TRange": tuple(_drange(args.TRangeStart, 
                                args.TRangeEnd, 
                                str(args.TRangeStepSize))),
        "pertSymbols": pertSymbols,
        "bVerbose": args.bVerbose,
        "initialGuesses": args.initialGuesses,
        "arg2Index": argumentMap,    
        "index2Arg": dict((v, k) for k, v in argumentMap.items())
    })
    
    
    minimizationResult = traceFreeEnergyMinimumInst.traceFreeEnergyMinimum(benchmark)
  
    filename = f"{args.resultsDirectory}/BM_{benchmark['bmNumber']}"
    
    from pathlib import Path
    Path(args.resultsDirectory).mkdir(parents = True, exist_ok = True)
    if args.bSave:
        if args.bVerbose:
            print(f"Saving {benchmark['bmNumber']} to {filename}")
        open(f"{filename}.json", "w").write(json.dumps(minimizationResult, indent = 4))
      
    if args.bPlot:
        if args.bVerbose:
            print(f"Plotting {benchmark['bmNumber']}")

        from ThreeHiggs.PlotResult import plotData
        plotData(minimizationResult, benchmark['bmNumber'], args.loopOrder, filename)

    if args.bProcessMin:
        if args.bVerbose:
            print(f"Processing {benchmark['bmNumber']} to {filename+'_interp'}")
        from ThreeHiggs.ProcessMinimization import interpretData
        open(f"{filename}_interp.json", "w").write(json.dumps(interpretData(minimizationResult,
                                                                        benchmark["bmNumber"],
                                                                        benchmark["bmInput"]),
                                                         indent = 4))
from ThreeHiggs.GetLines import getLines        
def minimization(args):
    parsedExpressions = json.load(open(args.parsedExpressionsFile, "r"))
    variableSymbols =  getLines( "Data/Variables/LagranianSymbols.json", mode = "json") 
    
    from ThreeHiggs.ParsedExpression import (ParsedExpressionSystem,
                                             MassMatrix,
                                             RotationMatrix)
    from ThreeHiggs.EffectivePotential import EffectivePotential, cNlopt
    nloptInst = cNlopt(config = {"nbrVars": len(variableSymbols["fieldSymbols"]), 
                                 "absGlobalTol" : args.absGlobalTolerance,
                                 "relGlobalTol" :args.relGlobalTolerance, 
                                 "absLocalTol" : args.absLocalTolerance, 
                                 "relLocalTol" : args.relLocalTolerance,
                                 "varLowerBounds" : args.varLowerBounds,
                                 "varUpperBounds" : args.varUpperBounds})
    
    effectivePotential = EffectivePotential(variableSymbols["fieldSymbols"],
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
        from ThreeHiggs.makeArgumentMap import makeArgumentMap
        minimizationDict = {"pertSymbols": frozenset(variableSymbols["fourPointSymbols"] + 
                                                     variableSymbols["yukawaSymbols"] + 
                                                     variableSymbols["gaugeSymbols"]), 
                            "effectivePotential": effectivePotential,
                            "dimensionalReduction": dimensionalReduction,
                            "argumentMap": makeArgumentMap(parsedExpressions),
                            "args": args} 
        if args.bPool:
            from pathos.multiprocessing import Pool
            with Pool(args.cores) as pool:
                from ijson import items
                pool.map(_doMinimization, (minimizationDict | {"benchmark": item} for item in items(benchmarkFile, "item", use_float = True)))

        else:
            for parameters in (minimizationDict | {"benchmark": item} for item in json.load(benchmarkFile)):
                _doMinimization(parameters)
