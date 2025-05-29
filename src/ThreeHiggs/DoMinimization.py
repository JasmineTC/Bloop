import json
from typing import Generator
import decimal

## This (sometimes) avoids floating point error in T gotten by np.arange or linspace
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
    betaFunction4DExpression = parameters["betaFunction4DExpression"]
    dimensionalReduction = parameters["dimensionalReduction"] 
    pertSymbols = parameters["pertSymbols"] 
    args = parameters["args"]
    allSymbols = parameters["allSymbols"]
    

    if args.verbose:
        print(f"Starting benchmark: {benchmark['bmNumber']}")

    if not args.firstBenchmark <= benchmark['bmNumber'] <= args.lastBenchmark:
        if args.verbose:
            print(f"Benchmark {benchmark['bmNumber']} has been rejected as outside benchmark range.")

        return
    from ThreeHiggs.TransitionFinder import TraceFreeEnergyMinimum
    traceFreeEnergyMinimumInst = TraceFreeEnergyMinimum(config = {"effectivePotential":effectivePotential, 
                                                                  "dimensionalReduction": dimensionalReduction, 
                                                                  "betaFunction4DExpression": betaFunction4DExpression,
                                                                  "TRange": tuple(_drange(args.TRangeStart, 
                                                                       args.TRangeEnd, 
                                                                       str(args.TRangeStepSize))),
                                                                  "pertSymbols": pertSymbols,
                                                                  "verbose": args.verbose,
                                                                  "initialGuesses": args.initialGuesses,
                                                                  "allSymbolsDict": {key: value for value, key in enumerate(allSymbols)}})
    if False:
        ##THIS IS FOR JASMINE TO MAKE PLOTS - IGNORE
        traceFreeEnergyMinimumInst.plotPotential(benchmark)
        exit()
    minimizationResult = traceFreeEnergyMinimumInst.traceFreeEnergyMinimum(benchmark)
  
    filename = f"{args.resultsDirectory}/BM_{benchmark['bmNumber']}"
    
    from pathlib import Path
    Path(args.resultsDirectory).mkdir(parents = True, exist_ok = True)
    if args.bSave:
        if args.verbose:
            print(f"Saving {benchmark['bmNumber']} to {filename}")
        open(f"{filename}.json", "w").write(json.dumps(minimizationResult, indent = 4))
      
    if args.bPlot:
        if args.verbose:
            print(f"Plotting {benchmark['bmNumber']}")

        from ThreeHiggs.PlotResult import plotData
        plotData(minimizationResult, benchmark['bmNumber'], args.loopOrder, filename)

    if args.bProcessMin:
        if args.verbose:
            print(f"Processing {benchmark['bmNumber']} to {filename+'_interp'}")
        from ThreeHiggs.ProcessMinimization import interpretData
        open(f"{filename}_interp.json", "w").write(json.dumps(interpretData(minimizationResult,
                                                                        benchmark["bmNumber"],
                                                                        benchmark["bmInput"]),
                                                         indent = 4))
from ThreeHiggs.GetLines import getLines        
def minimization(args):
    pythonisedExpressions = json.load(open(args.pythonisedExpressionsFile, "r"))
    allSymbols = pythonisedExpressions["allSymbols"]["allSymbols"]

    variableSymbols =  getLines( "Data/Variables/LagranianSymbols.json", mode = "json") 
    
    from ThreeHiggs.ParsedExpression import (ParsedExpressionSystem,
                                             ParsedExpressionSystemArray,
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
                                            args.verbose,
                                            nloptInst,
                                            ParsedExpressionSystemArray(
                                                pythonisedExpressions["vectorMassesSquared"]["expressions"],
                                                allSymbols,
                                                pythonisedExpressions["vectorMassesSquared"]["fileName"],
                                            ),
                                            ParsedExpressionSystemArray(
                                                pythonisedExpressions["vectorShortHands"]["expressions"],
                                                allSymbols,
                                                pythonisedExpressions["vectorShortHands"]["fileName"],
                                            ),
                                            pythonisedExpressions["scalarPermutationMatrix"],
                                            (MassMatrix(pythonisedExpressions["scalarMassMatrixUpperLeft"]["expressions"], 
                                                        pythonisedExpressions["scalarMassMatrixUpperLeft"]["fileName"]), 
                                             MassMatrix(pythonisedExpressions["scalarMassMatrixBottomRight"]["expressions"],
                                                        pythonisedExpressions["scalarMassMatrixBottomRight"]["fileName"])),
                                            RotationMatrix(pythonisedExpressions["scalarRotationMatrix"]["expressions"],
                                                           pythonisedExpressions["scalarRotationMatrix"]["fileName"]),
                                            ParsedExpressionSystem(pythonisedExpressions["veff"]["expressions"],
                                                                   pythonisedExpressions["veff"]["fileName"]),
                                            ParsedExpressionSystemArray(pythonisedExpressions["veffArray"]["expressions"], 
                                                                        allSymbols, 
                                                                        pythonisedExpressions["veffArray"]["fileName"]),
                                            allSymbols) 

    from ThreeHiggs.DimensionalReduction import DimensionalReduction
    from ThreeHiggs.ParsedExpression import ParsedExpressionSystemArray
    dimensionalReduction = DimensionalReduction(config = {
        "hardToSoft": ParsedExpressionSystemArray(
            pythonisedExpressions["hardToSoft"]["expressions"], 
            allSymbols,
            pythonisedExpressions["hardToSoft"]["fileName"],
        ),
        "softScaleRGE": ParsedExpressionSystemArray(
            pythonisedExpressions["softScaleRGE"]["expressions"],
            allSymbols,
            pythonisedExpressions["softScaleRGE"]["fileName"],
        ),
        "softToUltraSoft": ParsedExpressionSystemArray(
            pythonisedExpressions["softToUltraSoft"]["expressions"],
            allSymbols,
            pythonisedExpressions["softToUltraSoft"]["fileName"],
        )
    })
    
    with open(args.benchmarkFile) as benchmarkFile:
        minimizationDict = {"pertSymbols": frozenset(variableSymbols["fourPointSymbols"] + 
                                                     variableSymbols["yukawaSymbols"] + 
                                                     variableSymbols["gaugeSymbols"]), 
                            "effectivePotential": effectivePotential,
                            "dimensionalReduction": dimensionalReduction,
                            "betaFunction4DExpression": ParsedExpressionSystemArray(pythonisedExpressions["betaFunctions4D"]["expressions"], 
                                                                                    allSymbols, 
                                                                                    pythonisedExpressions["betaFunctions4D"]["fileName"]),
                            "args": args,
                            "allSymbols": allSymbols} 
        if args.bPool:
            from pathos.multiprocessing import Pool
            with Pool(args.cores) as pool:
                from ijson import items
                pool.map(_doMinimization, (minimizationDict | {"benchmark": item} for item in items(benchmarkFile, "item", use_float = True)))

        else:
            for parameters in (minimizationDict | {"benchmark": item} for item in json.load(benchmarkFile)):
                _doMinimization(parameters)
