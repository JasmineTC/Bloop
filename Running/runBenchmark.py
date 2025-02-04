import json 
from typing import Generator
import decimal

## This avoids floating point error in T gotten by np.arange or linspace
## However one must be careful as 1 = decimal.Decimal(1.000000000000001) 
def drange(start: float, end: float, jump: str) -> Generator:
    start =  decimal.Decimal(start) 
    while start <= end:
        yield float(start)
        start += decimal.Decimal(jump)

def doMinimization(parameters):
    benchmark = parameters["benchmark"] if "benchmark" in parameters else None
    effectivePotential = parameters["effectivePotential"] if "effectivePotential" in parameters else None
    dimensionalReduction = parameters["dimensionalReduction"] if "dimensionalReduction" in parameters else None
    pertSymbols = parameters["pertSymbols"] if "pertSymbols" in parameters else None

    if args.bVerbose:
        print(f"Starting benchmark: {benchmark['bmNumber']}")

    if not args.firstBenchmark <= benchmark['bmNumber'] <= args.lastBenchmark:
        if args.bVerbose:
            print(f"Benchmark {benchmark['bmNumber']} has been rejected as outside benchmark range.")

        return
    
    from ThreeHiggs.TransitionFinder import TraceFreeEnergyMinimum
    traceFreeEnergyMinimumInst = TraceFreeEnergyMinimum(config = {"effectivePotential":effectivePotential, 
                                                "dimensionalReduction": dimensionalReduction, 
                                                "TRange": tuple(drange(args.TRangeStart, 
                                                                       args.TRangeEnd, 
                                                                       str(args.TRangeStepSize))),
                                                "pertSymbols": pertSymbols,
                                                "bVerbose": args.bVerbose,
                                                "initialGuesses": args.initialGuesses})
    
    
    minimizationResult = traceFreeEnergyMinimumInst.traceFreeEnergyMinimum(benchmark)
  
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

## Adding a mode is not ideal but idk what else to do
def getLines(relativePathToResource, mode = "default"):
    ## fallback to hardcoded package name if the __package__ call fails
    packageName = __package__ or "ThreeHiggs"

    from importlib.resources import files
    path = files(packageName) / relativePathToResource
    
    if mode == "json":
        return json.load(open(path, "r"))
    return open(path, "r" , encoding = "utf-8").readlines()

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
    

    json.dump({"vectorMassesSquared": parseExpressionSystem(getLines(args.vectorMassesSquaredFile)),
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
        if args.bPool:
            from pathos.multiprocessing import Pool
            with Pool(args.cores) as pool:
                from ijson import items
                pool.map(doMinimization, ({"benchmark": item,
                                           "effectivePotential": effectivePotential,
                                           "dimensionalReduction": dimensionalReduction } for item in items(benchmarkFile, "item", use_float = True)))

        else:
            for parameters in ({"pertSymbols": frozenset(variableSymbols["fourPointSymbols"] + 
                                                   variableSymbols["yukawaSymbols"] + 
                                                   variableSymbols["gaugeSymbols"]), 
                                "benchmark": item,
                                "effectivePotential": effectivePotential,
                                "dimensionalReduction": dimensionalReduction} for item in json.load(benchmarkFile)):
                doMinimization(parameters)

from ThreeHiggs.UserInput import UserInput
userinput = UserInput()
args = userinput.parse()

from ThreeHiggs.UserInput import Stages
if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    convertMathematica(args)

if args.firstStage <= Stages.minimization <= args.lastStage:
   minimization(args)