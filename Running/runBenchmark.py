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
              "scalarMassMatrixUpperLeft": parseMassMatrix(getLines(args.scalarMassMatrixUpperLeftFile)),
              "scalarMassMatrixUpperLeftDefinitions": parseExpressionSystem(getLines(args.scalarMassMatrixUpperLeftDefinitionsFile)),
              "scalarMassMatrixBottomRight": parseMassMatrix(getLines(args.scalarMassMatrixBottomRightFile)),
              "scalarMassMatrixBottomRightDefinitions": parseExpressionSystem(getLines(args.scalarMassMatrixBottomRightDefinitionsFile)),
              "scalarRotationMatrix": parseRotationMatrix(getLines(args.scalarRotationFile)),
              "veff": parseExpressionSystem(veffLines)},
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

        effectivePotential = EffectivePotential(['v1', 'v2', 'v3'],
                                                True,
                                                ParsedExpressionSystem(parsedExpressions["vectorMassesSquared"]),
                                                ParsedExpressionSystem(parsedExpressions["vectorShortHands"]),
                                                parsedExpressions["scalarPermutationMatrix"]["matrix"],
                                                [MassMatrix(parsedExpressions["scalarMassMatrixUpperLeft"]["matrix"], 
                                                            ParsedExpressionSystem(parsedExpressions["scalarMassMatrixUpperLeftDefinitions"])),
                                                 MassMatrix(parsedExpressions["scalarMassMatrixBottomRight"]["matrix"], 
                                                            ParsedExpressionSystem(parsedExpressions["scalarMassMatrixBottomRightDefinitions"]))],
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
    from ThreeHiggs.ParameterMatching import ParameterMatching
    dimensionalReduction = DimensionalReduction(ParameterMatching(getLines(args.hardToSoftFile)),
                                                ParameterMatching(getLines(args.softScaleRGEFile)),
                                                ParameterMatching(getLines(args.softToUltraSoftFile)))

    ## Model object setup + load matching relations
    from ThreeHiggs.GenericModel import GenericModel
    model3HDM = GenericModel(effectivePotential, dimensionalReduction)
    
    #model3HDM.dimensionalReduction.setupHardToSoftMatching(hardToSoftFile, softScaleRGEFile)
    #model3HDM.dimensionalReduction.setupSoftToUltrasoftMatching(softToUltrasoftFile)
    
    # model3HDM.effectivePotential.minimizer.setTolerances(args.absGlobalTolerance,
    #                                                      args.relGlobalTolerance, 
    #                                                      args.absLocalTolerance, 
    #                                                      args.relLocalTolerance)
    
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

