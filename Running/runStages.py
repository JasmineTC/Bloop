
from ThreeHiggs.UserInput import UserInput, Stages
args = UserInput().parse()

if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    from ThreeHiggs.ConvertMathematica import convertMathematica
    convertMathematica(args)

if args.firstStage <= Stages.generateBenchmark <= args.lastStage:
    from ThreeHiggs.BmGenerator import generateBenchmarks
    generateBenchmarks(args)

if args.firstStage <= Stages.doMinimization <= args.lastStage:
    from ThreeHiggs.DoMinimization import minimization
    minimization(args)
    
