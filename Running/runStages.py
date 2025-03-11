
from ThreeHiggs.UserInput import UserInput, Stages
args = UserInput().parse()
# print(args)
# open("test.json", "w").write(json.dumps(vars(args), indent = 4))
# loadArgs = (json.load(open("test.json", "r")))
# print(vars(args) | loadArgs)
# print(UserInput().parse())
if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    from ThreeHiggs.ConvertMathematica import convertMathematica
    convertMathematica(args)

if args.firstStage <= Stages.generateBenchmark <= args.lastStage:
    from ThreeHiggs.BmGenerator import generateBenchmarks
    generateBenchmarks(args.benchmarkFile, args.benchmarkType, args.randomNum, args.prevResultDir)

if args.firstStage <= Stages.doMinimization <= args.lastStage:
    from ThreeHiggs.DoMinimization import minimization
    minimization(args)
    
