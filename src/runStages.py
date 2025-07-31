from ThreeHiggs.UserInput import UserInput, Stages

args = UserInput().parse()

if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    if args.verbose:
        print("Convert Mathematica stage started")

    from ThreeHiggs.PythoniseMathematica import pythoniseMathematica
    pythoniseMathematica(args)

if args.firstStage <= Stages.generateBenchmark <= args.lastStage:
    if args.verbose:
        print("Benchmark generation stage started")
        
    from ThreeHiggs.BmGenerator import generateBenchmarks
    generateBenchmarks(args)

if args.firstStage <= Stages.doMinimization <= args.lastStage:
    if args.verbose:
        print("Minimization stage started")
    from ThreeHiggs.DoMinimization import minimization
    minimization(args)
