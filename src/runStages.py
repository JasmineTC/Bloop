import sys
sys.path.append("../src/")

from ThreeHiggs.UserInput import UserInput, Stages
from Veff_generation import generate_veff_module, compile_veff_submodule

args = UserInput().parse()

if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    if args.verbose:
        print("Convert Mathematica stage started")

    from ThreeHiggs.PythoniseMathematica import pythoniseMathematica
    pythoniseMathematica(args)

    from ThreeHiggs.ConvertMathematica import convertMathematica
    generate_veff_module(args)
    compile_veff_submodule()
    convertMathematica(args)

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

