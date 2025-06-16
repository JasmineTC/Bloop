import sys
sys.path.append("../src/")

from ThreeHiggs.UserInput import UserInput, Stages
from Veff_generation import generate_veff_module, compile_veff_submodule

args = UserInput().parse()

if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    from ThreeHiggs.ConvertMathematica import convertMathematica
    generate_veff_module(args)
    compile_veff_submodule()
    convertMathematica(args)

if args.firstStage <= Stages.generateBenchmark <= args.lastStage:
    from ThreeHiggs.BmGenerator import generateBenchmarks
    generateBenchmarks(args.benchmarkFile, args.benchmarkType, args.randomNum, args.prevResultDir)

if args.firstStage <= Stages.doMinimization <= args.lastStage:
    from ThreeHiggs.DoMinimization import minimization
    minimization(args)
    
