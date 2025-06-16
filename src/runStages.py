<<<<<<< HEAD:src/runStages.py
from ThreeHiggs.UserInput import UserInput, Stages
=======
import sys
sys.path.append("../src/")

from ThreeHiggs.UserInput import UserInput, Stages
from Veff_generation import generate_veff_module, compile_veff_submodule
>>>>>>> bfe7c1b (Add proof of concept implementation):Running/runStages.py

args = UserInput().parse()

if args.firstStage <= Stages.convertMathematica <= args.lastStage:
<<<<<<< HEAD:src/runStages.py
    if args.verbose:
        print("Convert Mathematica stage started")

    from ThreeHiggs.PythoniseMathematica import pythoniseMathematica
    pythoniseMathematica(args)
=======
    from ThreeHiggs.ConvertMathematica import convertMathematica
    generate_veff_module(args)
    compile_veff_submodule()
    convertMathematica(args)
>>>>>>> bfe7c1b (Add proof of concept implementation):Running/runStages.py

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
