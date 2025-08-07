import argparse
import multiprocessing
from enum import IntEnum
import json
class Stages(IntEnum):
    convertMathematica = 0
    generateBenchmark = 1
    doMinimization = 2

    @staticmethod
    def fromString(*args, **kwargs):
        if args:
            return Stages[args[0]]

        if kwargs:
            return Stages[kwargs["default"]]

class UserInput(argparse.ArgumentParser):
    def __init__(self):
        super().__init__()
        self.add_argument('--config', action = 'store', dest = 'config', type = str,
                          help = "Str: Load cmd line args from json")
        
        self.add_argument('--loopOrder', action = 'store', default = 2, dest = 'loopOrder', type = int, choices = [1, 2],
                          help = "Int: Specify the order to compute the effective potential to")
        
        ## Should probably be made into a int rather than bool to allow for levels of verbosisty 
        self.add_argument('--verbose', action = 'store_true', default = False, dest = 'verbose',
                          help = 'Bool: If activated code will print as it progresses and do consistency checks')
        
        self.add_argument('--cores', action = 'store', default = 1, dest = 'cores', type = int, choices = list(range(1, multiprocessing.cpu_count() + 1)),
                          help = "Int: Specify how many cores pool uses to compute benchmarks",
                          metavar='')
        self.add_argument('--bSave', action = 'store_true', default=False, dest = 'bSave',  
                          help = "Bool: If activated the results of the minimisation will be saved")
        
        self.add_argument('--bPlot', action = 'store_true', default = False, dest = 'bPlot',
                          help = 'Bool: If activated a plot of the global min of the potential vs T is made')
        
        self.add_argument('--bProcessMin', action = 'store_true', default = False, dest = 'bProcessMin',
                          help = 'Bool: If activated phase transitions are identified from the minimisation results')
        
        self.add_argument('--resultsDirectory', action = 'store', default="Results", dest = 'resultsDirectory',  
                          help = "Str: Location to save files",
                          metavar='')
        
        self.add_argument('--absGlobalTolerance', action = 'store', default = 0.5, type = float, dest = 'absGlobalTolerance',
        metavar='')
        self.add_argument('--relGlobalTolerance', action = 'store', default = 0.5, type = float, dest = 'relGlobalTolerance',
        metavar='')
        self.add_argument('--absLocalTolerance', action = 'store', default = 1e-2, type = float, dest = 'absLocalTolerance',
        metavar='')
        self.add_argument('--relLocalTolerance', action = 'store', default = 1e-3, type = float, dest = 'relLocalTolerance',
        metavar='')
        
        self.add_argument('--varLowerBounds', nargs = "*", action = 'store', default = [-60, 1e-4, 1e-4], type = float, dest = 'varLowerBounds',
        metavar='')
        self.add_argument('--varUpperBounds', nargs = "*", action = 'store', default = [60, 60, 60], type = float, dest = 'varUpperBounds',
        metavar='')
        
        self.add_argument('--initialGuesses', 
                          nargs = "*", 
                          action = 'store', 
                          default = [[0.1,0.1,0.1], 
                                     [5,5,5],
                                     [-5,5,5], 
                                     [5,5,5],
                                     [-5,5,5],
                                     [0.1,0.1, 10],
                                     [0.1,0.1, 20],
                                     [40,40,40],
                                     [-40,40,40], 
                                     [59,59,59], 
                                     [-59,59,59]], 
                          type = tuple, dest = 'initialGuesses',
        metavar='')
        
        self.add_argument('--TRangeStart', action = 'store', default = 50, type = float, dest = 'TRangeStart',
        metavar='')
        self.add_argument('--TRangeEnd', action = 'store', default = 200, type = float, dest = 'TRangeEnd',
        metavar='')
        self.add_argument('--TRangeStepSize', action = 'store', default = 1, type = float, dest = 'TRangeStepSize',
        metavar='')

        self.add_argument('--firstStage', type = Stages.fromString, default = "convertMathematica", dest = 'firstStage',
        metavar='')
        self.add_argument('--lastStage', type = Stages.fromString, default = "doMinimization", dest = 'lastStage',
        metavar='')
        
        

        self.add_argument('--benchmarkFile', 
                  action = 'store', 
                  default = "ThreeHiggs/Data/Z2_3HDM/Benchmarks/handPicked.json",
                  dest = 'benchmarkFile',
                  metavar='')
        
        self.add_argument("--benchmarkType", action = "store",  dest = "benchmarkType", default = "handPicked",
                        choices = ["load", "handPicked", "random", "randomSSS"],
                        help = "Str: Specify the mode to generate bm with.")
        
        self.add_argument("--randomNum",type = int, action = "store",  dest = "randomNum",
                        help = "Int: Specify how many random bm to generate.")
        self.add_argument("--prevResultDir",type = str, action = "store",  dest = "prevResultDir",
                        help = "str: SLoad previous results to do a strong sub set with.")
        
        self.add_argument('--firstBenchmark', type = int, default = 0, dest = 'firstBenchmark',
        metavar='')
        from sys import maxsize
        self.add_argument('--lastBenchmark', type = int, default = maxsize, dest = 'lastBenchmark',
        metavar='')

        self.add_argument('--allSymbolsFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/Variables/allSymbols.json",
                          dest = 'allSymbolsFile',
                          metavar='')
        
        self.add_argument('--loFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/Veff_LO.txt",
                          dest = 'loFile',
                          metavar='')

        self.add_argument('--nloFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/Veff_NLO.txt",
                          dest = 'nloFile',
                          metavar='')

        self.add_argument('--nnloFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/Veff_NNLO.txt",
                          dest = 'nnloFile',
                          metavar='')

        self.add_argument('--betaFunctions4DFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/HardToSoft/BetaFunctions4D.txt",
                          dest = 'betaFunctions4DFile',
                          metavar='')

        self.add_argument('--vectorMassesSquaredFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/vectorMasses.txt",
                          dest = 'vectorMassesSquaredFile',
                          metavar='')

        self.add_argument('--vectorShortHandsFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/vectorShorthands.txt",
                          dest = 'vectorShortHandsFile',
                          metavar='')

        self.add_argument('--hardToSoftFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/HardToSoft/softScaleParams_NLO.txt",
                          dest = 'hardToSoftFile',
                          metavar='')

        self.add_argument('--softScaleRGEFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/HardToSoft/softScaleRGE.txt",
                          dest = 'softScaleRGEFile',
                          metavar='')

        self.add_argument('--softToUltraSoftFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/SoftToUltraSoft/ultrasoftScaleParams_NLO.txt",
                          dest = 'softToUltraSoftFile',
                          metavar='')

        self.add_argument('--scalarPermutationMatrixFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarPermutationMatrix.txt",
                          dest = 'scalarPermutationFile',
                          metavar='')

        self.add_argument('--scalarRotationMatrixFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarRotationMatrix.txt",
                          dest = 'scalarRotationFile',
                          metavar='')

        self.add_argument('--scalarMassMatrixUpperLeftFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarMassMatrix_upperLeft.txt",
                          dest = 'scalarMassMatrixUpperLeftFile',
                          metavar='')

        self.add_argument('--scalarMassMatrixUpperLeftDefinitionsFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarMassMatrix_upperLeft_definitions.txt",
                          dest = 'scalarMassMatrixUpperLeftDefinitionsFile',
                          metavar='')

        self.add_argument('--scalarMassMatrixBottomRightFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarMassMatrix_bottomRight.txt",
                          dest = 'scalarMassMatrixBottomRightFile',
                          metavar='')

        self.add_argument('--scalarMassMatrixBottomRightDefinitionsFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarMassMatrix_bottomRight_definitions.txt",
                          dest = 'scalarMassMatrixBottomRightDefinitionsFile',
                          metavar='')
        
        self.add_argument('--lagranianVariablesFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/Variables/LagranianSymbols.json",
                          dest = 'lagranianVariables',
                          metavar='')
        
        self.add_argument('--pythonisedExpressionsFile', 
                          action = 'store', 
                          default = "ThreeHiggs/Data/Z2_3HDM/pythonisedExpressionsFile.json",
                          dest = 'pythonisedExpressionsFile',
                          metavar='')
    def parse(self):
        userArg = super().parse_args()
        if userArg.config:
            userConfig = json.load(open(super().parse_args().config, "r"))
            unexpectedKeys = [userKey for userKey in set(userConfig.keys()) if userKey not in set(vars(userArg).keys())]
            if len(unexpectedKeys) > 0:
                print(f"User config file has unexpected key(s): {unexpectedKeys},\nExiting")
                exit(-1)
            self.set_defaults(**userConfig)
        return super().parse_args()

