import argparse
import multiprocessing
from enum import IntEnum
import json
from sys import maxsize

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
        self.add_argument('--config', 
                          action = 'store', 
                          type = str,
                          help = "Str: Load cmd line args from json"
                          )
        
        self.add_argument('--loopOrder', 
                          action = 'store', 
                          default = 1, 
                          type = int, 
                          choices = [1, 2],
                          help = "Int: Specify the order to compute the effective potential to"
                          )
        
        ## Should probably be made into a int rather than bool to allow for levels of verbosisty 
        self.add_argument('--verbose', 
                          action = 'store_true', 
                          default = False,
                          help = 'Bool: If activated code will print as it progresses'
                          )
        
        self.add_argument('--bPool', 
                          action = 'store_true', 
                          default = False,
                          help = 'Bool: If activated code will run in parallel using number of cores set by --cores'
                          )
        
        self.add_argument('--cores', 
                          action = 'store', 
                          default = 1, 
                          choices = list(range(1, multiprocessing.cpu_count() + 1)),
                          type = int,
                          help = "Int: Specify how many cores pool uses to compute benchmarks"
                          )
        
        self.add_argument('--bSave', 
                          action = 'store_true', 
                          default=False,  
                          help = "Bool: If activated the results of the minimisation will be saved"
                          )
        
        self.add_argument('--bPlot', 
                          action = 'store_true', 
                          default = False,
                          help = 'Bool: If activated a plot of the global min of the potential vs T is made'
                          )
        
        self.add_argument('--bProcessMin', 
                          action = 'store_true', 
                          default = False,
                          help = 'Bool: If activated phase transitions are identified from the minimisation results'
                          )
        
        self.add_argument('--resultsDirectory', 
                          action = 'store', 
                          default="Results",  
                          help = "Str: Location to save files"
                          )
        
        self.add_argument('--absGlobalTolerance', 
                          action = 'store', 
                          default = 0.5, 
                          type = float
                          )
        
        self.add_argument('--relGlobalTolerance', 
                          action = 'store', 
                          default = 0.5, 
                          type = float
                          )
        
        self.add_argument('--absLocalTolerance', 
                          action = 'store', 
                          default = 1e-2, 
                          type = float
                          )
        
        self.add_argument('--relLocalTolerance', 
                          action = 'store', 
                          default = 1e-3, 
                          type = float
                          )
        
        self.add_argument('--varLowerBounds', 
                          nargs = "*", 
                          action = 'store', 
                          default = [-60, 1e-4, 1e-4], 
                          type = float
                          )
        
        self.add_argument('--varUpperBounds', 
                          nargs = "*", 
                          action = 'store', 
                          default = [60, 60, 60], 
                          type = list
                          )
        
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
                          type = list
                          )
        
        self.add_argument('--TRangeStart', 
                          action = 'store', 
                          default = 50, 
                          type = float
                          )
        
        self.add_argument('--TRangeEnd', 
                          action = 'store', 
                          default = 200, 
                          type = float
                          )
        
        self.add_argument('--TRangeStepSize', 
                          action = 'store', 
                          default = 1,
                          type = float
                          )

        self.add_argument('--firstStage', 
                          default = "convertMathematica",
                          type = Stages.fromString
                          )
        
        self.add_argument('--lastStage', 
                          default = "doMinimization",
                          type = Stages.fromString
                          )
        
        self.add_argument('--benchmarkFile', 
                          action = 'store', 
                          default = "ThreeHiggs/Data/Z2_3HDM/Benchmarks/handPicked.json"
                          )
        
        self.add_argument("--benchmarkType",
                          action = "store", 
                          default = "handPicked",
                          choices = ["load", "handPicked", "random", "randomSSS"],
                          help = "Str: Specify the mode to generate bm with."
                          )
        
        self.add_argument("--randomNum",
                          type = int, 
                          action = "store",
                          help = "Int: Specify how many random bm to generate."
                          )
        
        self.add_argument("--prevResultDir",
                          action = "store",
                          help = "str: SLoad previous results to do a strong sub set with."
                          )
        
        self.add_argument('--firstBenchmark', 
                          default = 0,
                          type = int
                          )

        self.add_argument('--lastBenchmark', 
                          default = maxsize,
                          type = int
                          )

        self.add_argument('--BmGeneratorFile', 
                          action = 'store', 
                          default = "ThreeHiggs.Z2_ThreeHiggsBmGenerator",
                          )
        
        self.add_argument('--boundedConditions', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/Misc/bounded.txt",
                          )
 
        self.add_argument('--allSymbolsFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/Variables/allSymbols.json"
                          )
        
        self.add_argument('--loFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/Veff_LO.txt"
                          )

        self.add_argument('--nloFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/Veff_NLO.txt"
                          )

        self.add_argument('--nnloFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/Veff_NNLO.txt"
                          )

        self.add_argument('--betaFunctions4DFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/HardToSoft/BetaFunctions4D.txt"
                          )

        self.add_argument('--vectorMassesSquaredFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/vectorMasses.txt"
                          )

        self.add_argument('--vectorShortHandsFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/vectorShorthands.txt"
                          )

        self.add_argument('--hardToSoftFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/HardToSoft/softScaleParams_NLO.txt"
                          )

        self.add_argument('--softScaleRGEFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/HardToSoft/softScaleRGE.txt"
                          )

        self.add_argument('--softToUltraSoftFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/SoftToUltraSoft/ultrasoftScaleParams_NLO.txt"
                          )

        self.add_argument('--scalarPermutationMatrixFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarPermutationMatrix.txt"
                          )

        self.add_argument('--scalarRotationMatrixFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarRotationMatrix.txt"
                          )
        
        
        self.add_argument('--scalarMassMatricesFiles', 
                          nargs = "*", 
                          action = 'store', 
                          default = ["Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarMassMatrix_upperLeft.txt",
                                     "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarMassMatrix_bottomRight.txt"], 
                          type = list
                          )
        
        self.add_argument('--scalarMassMatricesDefinitionsFiles', 
                          nargs = "*", 
                          action = 'store', 
                          default = ["Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarMassMatrix_upperLeft_definitions.txt",
                                     "Data/Z2_3HDM/ModelFiles/EffectivePotential/scalarMassMatrix_bottomRight_definitions.txt"], 
                          type = list
                          )
        
        self.add_argument('--lagranianVariablesFile', 
                          action = 'store', 
                          default = "Data/Z2_3HDM/ModelFiles/Variables/LagranianSymbols.json"
                          )
        
        self.add_argument('--pythonisedExpressionsFile', 
                          action = 'store', 
                          default = "ThreeHiggs/Data/Z2_3HDM/pythonisedExpressionsFile.json"
                          )
        
    noMetaVar = {'store_true', 'store_false', 'help', 'version'} 
    def add_argument(self, *args, **kwargs):
        if (args[0].startswith('-') and 
        kwargs.get('action') not in self.noMetaVar):
            kwargs['metavar'] = ''
        return super().add_argument(*args, **kwargs)
    
    def parse(self):
        userArg = super().parse_args()
        if userArg.config:
            userConfig = json.load(open(super().parse_args().config, "r"))
            unexpectedKeys = [userKey for userKey in userConfig.keys() if userKey not in set(vars(userArg).keys())]
            if len(unexpectedKeys) > 0:
                print(f"User config file has unexpected key(s):\n {unexpectedKeys},\nExiting")
                exit(-1)
            self.set_defaults(**userConfig)
        return super().parse_args()

