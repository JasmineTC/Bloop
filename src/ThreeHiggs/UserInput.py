import argparse
import multiprocessing
from enum import IntEnum

class Stages(IntEnum):
    convertMathematica = 0
    minimization = 1
    plot = 2

    @staticmethod
    def fromString(*args, **kwargs):
        if args:
            return Stages[args[0]]

        if kwargs:
            return Stages[kwargs["default"]]

class UserInput(argparse.ArgumentParser):
    def __init__(self):
        super().__init__()
        self.add_argument('-l', '--loopOrder', action = 'store', default = 2, dest = 'loopOrder', type = int, choices = [0, 1, 2],
                          help = "Int: Specify the order to compute the effective potential to")
        
        ## Should probably be made into a int rather than bool to allow for levels of verbosisty 
        self.add_argument('-v', '--bVerbose', action = 'store_true', default = False, dest = 'bVerbose',
                          help = 'Bool: If activated code will print as it progresses and do consistency checks')
        
        
        
        self.add_argument('-p', '--bPool', action = 'store_true', default = False, dest = 'bPool', 
                          help = "Bool: Specify if pool should be used to compute benchmarks in parellel")
        
        self.add_argument('-c', '--cores', action = 'store', default = 1, dest = 'cores', type = int, choices = list(range(1, multiprocessing.cpu_count() + 1)),
                          help = "Int: Specify how many cores pool uses to compute benchmarks")



        self.add_argument('-s', '--bSave', action = 'store_true', default=False, dest = 'bSave',  
                          help = "Bool: If activated the results of the minimisation will be saved")
        
        self.add_argument('--bPlot', action = 'store_true', default = False, dest = 'bPlot',
                          help = 'Bool: If activated a plot of the global min of the potential vs T is made')
        
        self.add_argument('--bProcessMin', action = 'store_true', default = False, dest = 'bProcessMin',
                          help = 'Bool: If activated phase transitions are identified from the minimisation results')
        
        self.add_argument('--resultsDirectory', action = 'store', default="Results", dest = 'resultsDirectory',  
                          help = "Str: Location to save files")
        
        self.add_argument('--bAbsMass', action = 'store_true', default = False, dest = 'bAbsMass',
                          help = "Bool: If activated the abs value of masses is used")
        
        self.add_argument('--absGlobalTolerance', action = 'store', default = 0.5, type = float, dest = 'absGlobalTolerance')
        self.add_argument('--relGlobalTolerance', action = 'store', default = 0.5, type = float, dest = 'relGlobalTolerance')
        self.add_argument('--absLocalTolerance', action = 'store', default = 1e-2, type = float, dest = 'absLocalTolerance')
        self.add_argument('--relLocalTolerance', action = 'store', default = 1e-3, type = float, dest = 'relLocalTolerance')
        
        self.add_argument('--v1Bounds', nargs = 2, action = 'store', default = [-60, 60], type = float, dest = 'v1Bounds')
        self.add_argument('--v2Bounds', nargs = 2, action = 'store', default = [1e-4, 60], type = float, dest = 'v2Bounds')
        self.add_argument('--v3Bounds', nargs = 2, action = 'store', default = [1e-4, 60], type = float, dest = 'v3Bounds')
        
        self.add_argument('--TRangeStart', action = 'store', default = 50, type = float, dest = 'TRangeStart')
        self.add_argument('--TRangeEnd', action = 'store', default = 200, type = float, dest = 'TRangeEnd')
        self.add_argument('--TRangeStepSize', action = 'store', default = 1, type = float, dest = 'TRangeStepSize')
        
        self.add_argument('--diagAlgo', action = 'store', default = "numpy", type = str, dest = 'diagAlgo',
                          help = 'Str: Specify what diagonalisation algorithm used for mass matrices, \
                          options are numpy, scipy and mpmath')
        self.add_argument('--minimizationAlgo', action = 'store', default="combo", dest = 'minimizationAlgo',  
                          help = "Used to specify which algothrym the minimizer uses, options are combo, directGlobal, BOBYQA and scipy")    

        self.add_argument('--firstStage', type = Stages.fromString, default = "convertMathematica", dest = 'firstStage')
        self.add_argument('--lastStage', type = Stages.fromString, default = "plot", dest = 'lastStage')


        self.add_argument('--benchmarkFile', 
                  action = 'store', 
                  default = "Benchmarks/Benchmarks_3HDM.json",#"Benchmarks/ScanDict.json",
                  dest = 'benchmarkFile')
        
        self.add_argument('--firstBenchmark', type = int, default = 0, dest = 'firstBenchmark')
        from sys import maxsize
        self.add_argument('--lastBenchmark', type = int, default = maxsize, dest = 'lastBenchmark')
        
        self.add_argument('--loFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/Veff_LO.txt",
                          dest = 'loFile')

        self.add_argument('--nloFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/Veff_NLO.txt",
                          dest = 'nloFile')

        self.add_argument('--nnloFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/Veff_NNLO.txt",
                          dest = 'nnloFile')

        self.add_argument('--vectorMassesSquaredFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/vectorMasses.txt",
                          dest = 'vectorMassesSquaredFile')

        self.add_argument('--vectorShortHandsFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/vectorShorthands.txt",
                          dest = 'vectorShortHandsFile')

        self.add_argument('--hardToSoftFile', 
                          action = 'store', 
                          default = "Data/HardToSoft/softScaleParams_NLO.txt",
                          dest = 'hardToSoftFile')

        self.add_argument('--softScaleRGEFile', 
                          action = 'store', 
                          default = "Data/HardToSoft/softScaleRGE.txt",
                          dest = 'softScaleRGEFile')

        self.add_argument('--softToUltraSoftFile', 
                          action = 'store', 
                          default = "Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt",
                          dest = 'softToUltraSoftFile')

        self.add_argument('--scalarPermutationMatrixFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/scalarPermutationMatrix.txt",
                          dest = 'scalarPermutationFile')

        self.add_argument('--scalarRotationMatrixFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/scalarRotationMatrix.txt",
                          dest = 'scalarRotationFile')

        self.add_argument('--scalarMassMatrixUpperLeftFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/scalarMassMatrix_upperLeft.txt",
                          dest = 'scalarMassMatrixUpperLeftFile')

        self.add_argument('--scalarMassMatrixUpperLeftDefinitionsFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/scalarMassMatrix_upperLeft_definitions.txt",
                          dest = 'scalarMassMatrixUpperLeftDefinitionsFile')

        self.add_argument('--scalarMassMatrixBottomRightFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/scalarMassMatrix_bottomRight.txt",
                          dest = 'scalarMassMatrixBottomRightFile')

        self.add_argument('--scalarMassMatrixBottomRightDefinitionsFile', 
                          action = 'store', 
                          default = "Data/EffectivePotential_threeFields/scalarMassMatrix_bottomRight_definitions.txt",
                          dest = 'scalarMassMatrixBottomRightDefinitionsFile')

        self.add_argument('--parsedExpressionsFile', 
                  action = 'store', 
                  default = "parsedExpressions.json",
                  dest = 'parsedExpressionsFile')

    ##Used to check userinputs are valid, mostly done with the choice keyword above now though
    def parse(self):
        return super().parse_args()

