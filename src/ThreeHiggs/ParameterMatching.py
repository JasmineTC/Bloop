import numpy as np
import sympy
from sympy import pi
from sympy.parsing.mathematica import parse_mathematica

from typing import Callable

from CommonUtils import getPackagedDataPath
from CommonUtils import replaceGreekSymbols


## This would make a good abstract class
class ParameterMatching:


    ## stores lambdas for evaluating params in the soft scale theory
    matchingRelations: dict[str, Callable]

    ## This list specifies names for parameters that enter the matching relations. 
    # Needs to match symbol names in the file that defines matching relations (there can be extra symbols like RG scale)
    parameterNames = [ 'T', 'Lb', 'Lf', 'g1', 'g2', 'g3' ]  


    def __init__(self):


        self.__defineSymbols(self.parameterNames)
        
        self.__parsedMatchingRelations = self.parseMatchingRelations()

        self.defineMatchingRelations()


    def getMatchedParams(self, inputParams: dict[str, float]) -> dict[str, float]:

        ## put input params in an array so that matching relation lambdas can be evaluated
        inParamList = self._paramDictToOrderedList(inputParams)

        matchedParams = {}

        for key, value in self.matchingRelations.items():
            matchedParams[key] = value(*inParamList) ## Unpack because the lambdas don't take lists

        return matchedParams


    def defineMatchingRelations(self) -> None:

        ## I'll just have all matching relations take the same list of args
        argumentSymbols = self.parameterNames

        self.matchingRelations = {}

        self.matchingRelations["g1sq"] = self._exprToFunction(self.__parsedMatchingRelations['g13d^2'], argumentSymbols)
        self.matchingRelations["g2sq"] = self._exprToFunction(self.__parsedMatchingRelations['g23d^2'], argumentSymbols)

    
    ## Convert a dict of input params to list that has same ordering as that defined in defineMatchingRelations() 
    def _paramDictToOrderedList(self, params: dict[str, float]) -> list[float]:

        ## TODO would be much better to have a "params" class that handles all parameter stuff
        outList = [None] * len(self.parameterNames)
        
        for i in range(len(outList)):
            outList[i] = params[ self.parameterNames[i] ]

        
        return outList
        

    def __defineSymbols(self, symbolList: list[str]) -> None:

        self.symbols = {}

        ## TODO How to clear sympy definitions?!?!

        for s in symbolList:
            self.symbols[s] = sympy.Symbol(s)


    def parseMatchingRelations(self) -> dict[str, sympy.Expr]:
        
        ## Use this if changing this code to a proper package:
        #filePath = getPackagedDataPath("ThreeHiggs.Data", "softScaleParams.txt")
        
        filePath = "Data/softScaleParams.txt"

        with open(filePath, "r") as file:
            expressions = file.readlines()

        parsedExpressions = {}
        for line in expressions:
            try:
                lhs, rhs = map(str.strip, line.split("->"))

                parsedExpressions[ replaceGreekSymbols(lhs) ] = parse_mathematica( replaceGreekSymbols(rhs) )

            except ValueError:

                print(f"Error parsing line: {line}")

        return parsedExpressions



    ## argumentSymbols specifies which symbols the function depends on.
    # In principle you can include parameters that the function actually doesn't depend on
    def _exprToFunction(self, parsedExpression, argumentSymbols: list[str]) -> Callable:

        ## Find sympy symbols corresponding to user str input
        sympySymbols = []
        for s in argumentSymbols:
            sympySymbols.append(self.symbols[s])

        return sympy.lambdify(sympySymbols, parsedExpression, modules='numpy')