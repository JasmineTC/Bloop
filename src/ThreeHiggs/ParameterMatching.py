import numpy as np
import sympy
from sympy import pi
from sympy.parsing.mathematica import parse_mathematica

from typing import Callable, Tuple

from CommonUtils import getPackagedDataPath
from CommonUtils import replaceGreekSymbols

""" Notes.
Suppose you give DRalgo a model with a parameter called lam123. By default the parameter name after integrating out hard modes is then lam1233d, 
ie. DRalgo appends a '3d'. I don't want the '3d' part, but this is easy to remove after reading the file in a dict.  
"""


## This would make a good abstract class
class ParameterMatching:


    ## stores lambdas for evaluating params in the soft scale theory
    matchingRelations: dict[str, Callable]

    ## This list specifies names for parameters that enter the matching relations. 
    # Needs to match symbol names in the file that defines matching relations except that:
    # 1. "Unicode" symbols like λ are automatically converted to ANSI according to rules in CommonUtils.replaceGreekSymbols(). For ex. λ -> lam  
    # 2. You can include extra symbols that do not appear explicitly in the matching relations. Eg. "RGScale" 
    #parameterNames = [ 'T', 'Lb', 'Lf', 'g1', 'g2', 'g3',  ]  
    ## LN: not used currently since the symbols are read automatically from .txt


    def __init__(self):
        
        self.matchingRelations = {}
        


    def getMatchedParams(self, inputParams: dict[str, float]) -> dict[str, float]:

        ## put input params in an array so that matching relation lambdas can be evaluated
        inParamList = self._paramDictToOrderedList(inputParams)

        matchedParams = {}

        for key, value in self.matchingRelations.items():
            matchedParams[key] = value(*inParamList) ## Unpack because the lambdas don't take lists

        return matchedParams


    def createMatchingRelations(self, fileToRead: str = None) -> None:

        self.matchingRelations = {}

        ## Read from file if a filename was given
        if (fileToRead != None):
            symbolNames, parsedMatchingRelations = self.parseMatchingRelations(fileToRead)

            print("List of symbols in the matching relations:\n", symbolNames)
            self.__defineSymbols(symbolNames)
            
            ## I'll just have all matching relations take the same list of args
            self.parameterNames = symbolNames
            self.matchingRelations = self.lambdifyParsedExpressions(parsedMatchingRelations, self.parameterNames)


        ## No input file so use the manual matching implementation
        else: 
            self.implementMatchingRelations()
    ##

    def implementMatchingRelations(self):
        ### Manually implement matching here if you don't want to read them from file
        raise NotImplementedError
        """
        self.matchingRelations["g1sq"] = self._exprToFunction(self.__parsedMatchingRelations['g13d^2'], argumentSymbols)
        self.matchingRelations["g2sq"] = self._exprToFunction(self.__parsedMatchingRelations['g23d^2'], argumentSymbols)
        self.matchingRelations["g3sq"] = self._exprToFunction(self.__parsedMatchingRelations['g33d^2'], argumentSymbols)

        self.matchingRelations["lam11"] = self._exprToFunction(self.__parsedMatchingRelations["lam113d"], argumentSymbols)
        self.matchingRelations["lam22"] = self._exprToFunction(self.__parsedMatchingRelations["lam223d"], argumentSymbols)
        self.matchingRelations["lam33"] = self._exprToFunction(self.__parsedMatchingRelations["lam333d"], argumentSymbols)
        # etc
        """


    ## Convert a dict of input params to list that has same ordering as that defined in defineMatchingRelations().
    ## This is needed because the functions produced by sympy lambdify do not take dicts as argument 
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


    def parseMatchingRelations(self, filePath: str) -> Tuple[list[str], dict[str, sympy.Expr]]:
        
        ## Use this if changing this code to a proper package:
        #filePath = getPackagedDataPath("ThreeHiggs.Data", "softScaleParams.txt")
        
        #filePath = "Data/softScaleParams.txt"

        with open(filePath, "r") as file:
            expressions = file.readlines()

        parsedExpressions = {}
        parsedSymbols = [] ## Automatically find all symbols that appear in matching relations

        for line in expressions:
            try:
                lhs, rhs = map(str.strip, line.split("->"))

                ## Parsed right-hand side with Greek symbols replaced with ANSI text
                expr = parse_mathematica( replaceGreekSymbols(rhs) )

                parsedExpressions[ replaceGreekSymbols(lhs) ] = expr

                ## find symbols but store as string, not the sympy type  
                for symbol in expr.free_symbols:
                    ## NOTE this conversion will cause issues with pathological symbols like "A B" 
                    # https://stackoverflow.com/questions/59401738/convert-sympy-symbol-to-string-such-that-it-can-always-be-parsed
                    symbol_str = str(symbol)
                    if symbol_str not in parsedSymbols:
                        parsedSymbols.append(symbol_str)
                        
            except ValueError:

                print(f"Error parsing line: {line}")

        parsedSymbols.sort()
        return parsedSymbols, parsedExpressions
    

    ## Convert dict of sympy expressions to dict of lambdas, same keys. argumentNames (list of strings) 
    ## specifies what arguments the lambda function needs 
    def lambdifyParsedExpressions(self, parsedExpressions: dict, argumentNames: list[str]) -> dict[Callable]:

        convertedExpressions = {}
        
        for key, expr in parsedExpressions.items():
            convertedExpressions[key] = self._exprToFunction(expr, argumentNames)

        return convertedExpressions



    ## argumentSymbols specifies which symbols the function depends on. Give strings, we convert internally to sympy symbols.
    # In principle you can include parameters that the function actually doesn't depend on
    def _exprToFunction(self, parsedExpression, argumentSymbols: list[str]) -> Callable:

        ## Find sympy symbols corresponding to user str input
        sympySymbols = []
        for s in argumentSymbols:
            sympySymbols.append(self.symbols[s])

        return sympy.lambdify(sympySymbols, parsedExpression, modules='numpy')