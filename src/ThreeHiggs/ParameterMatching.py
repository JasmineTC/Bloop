import numpy as np

from typing import Callable, Tuple

from .ParsedExpression import ParsedExpression
from .CommonUtils import replaceGreekSymbols

""" Notes.
Suppose you give DRalgo a model with a parameter called lam123. By default the parameter name after integrating out hard modes is then lam1233d, 
ie. DRalgo appends a '3d'. I don't want the '3d' part, but this is easy to remove after reading the file in a dict.  
"""

## TODO replace this with ParsedExpressionSystem. It does almost the same things already.

## This would make a good abstract class
class ParameterMatching:
    ## This list specifies names for parameters that enter the matching relations. 
    # Needs to match symbol names in the file that defines matching relations except that:
    # 1. "Unicode" symbols like λ are automatically converted to ANSI according to rules in CommonUtils.replaceGreekSymbols(). For example λ -> lam  
    # 2. You can include extra symbols that do not appear explicitly in the matching relations. Eg. "RGScale" 
    #parameterNames = [ 'T', 'Lb', 'Lf', 'g1', 'g2', 'g3',  ]  
    ## LN: not used currently since the symbols are read automatically from .txt

    def __init__(self):
        self.matchingRelations = {}

    def getMatchedParams(self, inputParams):
        return {key: expr(inputParams) for key, expr in self.matchingRelations.items()}

    def __call__(self, inputParams):
        return self.getMatchedParams(inputParams)

    def createMatchingRelations(self, fileToRead):
        self.parameterNames, self.matchingRelations = self.parseMatchingRelations(fileToRead)

    def parseMatchingRelations(self, filePath: str) -> Tuple[list[str], dict[str, ParsedExpression]]:

        with open(filePath, "r", encoding="utf-8") as file:
            expressions = file.readlines()

        ## Dict for storing ParsedExpression objects
        parsedExpressions = {}
        parsedSymbols = [] ## Automatically find all symbols that appear in matching relations

        for line in expressions:
            lhs, rhs = map(str.strip, line.split("->"))

            expr = ParsedExpression(rhs, bReplaceGreekSymbols=True)

            ## Replace Greeks also on the LHS and store in dict as LHS : parsedRHS
            parsedExpressions[ replaceGreekSymbols(lhs) ] = expr

            ## find symbols but store as string, not the sympy type  
            for symbol in expr.symbols:
                ## NOTE this conversion will cause issues with pathological symbols like "A B" 
                # https://stackoverflow.com/questions/59401738/convert-sympy-symbol-to-string-such-that-it-can-always-be-parsed
                symbol_str = str(symbol)
                if symbol_str not in parsedSymbols:
                    parsedSymbols.append(symbol_str)

        parsedSymbols.sort()
        return parsedSymbols, parsedExpressions
    

    def lambdifyMatchingRelations(self, parsedExpressions: dict[str, ParsedExpression], argumentNames: list[str]) -> dict[str, Callable]:
        """Convert dict of ParsedExpressions to a dict of lambdas, same keys. argumentNames (list of strings)
        specifies what arguments the lambda functions need 
        """

        convertedExpressions = {}
        
        for key, expr in parsedExpressions.items():
            convertedExpressions[key] = expr

        return convertedExpressions
