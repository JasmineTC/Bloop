import numpy as np
import sympy
from sympy.parsing.mathematica import parse_mathematica

from typing import Callable, Tuple

from .CommonUtils import replaceGreekSymbols

pi = float(sympy.pi.evalf())
EulerGamma = float(sympy.EulerGamma.evalf())
Glaisher = 1.28242712910062
def log(x):
    return float(sympy.log(float(x)))

def sqrt(x):
    return float(sympy.sqrt(float(x)))

class ParsedExpression:
    def __init__(self, expression: str = None, bReplaceGreekSymbols=True):
        expression = replaceGreekSymbols(expression) 
        self.stringExpression = expression

        if ("->" in expression):
            lhs, rhs = map(str.strip, expression.split("->"))

            expression = rhs
            self.identifier = lhs

        self.stringExpression = expression
        self.sympyExpression = parse_mathematica(expression)
        self.symbols = self.sympyExpression.free_symbols
        self.lambdaExpression = compile(str(self.sympyExpression), "<string>", mode = "eval")

    def __call__(self, functionArguments: list[float]) -> float:
        return eval(self.lambdaExpression, 
                    functionArguments | {"log": log, 
                                         "sqrt": sqrt, 
                                         "pi": pi, 
                                         "EulerGamma": EulerGamma,
                                         "Glaisher": Glaisher})

""" class ParsedExpressionSystem -- Describes a collection of ParsedExpressions that are to be evaluated simultaneously with same input.
"""
class ParsedExpressionSystem:
    def __init__(self, fileName = None):
        self.parsedExpressions = list(map(lambda line: ParsedExpression(line, bReplaceGreekSymbols=True),
                                          open(fileName, "r", encoding="utf-8").readlines()))

    def evaluateSystemWithDict(self, inputDict: dict[str, float], bReturnDict=False) -> list[float]:
        """Optional argument is a hack
        """
        ## Collect inputs from the dict and put them in correct order. I do this by taking the right order from our first expression.
        ## This is fine since all our expressions use the same input list. 
        outList = [None] * len(self.parsedExpressions)
        for i in np.arange(len(outList)):
            outList[i] = self.parsedExpressions[i](inputDict)

        if not bReturnDict:
            return outList
        else:
            return  { self.parsedExpressions[i].identifier : outList[i] for i in np.arange(len(outList)) } 

    def __call__(self, arguments: list[float]) -> Tuple:
        """Just calls evaluateSystem()"""
        return self.evaluateSystemWithDict(arguments)

    def getExpressionNames(self) -> list[str]:
        return [ expr.identifier for expr in self.parsedExpressions ]

"""Class SystemOfEquations -- System of parsed expression that we interpret as a set of equation. 
Each expression is interpreted as an equation of form ``expr == 0``. We also distinguish between symbols 
that describe the unknowns versus symbols that are known inputs to the expressions.
"""
class SystemOfEquations(ParsedExpressionSystem):
    def __init__(self, fileName, unknownVariables):
        super().__init__(fileName)

        ## what we solve for
        self.unknownVariables = unknownVariables

        filteredArguments = [item for item in self.functionArguments if item not in self.unknownVariables]
        rearrangedArguments = self.unknownVariables + filteredArguments

        ## "known" inputs
        self.otherVariables = filteredArguments

    def getEquations(self):
        eqs = [None] * len(self.parsedExpressions)
        for i in range(len(eqs)):
            eqs[i] = self.parsedExpressions[i].lambdaExpression

        return eqs
        
