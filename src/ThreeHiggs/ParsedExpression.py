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
        self.sympyExpression, self.symbols = self._parseMathematicaExpression(expression)
        self.lambdaExpression = compile(str(self.sympyExpression), "<string>", mode = "eval")

    def __call__(self, functionArguments: list[float]) -> float:
        return eval(self.lambdaExpression, 
                    functionArguments | {"log": log, 
                                         "sqrt": sqrt, 
                                         "pi": pi, 
                                         "EulerGamma": EulerGamma,
                                         "Glaisher": Glaisher})

    def __str__(self) -> str:
        return self.identifier + " == " + self.stringExpression

    @staticmethod
    def _parseMathematicaExpression(expression: str, bSubstituteConstants: bool = True) -> Tuple[sympy.Expr, list[sympy.Symbol]]:
        """ Parses a string exported from Mathematica.
        bSubstituteConstants specifies if we replace numerical constants like Glaisher with their numerical values
        -----------
        Returns a tuple of sympy objects:
        expression, symbols.
        """

        sympyExpr = parse_mathematica(expression)
        symbols = sympyExpr.free_symbols
        return sympyExpr, symbols

""" class ParsedExpressionSystem -- Describes a collection of ParsedExpressions that are to be evaluated simultaneously with same input.
"""
class ParsedExpressionSystem:

    parsedExpressions: list[ParsedExpression]
    """Arguments needed to evaluate our lambda functions, in correct order."""
    functionArguments: list[str]

    def __init__(self, fileName: str = None):
        if (fileName):
            self.parseExpressions(fileName)

    def parseExpressions(self, fileName: str) -> Tuple:
        """Read system of expressions from file, one per line.
        The expressions will become functions of all symbols that appear. 
        """

        parsedSymbols = []
        self.parsedExpressions = []
        
        with open(fileName, "r", encoding="utf-8") as file:
            for line in file.readlines():
                expr = ParsedExpression(line, bReplaceGreekSymbols=True)
    
                self.parsedExpressions.append(expr)
    
                ## find symbols but store as string, not the sympy type  
                for symbol in expr.symbols:
                    ## NOTE this conversion will cause issues with pathological symbols like "A B" 
                    # https://stackoverflow.com/questions/59401738/convert-sympy-symbol-to-string-such-that-it-can-always-be-parsed
                    symbol_str = str(symbol)
                    if symbol_str not in parsedSymbols:
                        parsedSymbols.append(symbol_str)

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

    def __str__(self) -> str:
        return str([ str(expr) for expr in self.parsedExpressions ])

    def getExpressionNames(self) -> list[str]:
        return [ expr.identifier for expr in self.parsedExpressions ]

"""Class SystemOfEquations -- System of parsed expression that we interpret as a set of equation. 
Each expression is interpreted as an equation of form ``expr == 0``. We also distinguish between symbols 
that describe the unknowns versus symbols that are known inputs to the expressions.
"""
class SystemOfEquations(ParsedExpressionSystem):

    def __init__(self, fileName: str, unknownVariables: list[str]):
        
        super().__init__(fileName)

        ## what we solve for
        self.unknownVariables = unknownVariables

        ## Check that we have same number of eqs as unknowns 
        assert len(unknownVariables) == len(self.parsedExpressions)

        ## Check that the eqs actually contain our unknowns
        for s in unknownVariables:
            if not s in self.functionArguments:
                errorString = (f"Error parsing equations from file {fileName}! The expressions don't seem to depend on variable {s}, "
                    f"which was listed in our unknownParameters: {unknownVariables}.")
                raise RuntimeError(errorString)
            

        self._rearrangeSymbols()

    def _rearrangeSymbols(self) -> None:
        """This rearranges our internal symbol list and expression lambdas so that the unknown variables come before other variables.
        For internal use only.
        """
        ## First create a filtered list that contains all non-unknowns. With a crazy oneliner of course
        filteredArguments = [ item for item in self.functionArguments if item not in self.unknownVariables ]
        rearrangedArguments = self.unknownVariables + filteredArguments

        ## "known" inputs
        self.otherVariables = filteredArguments

    def getEquations(self) -> list[Callable]:
        """Returns a list of our lambdas.
        """
        eqs = [None] * len(self.parsedExpressions)
        for i in range(len(eqs)):
            eqs[i] = self.parsedExpressions[i].lambdaExpression

        return eqs
        
