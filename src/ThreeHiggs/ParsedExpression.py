import numpy as np
import sympy
from sympy.parsing.mathematica import parse_mathematica

from typing import Callable, Tuple

from .CommonUtils import replaceGreekSymbols

from math import log, sqrt
pi = float(sympy.pi.evalf())
EulerGamma = float(sympy.EulerGamma.evalf())
Glaisher = 1.28242712910062

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

    def __call__(self, inputDict: dict[str, float], bReturnDict=False) -> list[float]:
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
        return [expr.lambdaExpression for expr in self.parsedExpressions]

class ParsedMatrix:
    def __init__(self, matrixFile, matrixElementDefinitionsFile = None):
        """Need 2 files: the matrix is read from matrixFile, definitions for symbols in the matrix
        are read from matrixElementDefinitionsFile.
        If definitions are not given, this becomes a direct function of its symbolic matrix elements.
        """
        self.bHasExpressions = False
        sympyMatrix = self.parseMatrix(matrixFile)
        self.shape = sympyMatrix.shape

        if (matrixElementDefinitionsFile):
            self.matrixElementExpressions = ParsedExpressionSystem(matrixElementDefinitionsFile)
            self.bHasExpressions = True

        ## Make index map so that we know which symbol is which a_ij
        self.symbolMap: dict[str, Tuple[int, int]]  = {}
        for i in range(sympyMatrix.shape[0]):
            for j in range(sympyMatrix.shape[1]):
                element = sympyMatrix[i, j]

                if element.is_symbol:
                    self.symbolMap[str(element)] = i, j

        
        """Lambdify the matrix. For arguments: 
            1) If expressions were given, use same order as how the expressions are output
            2) If expressions were not given, use whatever order 
        """
        if (self.bHasExpressions):
            self.argumentSymbols = self.matrixElementExpressions.getExpressionNames()
        else:
            ## TODO use set
            self.argumentSymbols = [str(s) for s in sympyMatrix.free_symbols]

        self.matrix = sympy.lambdify(self.argumentSymbols, sympyMatrix, modules='numpy')


    @staticmethod
    def parseMatrix(fileName: str) -> sympy.Matrix:
        """Parses a symbolic matrix from text file.
        The matrix is assumed to be in format:
        {a11, a12, ...}
        {a21, a22, ...}
        etc, ie. one row per line. This is how Mathematica exports a 2D table.
        """

        with open(fileName, 'r') as file:
            lines = file.readlines()

        ## Interpret the lines
        matrixData = [[symbol.strip() for symbol in line.strip().lstrip('{').rstrip('}').split(',')] for line in lines]

        symbolicMatrix = sympy.Matrix(matrixData)

        return symbolicMatrix

    def __call__(self, arguments):
        """Evaluates the matrix element expressions and puts them in a 2D np.ndarray.
        The input dict needs to contain keys for all function arguments needed by the expressions. 
        """
        if (self.bHasExpressions):
            evaluatedExpressions = self.matrixElementExpressions(arguments)
            ## the above gives a list but we need a tuple, so unpack
            return self.matrix(*evaluatedExpressions) 
        else:
            raise RuntimeError("TODO! matrix mess...")


    def matchSymbols(self, numericalMatrix: np.ndarray) -> dict[str, float]:
        """Evaluates our symbols by plugging in numbers from the input numerical matrix.
        Returns a dict with symbols names as keys.
        """
        assert numericalMatrix.shape == self.shape, "Matrix substitute error, shapes don't match!"

        outDict = {}
        for symbol in self.argumentSymbols:
            idx = self.symbolMap[symbol]
            outDict[symbol] = numericalMatrix[idx]

        return outDict


    @staticmethod
    def parseConstantMatrix(fileName: str) -> np.ndarray:
        sympyMatrix = ParsedMatrix.parseMatrix(fileName)

        # should not have any free symbols
        assert len(sympyMatrix.free_symbols) == 0, "Failed to parse constant matrix: has free symbols!"
        # can't lambdify with no arguments, so do standard conversion to numpy array instead
        return np.array(sympyMatrix.tolist()).astype(np.float64)

