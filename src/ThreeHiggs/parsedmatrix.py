from .ParsedExpression import ParsedExpressionSystem

import sympy
import numpy as np
from typing import Callable, Tuple

class ParsedMatrix:
    """2D numpy array with symbolic coefficients.
    """


    def __init__(self, matrixFile: str, matrixElementDefinitionsFile: str):
        """Need 2 files: the matrix is read from matrixFile, definitions for symbols in the matrix
        are read from matrixElementDefinitionsFile.
        """
        self.readFromFile(matrixFile, matrixElementDefinitionsFile)


    def readFromFile(self, matrixFile: str, matrixElementDefinitionsFile: str) -> None:
        """"""
        self.matrixElementExpressions = ParsedExpressionSystem(matrixElementDefinitionsFile)
        self.matrix, self.functionArguments = self.parseMatrix(matrixFile)


    def parseMatrix(self, fileName: str) -> Tuple[Callable, list[str]] :
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

        ## Free sympy symbols
        #symbols = symbolicMatrix.free_symbols

        ## String valued names for the lambda function arguments. Use same order as how
        ## self.matrixElementExpressions outputs its expressions
        argumentSymbols = self.matrixElementExpressions.getExpressionNames()

        # Lambdify this to a 2D numpy array
        matrixLambda = sympy.lambdify(argumentSymbols, symbolicMatrix, modules='numpy')

        return matrixLambda, argumentSymbols


    def evaluate(self, arguments: list[str]) -> np.ndarray:
        """Evaluates the matrix element expressions and puts them in a 2D np.ndarray.
        Input list needs to be ordered as assumed by the expressions.
        """
        evaluatedExpressions = self.matrixElementExpressions.evaluateSystem(arguments)
        return self.matrix(evaluatedExpressions)
    
    def evaluateWithDict(self, arguments: dict[str, float]) -> np.ndarray:
        """Evaluates the matrix element expressions and puts them in a 2D np.ndarray.
        The input dict needs to contain keys for all function arguments needed by the expressions. 
        """
        evaluatedExpressions = self.matrixElementExpressions.evaluateSystemWithDict(arguments)
        ## the above gives a list but we need a tuple, so unpack
        return self.matrix(*evaluatedExpressions)