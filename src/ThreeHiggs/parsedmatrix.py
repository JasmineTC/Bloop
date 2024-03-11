from .ParsedExpression import ParsedExpressionSystem

import sympy
import numpy as np
from typing import Callable, Tuple
from dataclasses import dataclass

@dataclass(slots=True)
class MatrixDefinitionFiles:
    matrixFile: str
    expressionsFile: str


class ParsedMatrix:
    """2D numpy array with symbolic coefficients.
    """

    def __init__(self, matrixFile: str, matrixElementDefinitionsFile: str = None):
        """Need 2 files: the matrix is read from matrixFile, definitions for symbols in the matrix
        are read from matrixElementDefinitionsFile.
        If definitions are not given, this becomes a direct function of its symbolic matrix elements.
        """
        self.bHasExpressions = False
        self.readFromFile(matrixFile, matrixElementDefinitionsFile)


    def readFromFile(self, matrixFile: str, matrixElementDefinitionsFile: str = None) -> None:
        """"""
        sympyMatrix: sympy.Matrix = self.parseMatrix(matrixFile)

        assert len(sympyMatrix.shape) == 2, "Tried to parse a non-2D object as ParsedMatrix!"

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


    def evaluate(self, arguments: list[str]) -> np.ndarray:
        """Evaluates the matrix element expressions and puts them in a 2D np.ndarray.
        Input list needs to be ordered as assumed by the expressions.
        """
        if (self.bHasExpressions):
            evaluatedExpressions = self.matrixElementExpressions.evaluateSystem(arguments)
            return self.matrix(evaluatedExpressions)
        else:
            raise RuntimeError("TODO! matrix mess...")
    
    def evaluateWithDict(self, arguments: dict[str, float]) -> np.ndarray:
        """Evaluates the matrix element expressions and puts them in a 2D np.ndarray.
        The input dict needs to contain keys for all function arguments needed by the expressions. 
        """
        if (self.bHasExpressions):
            evaluatedExpressions = self.matrixElementExpressions.evaluateSystemWithDict(arguments)
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

