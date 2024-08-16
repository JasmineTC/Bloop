from numpy import pi, log, sqrt
EulerGamma = 0.5772156649015329
Glaisher = 1.28242712910062

class ParsedExpression:
    def __init__(self, parsedExpression):
        self.identifier = parsedExpression["identifier"]
        self.expression = parsedExpression["expression"]
        self.symbols = parsedExpression["symbols"]

        self.lambdaExpression = compile(self.expression, "<string>", mode = "eval")

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
    def __init__(self, parsedExpressionSystem):
        self.parsedExpressions = [ParsedExpression(parsedExpression) 
                                  for parsedExpression in parsedExpressionSystem]

    def __call__(self, inputDict: dict[str, float], bReturnDict=False) -> list[float]:
        """Optional argument is a hack
        """
        ## Collect inputs from the dict and put them in correct order. I do this by taking the right order from our first expression.
        ## This is fine since all our expressions use the same input list. 
        outList = [None] * len(self.parsedExpressions)
        for i in range(len(outList)):
            outList[i] = self.parsedExpressions[i](inputDict)

        if not bReturnDict:
            return outList
        else:
            return  { self.parsedExpressions[i].identifier : outList[i] for i in range(len(outList)) } 

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

class MassMatrix:
    def __init__(self, matrix, definitions):
        self.matrixElementExpressions = definitions
        self.matrix = matrix 

    def __call__(self, arguments):
        """Evaluates the matrix element expressions and puts them in a 2D np.ndarray.
        The input dict needs to contain keys for all function arguments needed by the expressions. 
        """
        arguments |= self.matrixElementExpressions(arguments, bReturnDict = True)
        return eval(self.matrix, arguments | {"log": log, 
                                              "sqrt": sqrt, 
                                              "pi": pi, 
                                              "EulerGamma": EulerGamma,
                                              "Glaisher": Glaisher})

class RotationMatrix:
    def __init__(self, symbolMap):
        self.symbolMap = symbolMap["matrix"]

    def __call__(self, numericalM):
        """Evaluates our symbols by plugging in numbers from the input numerical matrix.
        Returns a dict with symbols names as keys.
        """

        return {symbol: numericalM[indices[0]][indices[1]] for symbol, indices in self.symbolMap.items()}

from unittest import TestCase
class ParsedExpressionUnitTests(TestCase):
    def test_ParsedExpression(self):
        source = {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                  "identifier": "Identifier",
                  "symbols": ['lam', 'mssq']}

        reference = 5.400944901447568

        self.assertEqual(reference, ParsedExpression(source)({"lam": 100, "mssq": 100}))

    def test_ParsedExpressionSystem(self):
        source = [{"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                   "identifier": "Identifier",
                   "symbols": ['lam', 'mssq']},
                  {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                   "identifier": "Identifier",
                   "symbols": ['lam', 'mssq']},
                  {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                   "identifier": "Identifier",
                   "symbols": ['lam', 'mssq']}]

        reference = [5.400944901447568, 5.400944901447568, 5.400944901447568]

        self.assertEqual(reference, ParsedExpressionSystem(source)({"lam": 100, "mssq": 100}))

    def test_MassMatrix(self):
        source = [{"matrix": "[[1, 0], [0, mssq]]"}["matrix"],
                  ParsedExpressionSystem([{"identifier": "mssq", "symbols": [], "expression": "1"}])]

        reference = [[1, 0], [0, 1]]
        self.assertEqual(reference, MassMatrix(*source)({}))

    def test_RotationMatrix(self):
        source = {"matrix": {"mssq00": [0, 0], "mssq11": [1, 1]}}
        reference = {"mssq00": 1, "mssq11": -1}

        self.assertEqual(reference, RotationMatrix(source)([[1, 0], [0, -1]]))

