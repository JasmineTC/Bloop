from cmath import pi, log, sqrt
from numpy import euler_gamma
Glaisher = 1.28242712910062

class ParsedExpression:
    def __init__(self, parsedExpression):
        self.identifier = parsedExpression["identifier"]
        self.expression = parsedExpression["expression"]
        self.symbols = parsedExpression["symbols"]

        self.lambdaExpression = compile(self.expression, "<string>", mode = "eval")

    def evaluate(self, functionArguments: list[float]) -> float:
        return eval(self.lambdaExpression, 
                    functionArguments | {"log": log, 
                                         "sqrt": sqrt, 
                                         "pi": pi, 
                                         "EulerGamma": euler_gamma,
                                         "Glaisher": Glaisher})

""" class ParsedExpressionSystem -- Describes a collection of ParsedExpressions that are to be evaluated simultaneously with same input.
"""
class ParsedExpressionSystem:
    def __init__(self, parsedExpressionSystem):
        self.parsedExpressions = [ParsedExpression(parsedExpression) 
                                  for parsedExpression in parsedExpressionSystem]

    def evaluate(self, inputDict: dict[str, float], bReturnDict=False) -> list[float]:
        """Optional argument is a hack"""
        ## Collect inputs from the dict and put them in correct order. I do this by taking the right order from our first expression.
        ## This is fine since all our expressions use the same input list. 
        outList = [None] * len(self.parsedExpressions)    
        for i in range(len(outList)):
            outList[i] = self.parsedExpressions[i].evaluate(inputDict)

        if not bReturnDict:
            return outList
        else:
            return  { self.parsedExpressions[i].identifier : outList[i] for i in range(len(outList)) } 

    def getExpressionNames(self) -> list[str]:
        return [ expr.identifier for expr in self.parsedExpressions ]

class ParsedExpressionArray:
    def __init__(self, parsedExpression):
        self.identifier = parsedExpression["identifier"]
        self.expression = parsedExpression["expression"]
        self.symbols = parsedExpression["symbols"]

        self.lambdaExpression = compile(self.expression, "<string>", mode = "eval")

    def evaluate(self, params):
        return eval(self.lambdaExpression,  {"log": log, 
                                             "sqrt": sqrt, 
                                             "pi": pi, 
                                             "EulerGamma": euler_gamma,
                                             "Glaisher": Glaisher,
                                             "params": params})

class ParsedExpressionSystemArray:
    def __init__(self, parsedExpressionSystem):
        self.parsedExpressions = [ParsedExpressionArray(parsedExpression) 
                                  for parsedExpression in parsedExpressionSystem]

    def evaluate(self, params):
        return [expression.evaluate(params) for expression in self.parsedExpressions]

    def getExpressionNames(self) -> list[str]:
        return [ expr.identifier for expr in self.parsedExpressions ]

class MassMatrix:
    def __init__(self, massMatrix):
        self.definitions = ParsedExpressionSystem(massMatrix["definitions"])
        self.matrix = compile(massMatrix["matrix"], "<string>", mode = "eval")

    def evaluate(self, arguments):
        """Evaluates the matrix element expressions and puts them in a 2D np.ndarray.
        The input dict needs to contain keys for all function arguments needed by the expressions. 
        """
        arguments |= self.definitions.evaluate(arguments, bReturnDict = True)
        return eval(self.matrix, arguments | {"log": log, 
                                              "sqrt": sqrt, 
                                              "pi": pi, 
                                              "EulerGamma": euler_gamma,
                                              "Glaisher": Glaisher})

class RotationMatrix:
    def __init__(self, symbolMap):
        self.symbolMap = symbolMap["matrix"]

    def evaluate(self, numericalM):
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

        self.assertEqual(reference, ParsedExpression(source).evaluate({"lam": 100, "mssq": 100}))

    def test_ParsedExpressionComplex(self):
        source = {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                  "identifier": "Identifier",
                  "symbols": ['lam', 'mssq']}

        reference = complex(5.826048814042759, 1.1475471676948477)

        self.assertEqual(reference, 
                         ParsedExpression(source).evaluate({"lam": complex(100, 100), 
                                                            "mssq": complex(100, 100)}))

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

        self.assertEqual(reference, 
                         ParsedExpressionSystem(source).evaluate({"lam": 100, "mssq": 100}))

    def test_MassMatrix(self):
        source = {"definitions": [{"expression": "1",
                                   "identifier":"mssq",
                                   "symbols": []}],
                  "matrix": "[[1, 0], [0, mssq]]"}

        reference = [[1, 0], [0, 1]]
        self.assertEqual(reference, MassMatrix(source).evaluate({}))

    def test_RotationMatrix(self):
        source = {"matrix": {"mssq00": [0, 0], "mssq11": [1, 1]}}
        reference = {"mssq00": 1, "mssq11": -1}

        self.assertEqual(reference, RotationMatrix(source).evaluate([[1, 0], [0, -1]]))

