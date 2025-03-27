from cmath import pi, log, sqrt
from numpy import euler_gamma
Glaisher = 1.28242712910062

import copy

class ParsedExpression:
    def __init__(self, parsedExpression, fileName):
        self.identifier = parsedExpression["identifier"]
        self.expression = parsedExpression["expression"]
        self.symbols = parsedExpression["symbols"]
        self.fileName = fileName

        self.lambdaExpression = compile(self.expression, "<string>", mode = "eval")

    def evaluate(self, functionArguments: list[float]) -> float:
        return eval(self.lambdaExpression, 
                    functionArguments | {"log": log, 
                                         "sqrt": sqrt, 
                                         "pi": pi, 
                                         "EulerGamma": euler_gamma,
                                         "Glaisher": Glaisher})

class ParsedExpressionSystem:
    def __init__(self, parsedExpressionSystem, fileName):
        self.parsedExpressions = [ParsedExpression(parsedExpression, fileName) 
                                  for parsedExpression in parsedExpressionSystem]
        self.fileName = fileName
        

    def evaluate(self, inputDict: dict[str, float], bReturnDict=False) -> list[float]:
        """Optional argument is a hack"""
        outList = [expression.evaluate(inputDict) for expression in self.parsedExpressions] 
        if bReturnDict:
            return { self.parsedExpressions[i].identifier : outList[i] for i in range(len(outList)) }
        return  outList 

    def getExpressionNames(self) -> list[str]:
        return [ expr.identifier for expr in self.parsedExpressions ]

class ParsedExpressionArray:
    def __init__(self, parsedExpression, fileName):
        self.identifier = parsedExpression["identifier"]
        self.expression = parsedExpression["expression"]
        self.symbols = parsedExpression["symbols"]
        self.fileName = fileName

        self.lambdaExpression = compile(self.expression, "<string>", mode = "eval")

    def evaluate(self, params):
        return eval(self.lambdaExpression,  {"log": log, 
                                             "sqrt": sqrt, 
                                             "pi": pi, 
                                             "EulerGamma": euler_gamma,
                                             "Glaisher": Glaisher,
                                             "params": params})

class ParsedExpressionSystemArray:
    def __init__(self, parsedExpressionSystem, allSymbols, fileName):
        self.parsedExpressions = [(allSymbols.index(parsedExpression["identifier"]), 
                                   ParsedExpressionArray(parsedExpression, fileName))
                                  for parsedExpression in parsedExpressionSystem]

        self.allSymbols = allSymbols
        self.fileName = fileName
        

    def evaluate(self, params):
        ## Look into using copy.replace 3.13 feature
        newParams = copy.deepcopy(params)
        
        for expression in self.parsedExpressions:
            newParams[expression[0]] = expression[1].evaluate(params)

        return newParams

    def getParamSubset(self, params):
        return [params[index] for index, symbol in enumerate(self.allSymbols) if symbol in self.getExpressionNames()]

    def getExpressionNames(self) -> list[str]:
        return [ expr[1].identifier for expr in self.parsedExpressions ]

class MassMatrix:
    def __init__(self, massMatrix, fileName):
        self.definitions = ParsedExpressionSystem(massMatrix["definitions"], fileName)
        self.matrix = compile(massMatrix["matrix"], "<string>", mode = "eval")
        self.fileName = fileName

    def evaluate(self, arguments):
        arguments |= self.definitions.evaluate(arguments, bReturnDict = True)
        return eval(self.matrix, arguments | {"log": log, 
                                              "sqrt": sqrt, 
                                              "pi": pi, 
                                              "EulerGamma": euler_gamma,
                                              "Glaisher": Glaisher})

class RotationMatrix:
    def __init__(self, symbolMap, fileName):
        self.symbolMap = symbolMap["matrix"]
        self.fileName = fileName

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

        self.assertEqual(reference, ParsedExpression(source, None).evaluate({"lam": 100, "mssq": 100}))

    def test_ParsedExpressionComplex(self):
        source = {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                  "identifier": "Identifier",
                  "symbols": ['lam', 'mssq']}

        reference = complex(5.826048814042759, 1.1475471676948477)

        self.assertEqual(reference, 
                         ParsedExpression(source, None).evaluate({"lam": complex(100, 100), 
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
                         ParsedExpressionSystem(source, None).evaluate({"lam": 100, "mssq": 100}))

    def test_MassMatrix(self):
        source = {"definitions": [{"expression": "1",
                                   "identifier":"mssq",
                                   "symbols": []}],
                  "matrix": "[[1, 0], [0, mssq]]"}

        reference = [[1, 0], [0, 1]]
        self.assertEqual(reference, MassMatrix(source, None).evaluate({}))

    def test_RotationMatrix(self):
        source = {"matrix": {"mssq00": [0, 0], "mssq11": [1, 1]}}
        reference = {"mssq00": 1, "mssq11": -1}

        self.assertEqual(reference, RotationMatrix(source, None).evaluate([[1, 0], [0, -1]]))

