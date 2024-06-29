from math import pi, log, sqrt
EulerGamma = 0.5772156649015329
Glaisher = 1.28242712910062

class ParsedExpression:
    def __init__(self, line, bReplaceGreekSymbols=True):
        from ThreeHiggs.MathematicaParsers import parseExpression
        parsedExpression = parseExpression(line)

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
    def __init__(self, fileName = None):
        self.parsedExpressions = list(map(lambda line: ParsedExpression(line, bReplaceGreekSymbols=True),
                                          open(fileName, "r", encoding="utf-8").readlines()))

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
    def __init__(self, matrixFileName, definitionsFileName):
        self.matrixElementExpressions = ParsedExpressionSystem(definitionsFileName)
        
        from ThreeHiggs.MathematicaParsers import parseMassMatrix
        self.matrix = compile(str(parseMassMatrix(open(matrixFileName, 'r', encoding = "utf-8").readlines())["matrix"]),
                              "",
                              mode = "eval")

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
    def __init__(self, fileName):
        from ThreeHiggs.MathematicaParsers import parseRotationMatrix
        self.symbolMap = parseRotationMatrix(open(fileName, 'r', encoding = "utf-8").readlines())["matrix"]

    def __call__(self, numericalM):
        """Evaluates our symbols by plugging in numbers from the input numerical matrix.
        Returns a dict with symbols names as keys.
        """

        return {symbol: numericalM[indices] for symbol, indices in self.symbolMap.items()}

