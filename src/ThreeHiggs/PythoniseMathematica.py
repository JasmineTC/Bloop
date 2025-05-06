from ThreeHiggs.GetLines import getLines
from json import dump
from sympy import Matrix
from sympy.parsing.mathematica import parse_mathematica

def replaceGreekSymbols(string: str) -> str:
    ## TODO use unicodedata package here to do magic. 
    # Or bully DRalgo people to removing greek symbols
    return string.replace(u"\u03BB", "lam").replace(u"\u03BC", "mu")

def removeSuffices(string):
    return string.replace("^2", "sq")

def replaceSymbolsWithIndices(expression, symbols):
    expression = replaceGreekSymbols(expression)
    ## Reverse needed to deal with lam23 and lam23p i.e. substring replaces larger full string
    for idx, symbol in enumerate(sorted(symbols, reverse = True)):
        expression = expression.replace(symbol , f"params[{idx}]")

    return expression

def pythoniseExpressionArray(line, allSymbols):
    identifier, line = map(str.strip, line.split("->")) if ("->" in line) else ("anonymous", line)
    
    identifier = removeSuffices(replaceGreekSymbols(identifier))
    expression = parse_mathematica(replaceGreekSymbols(line))
    symbols = [str(symbol) for symbol in expression.free_symbols]

    return {"identifier": identifier, 
            "expression": replaceSymbolsWithIndices(str(expression), allSymbols), 
            "symbols": sorted(symbols)}

def pythoniseExpression(line):
    identifier, line = map(str.strip, line.split("->")) if ("->" in line) else ("anonymous", line)

    identifier = removeSuffices(replaceGreekSymbols(identifier))
    expression = parse_mathematica(replaceGreekSymbols(line))
    symbols = [str(symbol) for symbol in expression.free_symbols]

    return {"identifier": identifier, "expression": str(expression), "symbols": sorted(symbols)}

def pythoniseExpressionSystemArray(lines, allSymbols):
    return [pythoniseExpressionArray(line, allSymbols) for line in lines]

def pythoniseExpressionSystem(lines):
    return [pythoniseExpression(line) for line in lines]

def pythoniseMatrix(lines):
    return [[symbol.strip() for symbol in line.strip()
                                              .strip('}')
                                              .strip('{')
                                              .split(',')] for line in lines]

def pythoniseMassMatrix(definitionsLines, matrixLines):
    from sympy.core.sympify import sympify
    ##Need sympify that the parsed matrix line is ufuncable
    ##Need str to make compatable with json
    return {"definitions": pythoniseExpressionSystem(definitionsLines),
            "matrix": str(sympify(pythoniseMatrix(matrixLines)))}

def pythoniseRotationMatrix(lines):
    sympyMatrix = Matrix(pythoniseMatrix(lines))
    shape = sympyMatrix.shape
    
    symbolMap = {}
    for i in range(shape[0]):
        for j in range(shape[1]):
            element = sympyMatrix[i, j]

            if element.is_symbol:
                symbolMap[str(sympyMatrix[i, j])] = (i, j)
    return {"matrix": symbolMap}

def pythoniseMathematica(args):
    veffLines = getLines(args.loFile)
    veffLines += getLines(args.nloFile)
    if (args.loopOrder >= 2):
        veffLines += getLines(args.nnloFile)
    
    allSymbols = [replaceGreekSymbols(symbol) for symbol in getLines(args.allSymbolsFile, mode = "json")]
    ## Move get lines to the functions?
    ## Not ideal to have nested dicts but is future proof for when we move to arrays
    dump({"betaFunctions4D": {"expressions": pythoniseExpressionSystemArray(getLines(args.betaFunctions4DFile), allSymbols),
                              "fileName": args.betaFunctions4DFile},
          "hardToSoft": {"expressions":  pythoniseExpressionSystem(getLines(args.hardToSoftFile)),
                         "fileName": args.hardToSoftFile},
          "softScaleRGE": {"expressions": pythoniseExpressionSystem(getLines(args.softScaleRGEFile)),
                           "fileName": args.softScaleRGEFile},
          "softToUltraSoft": {"expressions": pythoniseExpressionSystem(getLines(args.softToUltraSoftFile)),
                              "fileName": args.softToUltraSoftFile},
          "vectorMassesSquared": {"expressions": pythoniseExpressionSystem(getLines(args.vectorMassesSquaredFile)),
                                  "fileName": args.vectorMassesSquaredFile},
          "vectorShortHands": {"expressions": pythoniseExpressionSystem(getLines(args.vectorShortHandsFile)),
                               "fileName": args.vectorShortHandsFile},
          "veff": {"expressions": pythoniseExpressionSystem(veffLines),
                     "fileName": "Combined Veff files"},
          "scalarMassMatrixUpperLeft": {"expressions": pythoniseMassMatrix(getLines(args.scalarMassMatrixUpperLeftDefinitionsFile),
                                                                       getLines(args.scalarMassMatrixUpperLeftFile)),
                                        "fileName": (args.scalarMassMatrixUpperLeftDefinitionsFile,
                                                     args.scalarMassMatrixBottomRightFile)},
          "scalarMassMatrixBottomRight": {"expressions": pythoniseMassMatrix(getLines(args.scalarMassMatrixBottomRightDefinitionsFile),
                                                                         getLines(args.scalarMassMatrixBottomRightFile)),
                                        "fileName": (args.scalarMassMatrixBottomRightDefinitionsFile,
                                                     args.scalarMassMatrixBottomRightFile)},
          "scalarRotationMatrix": {"expressions": pythoniseRotationMatrix(getLines(args.scalarRotationFile)),
                                   "fileName": args.scalarRotationFile},
          "scalarPermutationMatrix": getLines(args.scalarPermutationFile, mode="json")},
         open(args.parsedExpressionsFile, "w"),
         indent = 4)

from unittest import TestCase
class PythoniseMathematicaUnitTests(TestCase):
    def test_replaceGreekSymbols(self):
        reference = ["lam", "lam lam", "mu", "mu mu", "lam mu", "mu lam"]
        source = ["λ", "λ λ", "μ", "μ μ", "λ μ", "μ λ"]

        self.assertEqual(reference, [replaceGreekSymbols(sourceString) for sourceString in source])

    def test_removeSuffices(self):
        reference = ["myVarsq", "sqmyVar"]

        source = ["myVar^2", "^2myVar"]

        self.assertEqual(reference, [removeSuffices(sourceString) for sourceString in source])

    def test_pythoniseExpression(self):
        reference = {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                     "identifier": "Identifier",
                     "symbols": ['lam', 'mssq']}

        source = "Identifier -> Sqrt[λ] / (4 * Pi) + Log[mssq]"

        self.assertEqual(reference, pythoniseExpression(source))

    def test_paseExpressionSystem(self):
        reference = [{"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                      "identifier": "Identifier",
                      "symbols": ['lam', 'mssq']},
                     {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                      "identifier": "Identifier",
                      "symbols": ['lam', 'mssq']},
                     {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                      "identifier": "Identifier",
                      "symbols": ['lam', 'mssq']}]

        source = ["Identifier -> Sqrt[λ] / (4 * Pi) + Log[mssq]",
                  "Identifier -> Sqrt[λ] / (4 * Pi) + Log[mssq]",
                  "Identifier -> Sqrt[λ] / (4 * Pi) + Log[mssq]"]

        self.assertEqual(reference, pythoniseExpressionSystem(source))

    def test_pythoniseMatrix(self):
        reference = [["1", "0"], ["0", "0"]]
        source = ["{1, 0}", "{0, 0}"]

        self.assertEqual(reference, pythoniseMatrix(source))

    def test_pythoniseMassMatrix(self):
        reference = {'definitions': [], 'matrix': "[[1, 0], [0, mssq]]"}
        source = ["{1, 0}", "{0, mssq}"]

        self.assertEqual(reference, pythoniseMassMatrix([], source))

    def test_pythoniseRotationMatrix(self):
        reference = {"matrix": {"mssq00": (0, 0), "mssq11": (1, 1)}}
        source = ["{mssq00, 0}", "{0, mssq11}"]

        self.assertEqual(reference, pythoniseRotationMatrix(source))