from ThreeHiggs.GetLines import getLines
from json import dump
from sympy import Matrix
from sympy.parsing.mathematica import parse_mathematica
from numpy import euler_gamma, pi
from pathlib import Path

def replaceGreekSymbols(string: str) -> str:
    ## TODO use unicodedata package here to do magic. 
    # Or bully DRalgo people to removing greek symbols
    return string.replace(u"\u03BB", "lam").replace(u"\u03BC", "mu")

def replaceSymbolsConst(string):
    return string.replace("Pi", str(pi)) \
                 .replace("EulerGamma", str(euler_gamma)) \
                 .replace("Glaisher", "1.28242712910062") 

def removeSuffices(string):
    return string.replace("^2", "sq")

def replaceSymbolsWithIndices(expression, symbols):
    expression = replaceGreekSymbols(expression)
    ## Reverse needed to deal with lam23 and lam23p i.e. substring replaces larger full string
    for idx, symbol in enumerate(sorted(symbols, reverse = True)):
        expression = expression.replace(symbol , f"params[{idx}]")

    return expression

def pythoniseExpressionArray(line, allSymbols):
    identifier, line = map(str.strip, line.split("->")) if ("->" in line) else ("missing", line)
    
    identifier = removeSuffices(replaceGreekSymbols(identifier))
    expression = parse_mathematica(replaceSymbolsConst(replaceGreekSymbols(line)))
    symbols = [str(symbol) for symbol in expression.free_symbols]

    return {"identifier": identifier, 
            "expression": replaceSymbolsWithIndices(str(expression), allSymbols), 
            "symbols": sorted(symbols)}

def pythoniseExpression(line):
    identifier, line = map(str.strip, line.split("->")) if ("->" in line) else ("missing", line)

    identifier = removeSuffices(replaceGreekSymbols(identifier))
    expression = parse_mathematica(replaceSymbolsConst(replaceGreekSymbols(line)))
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
    
    allSymbols = getLines(args.allSymbolsFile, mode = "json") + ["missing"]
    allSymbols = sorted([replaceGreekSymbols(symbol) for symbol in allSymbols], 
                        reverse = True)

    (outputFile := Path(args.pythonisedExpressionsFile)).parent.mkdir(exist_ok=True, 
                                                             parents=True)  
    ## Move get lines to the functions? -- Would need to rework veffLines in this case
    ## Not ideal to have nested dicts but is future proof for when we move to arrays
    dump(
        {
            "bounded": {
                "expressions": pythoniseExpressionSystemArray(getLines(args.boundedConditions), allSymbols),
                "fileName": "bounded",
            },
            "betaFunctions4D": {
                "expressions": pythoniseExpressionSystemArray(getLines(args.betaFunctions4DFile), allSymbols),
                "fileName": args.betaFunctions4DFile,
            },
            "hardToSoft": {
                "expressions":  pythoniseExpressionSystemArray(getLines(args.hardToSoftFile), allSymbols),
                "fileName": args.hardToSoftFile,
            },
            "softScaleRGE": {
                "expressions": pythoniseExpressionSystemArray(getLines(args.softScaleRGEFile), allSymbols),
                "fileName": args.softScaleRGEFile,
            },
            "softToUltraSoft": {
                "expressions": pythoniseExpressionSystemArray(getLines(args.softToUltraSoftFile), allSymbols),
                "fileName": args.softToUltraSoftFile,
            },
            "vectorMassesSquared": {
                "expressions": pythoniseExpressionSystemArray(getLines(args.vectorMassesSquaredFile), allSymbols),
                "fileName": args.vectorMassesSquaredFile,
            },
            "vectorShortHands": {
                "expressions": pythoniseExpressionSystemArray(getLines(args.vectorShortHandsFile), allSymbols),
                "fileName": args.vectorShortHandsFile,
            },
            "veff": {
                "expressions": pythoniseExpressionSystem(veffLines),
                "fileName": "Combined Veff files",
            },
            "veffArray": {
                "expressions": pythoniseExpressionSystemArray(veffLines, allSymbols),
                "fileName": "Combined Veff files",
            },
            "scalarMassMatrixUpperLeft": {
                "expressions": pythoniseMassMatrix(
                    getLines(args.scalarMassMatrixUpperLeftDefinitionsFile),
                    getLines(args.scalarMassMatrixUpperLeftFile),
                ),
                "fileName": (
                    args.scalarMassMatrixUpperLeftDefinitionsFile,
                    args.scalarMassMatrixBottomRightFile),
            },
            "scalarMassMatrixBottomRight": {
                "expressions": pythoniseMassMatrix(
                    getLines(args.scalarMassMatrixBottomRightDefinitionsFile),
                    getLines(args.scalarMassMatrixBottomRightFile),
                ),
                "fileName": (
                    args.scalarMassMatrixBottomRightDefinitionsFile,
                    args.scalarMassMatrixBottomRightFile,
                ),
            },
            "scalarRotationMatrix": {
                "expressions": pythoniseRotationMatrix(getLines(args.scalarRotationFile)),
                "fileName": args.scalarRotationFile,
            },
            "scalarPermutationMatrix": getLines(args.scalarPermutationFile, mode="json"),
            "allSymbols": {
                "allSymbols": allSymbols,
                "fileName": args.allSymbolsFile,
            },
        },
        open(outputFile, "w"),
        indent = 4
    )

    from Veff_generation import generate_veff_module, compile_veff_submodule
    generate_veff_module(args, allSymbols)
    compile_veff_submodule()

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
        reference = {"expression": "0.07957747154594767*sqrt(lam) + log(mssq)",
                     "identifier": "Identifier",
                     "symbols": ['lam', 'mssq']}

        source = "Identifier -> Sqrt[λ] / (4 * Pi) + Log[mssq]"

        self.assertEqual(reference, pythoniseExpression(source))

    def test_paseExpressionSystem(self):
        reference = [{"expression": "0.07957747154594767*sqrt(lam) + log(mssq)",
                      "identifier": "Identifier",
                      "symbols": ['lam', 'mssq']},
                     {"expression": "0.07957747154594767*sqrt(lam) + log(mssq)",
                      "identifier": "Identifier",
                      "symbols": ['lam', 'mssq']},
                     {"expression": "0.07957747154594767*sqrt(lam) + log(mssq)",
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
