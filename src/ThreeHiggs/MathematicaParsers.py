from sympy import Matrix
from sympy.parsing.mathematica import parse_mathematica
from numpy import euler_gamma, pi
Glaisher = "1.28242712910062"

def replaceSymbolsConst(string):
    ## Change expressions to use either pi or Pi but not both!!!
    return string.replace("pi", str(pi)).replace("Pi", str(pi)) \
                 .replace("EulerGamma", str(euler_gamma)) \
                 .replace("Glaisher", Glaisher) 

def replaceGreekSymbols(string):
    """ TODO use unicodedata package here to do magic."""
    return string.replace(u"\u03BB", "lam").replace(u"\u03BC", "mu")

def removeSuffices(string):
    return string.replace("^2", "sq")

def replaceSymbolsWithIndices(expression, symbols):
    ## Reverse needed to deal with lam23 and lam23p i.e. substring replaces larger full string
    for idx, symbol in enumerate(sorted(symbols, reverse = True, key = lambda symbol: replaceGreekSymbols(symbol))):
        expression = expression.replace(replaceGreekSymbols(symbol) , f"params[{idx}]")

    return expression

def parseExpressionArray(line, allSymbols):
    ## Moving removeSuffices here breaks code (you get params[1]sq)
    line = replaceSymbolsConst(replaceGreekSymbols(line))
    identifier, expression = map(str.strip, line.split("->")) if "->" in line else ("anonymous", line)

    expression = parse_mathematica(expression)

    return {"identifier": removeSuffices(identifier), 
            "expression": replaceSymbolsWithIndices(str(expression), allSymbols), 
            "symbols": sorted([str(symbol) for symbol in expression.free_symbols])}

def parseExpression(line, remove3DSuffices = False):
    ## Moving removeSuffices here breaks code (you get params[1]sq)
    line = replaceSymbolsConst(replaceGreekSymbols(line))
    identifier, expression = map(str.strip, line.split("->")) if "->" in line else ("anonymous", line)

    expression = parse_mathematica(expression)

    return {"identifier": removeSuffices(identifier), 
            "expression": str(expression),
            "symbols": sorted([str(symbol) for symbol in expression.free_symbols])}

def parseExpressionSystemArray(lines, allSymbols, remove3DSuffices = False):
    return [parseExpressionArray(line, allSymbols) for line in lines]

def parseExpressionSystem(lines, remove3DSuffices = False):
    return [parseExpression(line, remove3DSuffices) for line in lines]

def parseMatrix(lines):
    return [[symbol.strip() for symbol in line.strip()
                                              .strip('}')
                                              .strip('{')
                                              .split(',')] for line in lines]

def parseConstantMatrix(lines):
    from numpy import array
    ## Convert to array to specify type, convert to list to be compatable with json
    return {"matrix": array(Matrix(parseMatrix(lines)), dtype = int).tolist()}

def parseMassMatrix(definitionsLines, matrixLines):
    from sympy.core.sympify import sympify
    ##Need sympify that the parsed matrix line is ufuncable
    ##Need str to make compatable with json 
    return {"definitions": parseExpressionSystem(definitionsLines),
            "matrix": str(sympify(parseMatrix(matrixLines)))}

def parseRotationMatrix(lines):
    sympyMatrix = Matrix(parseMatrix(lines))
    shape = sympyMatrix.shape

    symbolMap = {}
    for i in range(shape[0]):
        for j in range(shape[1]):
            element = sympyMatrix[i, j]

            if element.is_symbol:
                symbolMap[str(sympyMatrix[i, j])] = [i, j]

    return {"matrix": symbolMap}

from unittest import TestCase
class MathematicaParsersUnitTests(TestCase):
    def test_replaceGreekSymbols(self):
        reference = ["lam", "lam lam", "mu", "mu mu", "lam mu", "mu lam"]
        source = ["λ", "λ λ", "μ", "μ μ", "λ μ", "μ λ"]

        self.assertEqual(reference, [replaceGreekSymbols(sourceString) for sourceString in source])

    def test_removeSuffices(self):
        reference = ["myVarsq", "sqmyVar"]

        source = ["myVar^2", "^2myVar"]

        self.assertEqual(reference, [removeSuffices(sourceString) for sourceString in source])

    def test_parseExpression(self):
        reference = {"expression": "sqrt(lam) + log(mssq)",
                     "identifier": "Identifier",
                     "symbols": ['lam', 'mssq']}

        source = "Identifier -> Sqrt[λ] + Log[mssq]"

        self.assertEqual(reference, parseExpression(source))

    def test_paseExpressionSystem(self):
        reference = [{"expression": "sqrt(lam) + log(mssq)",
                      "identifier": "Identifier",
                      "symbols": ['lam', 'mssq']},
                     {"expression": "sqrt(2)*sqrt(lam) + log(mssq)",
                      "identifier": "Identifier",
                      "symbols": ['lam', 'mssq']}]

        source = ["Identifier -> Sqrt[λ] + Log[mssq]",
                  "Identifier -> Sqrt[2]*Sqrt[λ]+ Log[mssq]"]

        self.assertEqual(reference, parseExpressionSystem(source))

    def test_parseMatrix(self):
        reference = [["1", "0"], ["0", "0"]]
        source = ["{1, 0}", "{0, 0}"]

        self.assertEqual(reference, parseMatrix(source))

    def test_parseConstantMatrix(self):
        reference = {"matrix": [[1.0, 0.0], [0.0, 0.0]]}
        source = ["{1, 0}", "{0, 0}"]

        self.assertEqual(reference, parseConstantMatrix(source))

    def test_parseMassMatrix(self):
        reference = {'definitions': [], 'matrix': "[[1, 0], [0, mssq]]"}
        source = ["{1, 0}", "{0, mssq}"]

        self.assertEqual(reference, parseMassMatrix([], source))

    def test_parseRotationMatrix(self):
        reference = {"matrix": {"mssq00": [0, 0], "mssq11": [1, 1]}}
        source = ["{mssq00, 0}", "{0, mssq11}"]

        self.assertEqual(reference, parseRotationMatrix(source))

