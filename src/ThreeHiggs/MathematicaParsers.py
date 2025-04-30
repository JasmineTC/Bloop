from sympy import Matrix
def replaceGreekSymbols(string: str) -> str:
    #import unicodedata
    ## Unicode magic, this is definitely not ideal
    lowerCaseMu = u"\u03BC"
    lowerCaseLambda = u"\u03BB"

    newString = string 

    newString = newString.replace(lowerCaseLambda, "lam")
    newString = newString.replace(lowerCaseMu, "mu")
    
    """ TODO use unicodedata package here to do magic.
    """
    # NOTE: Manual replacements are definitely not a general solution. Consider problematic case: expression that contains both unicode lambda and separate symbol "lam" 
    # So tbh I'd like to keep the symbols are they are. But parse_mathematica from sympy does not seem to manage greek symbols at all!! 
    return newString

def removeSuffices(string):
    return string.replace("^2", "sq")

def replaceSymbolsWithIndices(expression, symbols):
    expression = replaceGreekSymbols(expression)
    ## Reverse needed to deal with lam23 and lam23p i.e. substring replaces larger full string
    for idx, symbol in enumerate(sorted(symbols, reverse = True, key = lambda symbol: replaceGreekSymbols(symbol))):
        expression = expression.replace(replaceGreekSymbols(symbol) , f"params[{idx}]")

    return expression

def parseExpressionArray(line, allSymbols, remove3DSuffices = False):
    identifier = "anonymous"
    if ("->" in line):
        identifier, line = map(str.strip, line.split("->"))
    from sympy.parsing.mathematica import parse_mathematica
    identifier = removeSuffices(replaceGreekSymbols(identifier))
    expression = parse_mathematica(replaceGreekSymbols(line).replace("3d", ""))
    symbols = [str(symbol) for symbol in expression.free_symbols]

    return {"identifier": identifier, 
            "expression": replaceSymbolsWithIndices(str(expression), allSymbols), 
            "symbols": sorted(symbols)}

def parseExpression(line, remove3DSuffices = False):
    identifier = "anonymous"

    if ("->" in line):
        identifier, line = map(str.strip, line.split("->"))

    from sympy.parsing.mathematica import parse_mathematica
    identifier = removeSuffices(replaceGreekSymbols(identifier))
    expression = parse_mathematica(replaceGreekSymbols(line))
    symbols = [str(symbol) for symbol in expression.free_symbols]

    return {"identifier": identifier, "expression": str(expression), "symbols": sorted(symbols)}

def parseExpressionSystemArray(lines, allSymbols, remove3DSuffices = False):
    return [parseExpressionArray(line, allSymbols, remove3DSuffices) for line in lines]

def parseExpressionSystem(lines, remove3DSuffices = False):
    return [parseExpression(line, remove3DSuffices) for line in lines]

def parseMatrix(lines):
    return [[symbol.strip() for symbol in line.strip()
                                              .strip('}')
                                              .strip('{')
                                              .split(',')] for line in lines]

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
        reference = {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                     "identifier": "Identifier",
                     "symbols": ['lam', 'mssq']}

        source = "Identifier -> Sqrt[λ] / (4 * Pi) + Log[mssq]"

        self.assertEqual(reference, parseExpression(source))

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

        self.assertEqual(reference, parseExpressionSystem(source))

    def test_parseMatrix(self):
        reference = [["1", "0"], ["0", "0"]]
        source = ["{1, 0}", "{0, 0}"]

        self.assertEqual(reference, parseMatrix(source))

    def test_parseMassMatrix(self):
        reference = {'definitions': [], 'matrix': "[[1, 0], [0, mssq]]"}
        source = ["{1, 0}", "{0, mssq}"]

        self.assertEqual(reference, parseMassMatrix([], source))

    def test_parseRotationMatrix(self):
        reference = {"matrix": {"mssq00": [0, 0], "mssq11": [1, 1]}}
        source = ["{mssq00, 0}", "{0, mssq11}"]

        self.assertEqual(reference, parseRotationMatrix(source))

