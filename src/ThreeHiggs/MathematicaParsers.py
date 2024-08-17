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
    string = string.replace("^2", "sq")

    #""" For ultrasoft theory DRalgo appends "US" => remove that too. Gauge couplings again need special treatment."""
    string = string[:-len("US")] if string.endswith("US") else string
    string = string.replace("USsq", "sq") if string.endswith("USsq") else string

    ## Remove "3d" suffix with even crazier oneliner (suffix meaning that it's removed only from end of the string)
    string = string[:-len("3d")] if string.endswith("3d") else string
    ## Gauge couplings are originally of form g3d^2 so account for that too 
    string = string.replace("3dsq", "sq") if string.endswith("3dsq") else string
    
    return string

def parseExpression(line, remove3DSuffices = False):
    identifier = "anonymous"

    if ("->" in line):
        identifier, line = map(str.strip, line.split("->"))

    from sympy.parsing.mathematica import parse_mathematica
    identifier = removeSuffices(replaceGreekSymbols(identifier))
    expression = parse_mathematica(replaceGreekSymbols(line).replace("3d", ""))
    symbols = [str(symbol) for symbol in expression.free_symbols]

    return {"identifier": identifier, "expression": str(expression), "symbols": sorted(symbols)}

def parseExpressionSystem(lines, remove3DSuffices = False):
    return [parseExpression(line, remove3DSuffices) for line in lines]

def parseMatrix(lines):
    return [[symbol.strip() for symbol in line.strip()
                                              .strip('}')
                                              .strip('{')
                                              .split(',')] for line in lines]

def parseConstantMatrix(lines):
    matrix = parseMatrix(lines)

    from sympy import Matrix
    sympyMatrix = Matrix(matrix)

    from numpy import array, float64
    return {"matrix": array(sympyMatrix.tolist()).astype(float64).tolist()}

def parseMassMatrix(definitionsLines, matrixLines):
    matrix = parseMatrix(matrixLines)

    from sympy import Matrix
    sympyMatrix = Matrix(matrix)

    from numpy import array, float64
    return {"definitions": parseExpressionSystem(definitionsLines),
            "matrix": str(array(sympyMatrix.tolist()).tolist())}

def parseRotationMatrix(lines):
    from sympy import Matrix
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
    def test_removeSuffices(self):
        reference = ["myVarsq", "sqmyVar",
                     "myVarsq", "myVar", "3dsqmyVar", "3dmyVar",
                     "myVarsq", "myVar", "USsqmyVar", "USmyVar"]

        source = ["myVar^2", "^2myVar", 
                  "myVar3dsq", "myVar3d", "3dsqmyVar", "3dmyVar",
                  "myVarUSsq", "myVarUS", "USsqmyVar", "USmyVar"]

        self.assertEqual(reference, [removeSuffices(sourceString) for sourceString in source])

    def test_parseExpression(self):
        reference = {"expression": "sqrt(lam)/(4*pi) + log(mssq)",
                     "identifier": "Identifier",
                     "symbols": ['lam', 'mssq']}

        source = "Identifier -> Sqrt[λ] / (4 * Pi) + Log[mssq]"

        from ThreeHiggs.MathematicaParsers import parseExpression
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

        from ThreeHiggs.MathematicaParsers import parseExpressionSystem
        self.assertEqual(reference, parseExpressionSystem(source))

    def test_parseMatrix(self):
        reference = [["1", "0"], ["0", "0"]]
        source = ["{1, 0}", "{0, 0}"]

        from ThreeHiggs.MathematicaParsers import parseMatrix
        self.assertEqual(reference, parseMatrix(source))

    def test_parseConstantMatrix(self):
        reference = {"matrix": [[1.0, 0.0], [0.0, 0.0]]}
        source = ["{1, 0}", "{0, 0}"]

        from ThreeHiggs.MathematicaParsers import parseConstantMatrix
        self.assertEqual(reference, parseConstantMatrix(source))

    def test_parseMassMatrix(self):
        reference = {'definitions': [], 'matrix': "[[1, 0], [0, mssq]]"}
        source = ["{1, 0}", "{0, mssq}"]

        from ThreeHiggs.MathematicaParsers import parseMassMatrix
        self.assertEqual(reference, parseMassMatrix([], source))

    def test_parseRotationMatrix(self):
        reference = {"matrix": {"mssq00": [0, 0], "mssq11": [1, 1]}}
        source = ["{mssq00, 0}", "{0, mssq11}"]

        from ThreeHiggs.MathematicaParsers import parseRotationMatrix
        self.assertEqual(reference, parseRotationMatrix(source))

