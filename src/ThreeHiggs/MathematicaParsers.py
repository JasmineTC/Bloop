def parseExpression(line):
    identifier = "anonymous"

    if ("->" in line):
        identifier, line = map(str.strip, line.split("->"))

    from sympy.parsing.mathematica import parse_mathematica
    from .CommonUtils import replaceGreekSymbols
    expression = parse_mathematica(replaceGreekSymbols(line))
    symbols = [str(symbol) for symbol in expression.free_symbols]

    return {"identifier": identifier, "expression": str(expression), "symbols": sorted(symbols)}

def parseExpressionSystem(lines):
    return [parseExpression(line) for line in lines]

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

def parseMassMatrix(lines):
    from sympy import Matrix
    return {"matrix": str(Matrix(parseMatrix(lines)).tolist())}

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
        reference = {"matrix": "[[1, 0], [0, mssq]]"}
        source = ["{1, 0}", "{0, mssq}"]

        from ThreeHiggs.MathematicaParsers import parseMassMatrix
        self.assertEqual(reference, parseMassMatrix(source))

    def test_parseRotationMatrix(self):
        reference = {"matrix": {"mssq00": [0, 0], "mssq11": [1, 1]}}
        source = ["{mssq00, 0}", "{0, mssq11}"]

        from ThreeHiggs.MathematicaParsers import parseRotationMatrix
        self.assertEqual(reference, parseRotationMatrix(source))

