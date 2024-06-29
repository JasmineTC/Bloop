from unittest import TestCase
from unittest import main

class MathematicaParsersTestCase(TestCase):
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
        reference = {"matrix": {"mssq00": (0, 0), "mssq11": (1, 1)}}
        source = ["{mssq00, 0}", "{0, mssq11}"]

        from ThreeHiggs.MathematicaParsers import parseRotationMatrix
        self.assertEqual(reference, parseRotationMatrix(source))

if __name__ == "__main__":
    main()

