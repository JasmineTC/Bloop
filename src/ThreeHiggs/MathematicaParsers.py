def parseConstantMatrix(lines):
    symbols = [[symbol for symbol in line.rstrip()
                                         .rstrip('}')
                                         .lstrip('{')
                                         .split(',')] for line in lines]

    from sympy import Matrix
    sympyMatrix = Matrix(symbols)

    from numpy import array, float64
    return array(sympyMatrix.tolist()).astype(float64)

def parseExpression(line):
    identifier = "anonymous"

    if ("->" in line):
        identifier, line = map(str.strip, line.split("->"))

    from sympy.parsing.mathematica import parse_mathematica
    from .CommonUtils import replaceGreekSymbols
    expression = parse_mathematica(replaceGreekSymbols(line))
    symbols = list(expression.free_symbols)

    return {"identifier": identifier, "expression": str(expression), "symbols": symbols}

