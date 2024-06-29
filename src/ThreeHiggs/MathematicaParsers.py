def parseExpression(line):
    identifier = "anonymous"

    if ("->" in line):
        identifier, line = map(str.strip, line.split("->"))

    from sympy.parsing.mathematica import parse_mathematica
    from .CommonUtils import replaceGreekSymbols
    expression = parse_mathematica(replaceGreekSymbols(line))
    symbols = list(expression.free_symbols)

    return {"identifier": identifier, "expression": str(expression), "symbols": symbols}

def parseMatrix(lines):
    return [[symbol for symbol in line.rstrip()
                                      .rstrip('}')
                                      .lstrip('{')
                                      .split(',')] for line in lines]

def parseConstantMatrix(lines):
    symbols = parseMatrix(lines)

    from sympy import Matrix
    sympyMatrix = Matrix(symbols)

    from numpy import array, float64
    return array(sympyMatrix.tolist()).astype(float64)

def parseMassMatrix(lines):
    from sympy import Matrix
    sympyMatrix = Matrix(parseMatrix(lines))

    return str(sympyMatrix.tolist())

def parseRotationMatrix(lines):
    from sympy import Matrix
    sympyMatrix = Matrix(parseMatrix(lines))
    shape = sympyMatrix.shape

    symbolMap = {}
    for i in range(shape[0]):
        for j in range(shape[1]):
            element = sympyMatrix[i, j]

            if element.is_symbol:
                symbolMap[str(sympyMatrix[i, j])] = i, j

    return symbolMap

