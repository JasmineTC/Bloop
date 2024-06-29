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
                symbolMap[str(sympyMatrix[i, j])] = i, j

    return {"matrix": symbolMap}

