def getSymbols(container):
    symbols = set()

    if type(container) is dict:
        if "symbols" in container:
            symbols |= set(container["symbols"])

        for key, value in container.items():
            symbols |= set(getSymbols(value))

    if type(container) is list:
        for item in container:
            symbols |= set(getSymbols(item))

    return symbols

def makeArgumentMap(parsedExpressions):
    symbols = set()

    for key, value in parsedExpressions.items():
        symbols |= getSymbols(value)

    return {symbol: index for symbol, index in enumerate(list(symbols))}

if __name__ == "__main__":
    from json import load
    from sys import argv
    print(makeArgumentMap(load(open(argv[1]))))

