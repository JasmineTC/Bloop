from ThreeHiggs.GetLines import getLines
from json import dump
def convertMathematica(args):
    veffLines = getLines(args.loFile)
    if (args.loopOrder >= 1):
        veffLines += getLines(args.nloFile)
    if (args.loopOrder >= 2):
        veffLines += getLines(args.nnloFile)
    from ThreeHiggs.MathematicaParsers import (parseExpressionSystem,
                                               parseExpressionSystemArray,
                                               parseMassMatrix,
                                               parseRotationMatrix,
                                               replaceGreekSymbols)
    
    allSymbols = [replaceGreekSymbols(symbol) for symbol in getLines(args.allSymbolsFile, mode = "json")]

    dump({"vectorMassesSquared": parseExpressionSystem(getLines(args.vectorMassesSquaredFile)),
          "vectorShortHands": parseExpressionSystem(getLines(args.vectorShortHandsFile)),
          "scalarPermutationMatrix": getLines(args.scalarPermutationFile, mode="json"),
          "scalarMassMatrixUpperLeft": parseMassMatrix(getLines(args.scalarMassMatrixUpperLeftDefinitionsFile),
                                                       getLines(args.scalarMassMatrixUpperLeftFile)),
          "scalarMassMatrixBottomRight": parseMassMatrix(getLines(args.scalarMassMatrixBottomRightDefinitionsFile),
                                                         getLines(args.scalarMassMatrixBottomRightFile)),
          "scalarRotationMatrix": parseRotationMatrix(getLines(args.scalarRotationFile)),
          "veff": parseExpressionSystem(veffLines),
          "betaFunctions4D": parseExpressionSystemArray(getLines(args.betaFunctions4DFile), allSymbols),
          "hardToSoft": parseExpressionSystem(getLines(args.hardToSoftFile)),
          "softScaleRGE": parseExpressionSystem(getLines(args.softScaleRGEFile)),
          "softToUltraSoft": parseExpressionSystem(getLines(args.softToUltraSoftFile))},
         open(args.parsedExpressionsFile, "w"),
         indent = 4)
