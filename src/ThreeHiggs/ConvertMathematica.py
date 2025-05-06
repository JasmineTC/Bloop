from ThreeHiggs.GetLines import getLines
from json import dump
def convertMathematica(args):
    veffLines = getLines(args.loFile)
    veffLines += getLines(args.nloFile)
    if (args.loopOrder >= 2):
        veffLines += getLines(args.nnloFile)
    from ThreeHiggs.MathematicaParsers import (parseExpressionSystem,
                                               parseExpressionSystemArray,
                                               parseMassMatrix,
                                               parseRotationMatrix,
                                               replaceGreekSymbols)
    
    allSymbols = [replaceGreekSymbols(symbol) for symbol in getLines(args.allSymbolsFile, mode = "json")]

    dump({"betaFunctions4D": {"expressions": parseExpressionSystemArray(getLines(args.betaFunctions4DFile), allSymbols),
                              "fileName": args.betaFunctions4DFile},
          "hardToSoft": {"expressions":  parseExpressionSystem(getLines(args.hardToSoftFile)),
                         "fileName": args.hardToSoftFile},
          "softScaleRGE": {"expressions": parseExpressionSystem(getLines(args.softScaleRGEFile)),
                           "fileName": args.softScaleRGEFile},
          "softToUltraSoft": {"expressions": parseExpressionSystem(getLines(args.softToUltraSoftFile)),
                              "fileName": args.softToUltraSoftFile},
          "vectorMassesSquared": {"expressions": parseExpressionSystem(getLines(args.vectorMassesSquaredFile)),
                                  "fileName": args.vectorMassesSquaredFile},
          "vectorShortHands": {"expressions": parseExpressionSystem(getLines(args.vectorShortHandsFile)),
                               "fileName": args.vectorShortHandsFile},
          "veff": {"expressions": parseExpressionSystem(veffLines),
                     "fileName": "Combined Veff files"},
          # "": {"expressions": ,
          #            "fileName": },
          # "": {"expressions": ,
          #            "fileName": },
        
        
          "scalarPermutationMatrix": getLines(args.scalarPermutationFile, mode="json"),
          "scalarMassMatrixUpperLeft": parseMassMatrix(getLines(args.scalarMassMatrixUpperLeftDefinitionsFile),
                                                       getLines(args.scalarMassMatrixUpperLeftFile)),
          "scalarMassMatrixBottomRight": parseMassMatrix(getLines(args.scalarMassMatrixBottomRightDefinitionsFile),
                                                         getLines(args.scalarMassMatrixBottomRightFile)),
          "scalarRotationMatrix": parseRotationMatrix(getLines(args.scalarRotationFile))},
         open(args.parsedExpressionsFile, "w"),
         indent = 4)
