import ThreeHiggs.getLines

def convertMathematica(args):
    veffLines = ThreeHiggs.getLines.getLines(args.loFile)
    if (args.loopOrder >= 1):
        veffLines += ThreeHiggs.getLines.getLines(args.nloFile)
    if (args.loopOrder >= 2):
        veffLines += ThreeHiggs.getLines.getLines(args.nnloFile)

    from ThreeHiggs.MathematicaParsers import (parseExpressionSystem,
                                               parseConstantMatrix,
                                               parseMassMatrix,
                                               parseRotationMatrix)


    from json import dump
    dump({"vectorMassesSquared": parseExpressionSystem(ThreeHiggs.getLines.getLines(args.vectorMassesSquaredFile)),
          "vectorShortHands": parseExpressionSystem(ThreeHiggs.getLines.getLines(args.vectorShortHandsFile)),
          "scalarPermutationMatrix": parseConstantMatrix(ThreeHiggs.getLines.getLines(args.scalarPermutationFile)),
          "scalarMassMatrixUpperLeft": parseMassMatrix(ThreeHiggs.getLines.getLines(args.scalarMassMatrixUpperLeftDefinitionsFile),
                                                       ThreeHiggs.getLines.getLines(args.scalarMassMatrixUpperLeftFile)),
          "scalarMassMatrixBottomRight": parseMassMatrix(ThreeHiggs.getLines.getLines(args.scalarMassMatrixBottomRightDefinitionsFile),
                                                         ThreeHiggs.getLines.getLines(args.scalarMassMatrixBottomRightFile)),
          "scalarRotationMatrix": parseRotationMatrix(ThreeHiggs.getLines.getLines(args.scalarRotationFile)),
          "veff": parseExpressionSystem(veffLines),
          "hardToSoft": parseExpressionSystem(ThreeHiggs.getLines.getLines(args.hardToSoftFile)),
          "softScaleRGE": parseExpressionSystem(ThreeHiggs.getLines.getLines(args.softScaleRGEFile)),
          "softToUltraSoft": parseExpressionSystem(ThreeHiggs.getLines.getLines(args.softToUltraSoftFile))},
         open(args.parsedExpressionsFile, "w"),
         indent = 4)

