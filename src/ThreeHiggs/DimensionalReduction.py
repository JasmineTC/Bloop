from dataclasses import dataclass, InitVar
@dataclass(frozen=True)
class DimensionalReduction():
    hardToSoft: callable = 0
    softScaleRGE: callable = 0 
    softToUltraSoft: callable = 0 
    config: InitVar[dict] = None
    
    def __post_init__(self, config: dict):
        if config:
            self.__init__(**config)

    def getUltraSoftParams(self, paramsForMatching: dict[str, float], goalRGScale: float) -> dict[str, float]:
        outParams = self.hardToSoft.evaluate(paramsForMatching, bReturnDict = True)
        ## TODO Talk to someone about this RGScale stuff!!!!!
        outParams |= {"RGScale": paramsForMatching["RGScale"],
                      "goalScale": goalRGScale,
                      "startScale": paramsForMatching["RGScale"]}
        
        outParams |= self.softScaleRGE.evaluate(outParams, bReturnDict = True)
        ## For reasons unknown the RG scale is different for soft and ultra soft physics
        outParams["RGScale"] = goalRGScale
        ## HACK this is the RG scale name in Veff
        return self.softToUltraSoft.evaluate(outParams, bReturnDict = True)




from ThreeHiggs.GetLines import getLines
    
from unittest import TestCase
class DimensionalReductionUnitTest(TestCase):
    def test_getEFTParams(self):
        from ThreeHiggs.PythoniseMathematica import pythoniseExpressionSystem
        from ThreeHiggs.ParsedExpression import ParsedExpressionSystem
        dimensionalReduction = DimensionalReduction(ParsedExpressionSystem(pythoniseExpressionSystem(["b -> a", "T->T"])),
                                                    ParsedExpressionSystem(pythoniseExpressionSystem(["c -> b"])), 
                                                    ParsedExpressionSystem(pythoniseExpressionSystem(["a -> c", "mu3US -> T"])))
        reference = {'a': 1, 'mu3US': 0}


        self.assertEqual(reference,
                         dimensionalReduction.getUltraSoftParams({"a": 1, "RGScale": -1, "T":0}, 0))
        
    def test_getEFTParamsFull(self):
        from ThreeHiggs.PythoniseMathematica import pythoniseExpressionSystem
        from ThreeHiggs.ParsedExpression import ParsedExpressionSystem
        
        dimensionalReduction = DimensionalReduction(ParsedExpressionSystem(pythoniseExpressionSystem(getLines("Data/HardToSoft/softScaleParams_NLO.txt")), None),
                                                    ParsedExpressionSystem(pythoniseExpressionSystem(getLines("Data/HardToSoft/softScaleRGE.txt")), None), 
                                                    ParsedExpressionSystem(pythoniseExpressionSystem(getLines("Data/SoftToUltrasoft/ultrasoftScaleParams_NLO.txt")), None))
        source = {'yt3': 0.9270871707819464, 'g1': 0.3524790946562762, 'g2': 0.6462351552111576, 'g3': 1.1232276401951706, 'lam1Re': 0.09958019847087743, 'lam1Im': 0.0, 'lam2Re': 0.0, 'lam2Im': 0.0, 'lam11': 0.11179799168185336, 'lam22': 0.12191841915170326, 'lam12': 0.13403604439831854, 'lam12p': 0.14210988273784536, 'lam23': 0.0029276560925431193, 'lam23p': 0.0013512143449770024, 'lam3Re': 0.0, 'lam3Im': 0.0, 'lam31': 0.0029258903606322086, 'lam31p': 0.0013509802864778528, 'lam33': 0.09672085995302131, 'mu12sqRe': 0.0, 'mu12sqIm': 0.0, 'mu2sq': -90136.57045751583, 'mu3sq': 8127.606718757464, 'mu1sq': -90044.04861303422, 'RGScale': 352.7753977724091, 'T': 50.0, 'Lb': -4.440892098500626e-16, 'Lf': 2.7725887222397807}
        reference = {'mu3US': 50, 'lam11': (5.523844829176312+0j), 'lam12': (6.621577064652674+0j), 'lam12p': (7.053488262004727+0j), 'lam1Im': 0.0, 'lam1Re': 4.979009923543871, 'lam22': (6.029744638552994+0j), 'lam23': (0.06575174120618288+0j), 'lam23p': (0.011891632008381406+0j), 'lam2Im': 0.0, 'lam2Re': 0.0, 'lam31': (0.06576664548054172+0j), 'lam31p': (0.011898384154257216+0j), 'lam33': (6.281251508922814+0j), 'lam3Im': 0.0, 'lam3Re': 0.0, 'g1': (2.4742142303484282+0j), 'g2': (4.493373011738393+0j), 'g3': (7.560030477241674+0j), 'mu12sqIm': 0.0, 'mu12sqRe': 0.0, 'mu1sq': (-90440.68446664244+0j), 'mu2sq': (-90545.99168869366+0j), 'mu3sq': (6816.682950673874+0j)}

        self.assertEqual(reference, dimensionalReduction.getUltraSoftParams(source, 50))

