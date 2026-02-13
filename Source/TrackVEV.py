import numpy as np
import scipy
import nlopt
from dataclasses import dataclass, InitVar

from PythoniseMathematica import replaceGreekSymbols
from ParsedExpression import ParsedExpression, ParsedExpressionSystem

class cNlopt:
    def __init__(self,
        nloptConfig
    ):
        optG = nlopt.opt(nlopt.GN_DIRECT_NOSCAL, nloptConfig["nbrVars"])
        optG.set_lower_bounds(nloptConfig["varLowerBounds"])
        optG.set_upper_bounds(nloptConfig["varUpperBounds"])
        optG.set_xtol_abs(nloptConfig["absGlobalTol"])
        optG.set_xtol_rel(nloptConfig["relGlobalTol"])
        self.optG = optG

        optL = nlopt.opt(nlopt.LN_BOBYQA, nloptConfig["nbrVars"])
        optL.set_lower_bounds(nloptConfig["varLowerBounds"])
        optL.set_upper_bounds(nloptConfig["varUpperBounds"])
        optL.set_xtol_abs(nloptConfig["absLocalTol"])
        optL.set_xtol_rel(nloptConfig["relLocalTol"])
        self.optL = optL
    
    def funcInit(self, func):
        self.optG.set_min_objective(func)
        self.optL.set_min_objective(func)
    
    def nloptGlobal(self, initialGuess: list[float]):
        return self.nloptLocal(self.optG.optimize(initialGuess))

    def nloptLocal(self, initialGuess: list[float]):
        return self.optL.optimize(initialGuess), self.optL.last_optimum_value()

def bIsPerturbative(params, pertSymbols, allSymbols):
    for pertSymbol in pertSymbols:
        if abs(params[allSymbols.index(pertSymbol)]) > 4 * np.pi:
            return False

    return True

class TrackVEV:
    def __init__(self, TRange,
             initialGuesses,
             verbose,
             pythonisedExpressions,
             nloptConfig,
        ):
        
        self.verbose = verbose
        self.TRange = TRange
        self.initialGuesses = initialGuesses
        self.nloptInst = cNlopt(nloptConfig)
        self.verbose = verbose
        
        self.allSymbols = pythonisedExpressions["allSymbols"]["allSymbols"]
        
        self.pertSymbols = {replaceGreekSymbols(symbol) 
                            for symbolSet in ("fourPointSymbols", "yukawaSymbols", "gaugeSymbols")
                            for symbol in pythonisedExpressions["lagranianVariables"]["lagranianVariables"][symbolSet]
                            }
        
        self.hardToSoft = ParsedExpressionSystem(
                             pythonisedExpressions["hardToSoft"],
                             self.allSymbols,
                         )
        
        self.hardScale = ParsedExpression(
                             pythonisedExpressions["hardScale"]["expressions"]["expression"],
                             pythonisedExpressions["hardScale"]["filePath"],
                         )
        
        self.softScaleRGE = ParsedExpressionSystem(
                             pythonisedExpressions["softScaleRGE"],
                             self.allSymbols,
                         )
        if pythonisedExpressions["softToUltraSoft"] == "none":
            self.softToUltraSoft = None
        else:
            self.softToUltraSoft = ParsedExpressionSystem(
                             pythonisedExpressions["softToUltraSoft"],
                             self.allSymbols,
                         )
        
        self.betaFunction4DExpression = ParsedExpressionSystem(
                             pythonisedExpressions["betaFunctions4D"],
                             self.allSymbols,
                         )
        
        self.bounded = ParsedExpressionSystem(
                             pythonisedExpressions["bounded"],
                             self.allSymbols,
                         )
        
    def trackVEV(self, benchmark):
        minimizationResults = {
            "T": [],
            "vevDepthReal": [],
            "vevDepthImag": [],
            "vevLocation": [],
            "bIsPerturbative": [],
            "failureReason": False,
        }

        params = np.zeros(len(self.allSymbols), dtype="float64")
        for key, value in benchmark["lagranianParameters"].items():
            if key == "RGScale":
                continue
            params[self.allSymbols.index(key)] = value

        ## What to do if user RGScale > 7.3TMax? Idk why someone might do this though

        muRange = np.linspace(
            benchmark["lagranianParameters"]["RGScale"],
            7.3 * self.TRange[-1],
            len(self.TRange) * 10,
        )

        ## Dev note on solve_ivp numerical instability:
        ## Any none-zero value in initalConditions will be updated even if not
        ## included in beta function 
        ## The number of symbols in all symbols also affects the result
        ## These differences are at or below the tol of sol_ivp
        def betaFunction(
                mu, 
                initialConditions
            ):
                return self.betaFunction4DExpression.evaluate(initialConditions) / mu
                
        solvedBetaFunction = scipy.integrate.solve_ivp(
            betaFunction,
            (muRange[0], muRange[-1]),
            params,
            t_eval=muRange
        )
        
        if not solvedBetaFunction.success:
            return minimizationResults | {"failureReason":  solvedBetaFunction.message}
        
        betaSpline4D = {
            symbol: scipy.interpolate.CubicSpline(muRange, solvedBetaFunction.y[idx])
            for idx, symbol in enumerate(self.allSymbols)
            if np.any(solvedBetaFunction.y[idx] != solvedBetaFunction.y[idx][0])
        }

        counter = 0
        ## Initialise vevLocation to feed into the minimisation algo so it can
        ## use the location of the previous minimum as a guess for the next
        ## Not ideal as the code has to repeat an initial guess on first T
        vevLocation = np.array(self.initialGuesses[0])

        for T in self.TRange:
            if self.verbose:
                print(f"Start of temp = {T} loop")

            params = np.zeros(len(self.allSymbols), dtype="float64")
            params[self.allSymbols.index("T")] = T
            
            for key, spline in betaSpline4D.items():
                params[self.allSymbols.index(key)] = spline(self.hardScale.evaluate(params))
            
            if not np.all(self.bounded.evaluateUnordered(params)):
                return minimizationResults | {"failureReason": "unBounded"}
            
            minimizationResults["bIsPerturbative"].append(bIsPerturbative(params, 
                                                                          self.pertSymbols, 
                                                                          self.allSymbols)
                                                          )

            params = self.hardToSoft.evaluate(params)
            params = self.softScaleRGE.evaluate(params)
            if self.softToUltraSoft:
                params = self.softToUltraSoft.evaluate(params)
            
            ## Round needed because nlopt result sometimes fp out of bounds
            ## See https://github.com/stevengj/nlopt/issues/625
            vevLocation, vevDepth = self.findGlobalMinimum(
                params, self.initialGuesses + [np.round(vevLocation, 8)]
            )
           
            minimizationResults["T"].append(T)
            minimizationResults["vevDepthReal"].append(vevDepth.real)
            minimizationResults["vevDepthImag"].append(vevDepth.imag)
            minimizationResults["vevLocation"].append(vevLocation)

            if np.all(np.abs(vevLocation) < 0.5):
                if self.verbose:
                    print(f"Symmetric phase found at temp {T}")

                if counter == 3:
                    break

                counter += 1

        minimizationResults["vevLocation"] = np.transpose(
            minimizationResults["vevLocation"]
        ).tolist()

        return minimizationResults
    
    def findGlobalMinimum(self, params, minimumCandidates):
        from evaluatePotential  import evaluatePotential
        """For physics reasons we only minimise the real part,
        for nlopt reasons we need to give a redunant grad arg"""
        def veffWrapper(fields, grad):
            return np.real(
                    evaluatePotential(fields, params)
                )

        self.nloptInst.funcInit(veffWrapper)     
        bestResult = self.nloptInst.nloptGlobal(minimumCandidates[0])

        for candidate in minimumCandidates:
            result = self.nloptInst.nloptLocal(candidate)
            if result[1] < bestResult[1]:
                bestResult = result

        ## Potential computed again in case its complex
        return bestResult[0], evaluatePotential(bestResult[0], params)
    


from unittest import TestCase
class TrackVEVUnitTests(TestCase):
    def test_bIsPerturbativeTrue(self):
        reference = True
        source = [0.7, -0.8, 0]
        pertSymbols = {"lam11", "lam12", "lam12p"}
        allSymbols = ["lam11", "lam12", "lam12p"]

        self.assertEqual(reference, bIsPerturbative(source, pertSymbols, allSymbols))

    def test_bIsPerturbativeFalse(self):
        reference = False
        source = [-12.57, 0, 0]
        pertSymbols = {"lam11", "lam12", "lam12p"}
        allSymbols = ["lam11", "lam12", "lam12p"]

        self.assertEqual(reference, bIsPerturbative(source, pertSymbols, allSymbols))
