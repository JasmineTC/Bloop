import numpy as np
## This is for the hard coded - remove later
from numpy import pi as Pi
import scipy

class BetaFunctions4D():
    def __init__(self):
        self.pi16 = 16.*Pi**2
        
    def unpackParamDict(self, params: dict[str, float]) -> np.ndarray:
        """Puts a 3HDM parameter dict in array format that odeint understands.
        Also produces a name <-> index mapping for easier access to the params in beta function expressions."""
        indexMapping = {}
        keyMapping = []
        paramsList = []
        for idx, (key, val) in enumerate(params.items()):
            ## Do NOT include the RG scale in this
            if (key == "RGScale"): 
                continue

            paramsList.append(val)
            indexMapping[key] = idx
            keyMapping.append(key)
        self._keyMapping = keyMapping
        return paramsList, indexMapping
        
    def constructSplineDict(self, muRange, initialParams: dict[str, float]):
        ## Need to unpack the param dict for odeint
        paramsList, indexMapping = self.unpackParamDict(initialParams)
        solution = np.transpose( scipy.integrate.odeint(self.hardCodeBetaFunction, 
                                        paramsList, 
                                        muRange,
                                        args = ( indexMapping, ) ) )

        
        ## This makes it non-trivial to make this a frozen class
        self.interpDict = {}
        for i, key in enumerate(self._keyMapping):
            self.interpDict[key] =  scipy.interpolate.CubicSpline(muRange, solution[i], extrapolate = False)
        
    def hardCodeBetaFunction(self, InitialConditions: np.ndarray, mu: float, indexMapping) -> np.ndarray:

        ## Pick params from the input array since they are appear as hardcoded symbols below
        g1 = InitialConditions[ indexMapping["g1"] ]
        g2 = InitialConditions[ indexMapping["g2"] ]
        g3 = InitialConditions[ indexMapping["g3"] ]
        lam11 = InitialConditions[ indexMapping["lam11"] ]
        lam12p = InitialConditions[ indexMapping["lam12p"] ]
        lam12 = InitialConditions[ indexMapping["lam12"] ]
        lam1Im = InitialConditions[ indexMapping["lam1Im"] ]
        lam1Re = InitialConditions[ indexMapping["lam1Re"] ]
        lam22 = InitialConditions[ indexMapping["lam22"] ]
        lam23p = InitialConditions[ indexMapping["lam23p"] ]
        lam23 = InitialConditions[ indexMapping["lam23"] ]
        lam2Im = InitialConditions[ indexMapping["lam2Im"] ]
        lam2Re = InitialConditions[ indexMapping["lam2Re"] ]
        lam31p = InitialConditions[ indexMapping["lam31p"] ]
        lam31 = InitialConditions[ indexMapping["lam31"] ]
        lam33 = InitialConditions[ indexMapping["lam33"] ]
        lam3Im = InitialConditions[ indexMapping["lam3Im"] ]
        lam3Re = InitialConditions[ indexMapping["lam3Re"] ]
        yt3 = InitialConditions[ indexMapping["yt3"] ]
        mu12sqIm = InitialConditions[ indexMapping["mu12sqIm"] ]
        mu12sqRe = InitialConditions[ indexMapping["mu12sqRe"] ]
        mu1sq = InitialConditions[ indexMapping["mu1sq"] ]
        mu2sq = InitialConditions[ indexMapping["mu2sq"] ]
        mu3sq = InitialConditions[ indexMapping["mu3sq"] ]
        
        ## Each differential equation is copy and pasted from the DRalgo file - BetaFunctions4D[]//FortranForm, except
        ## Except for the gauge couplings which need to be divided by 2 g as DRalgo gives the beta function as dg^2/dmu and odeint assumes dg/dmu
        betaArray = np.zeros(len(InitialConditions))
        betaArray[ indexMapping["g1"] ] = (43*g1**3)/(self.pi16*6)
        betaArray[ indexMapping["g2"] ] = (-17*g2**3)/(self.pi16*6)
        betaArray[ indexMapping["g3"] ] = (-7*g3**3)/(self.pi16)
        betaArray[ indexMapping["lam11"] ] = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lam11) - 72*g2**2*lam11 + 8*(24*lam11**2 + 2*lam12**2 + 2*lam12*lam12p + lam12p**2 + 4*lam1Im**2 + 4*lam1Re**2 + 2*lam31**2 + 2*lam31*lam31p + lam31p**2 + 4*(lam3Im**2 + lam3Re**2)))/(8*self.pi16)
        betaArray[ indexMapping["lam12p"] ] = (3*g1**2*(g2**2 - lam12p) - 9*g2**2*lam12p + 32*(lam1Im**2 + lam1Re**2) + 4*lam12p*(lam11 + 2*lam12 + lam12p + lam22) + 2*lam23p*lam31p)/(self.pi16)
        betaArray[ indexMapping["lam12"] ] = (3*g1**4 + 9*g2**4 - 36*g2**2*lam12 - 6*g1**2*(g2**2 + 2*lam12) + 8*(6*lam11*lam12 + 2*lam12**2 + 2*lam11*lam12p + lam12p**2 + 4*lam1Im**2 + 4*lam1Re**2 + 6*lam12*lam22 + 2*lam12p*lam22 + 2*lam23*lam31 + lam23p*lam31 + lam23*lam31p))/(4*self.pi16)
        betaArray[ indexMapping["lam1Im"] ] = -0.0625*(3*g1**2*lam1Im + 9*g2**2*lam1Im - 4*lam1Im*(lam11 + 2*lam12 + 3*lam12p + lam22) + 4*lam2Re*lam3Im + 4*lam2Im*lam3Re)/Pi**2
        betaArray[ indexMapping["lam1Re"] ] = (-3*g1**2*lam1Re - 9*g2**2*lam1Re + 4*lam1Re*(lam11 + 2*lam12 + 3*lam12p + lam22) - 4*lam2Im*lam3Im + 4*lam2Re*lam3Re)/(self.pi16)
        betaArray[ indexMapping["lam22"] ] = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lam22) - 72*g2**2*lam22 + 8*(2*lam12**2 + 2*lam12*lam12p + lam12p**2 + 4*lam1Im**2 + 4*lam1Re**2 + 24*lam22**2 + 2*lam23**2 + 2*lam23*lam23p + lam23p**2 + 4*(lam2Im**2 + lam2Re**2)))/(8*self.pi16)
        betaArray[ indexMapping["lam23p"] ] = (3*g1**2*(g2**2 - lam23p) - 9*g2**2*lam23p + 6*yt3**2*lam23p + 32*(lam2Im**2 + lam2Re**2) + 2*lam12p*lam31p + 4*lam23p*(lam22 + 2*lam23 + lam23p + lam33))/(self.pi16)
        betaArray[ indexMapping["lam23"] ] = (3*g1**4 + 9*g2**4 - 36*g2**2*lam23 - 6*g1**2*(g2**2 + 2*lam23) + 8*(3*yt3**2*lam23 + 6*lam22*lam23 + 2*lam23**2 + 2*lam22*lam23p + lam23p**2 + 4*lam2Im**2 + 4*lam2Re**2 + 2*lam12*lam31 + lam12p*lam31 + lam12*lam31p + 6*lam23*lam33 + 2*lam23p*lam33))/(4.*self.pi16)
        betaArray[ indexMapping["lam2Im"] ] = (lam2Im*(-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam22 + 2*lam23 + 3*lam23p + lam33)) - 4*(lam1Re*lam3Im + lam1Im*lam3Re))/(self.pi16)
        betaArray[ indexMapping["lam2Re"] ] = (lam2Re*(-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam22 + 2*lam23 + 3*lam23p + lam33)) - 4*lam1Im*lam3Im + 4*lam1Re*lam3Re)/(self.pi16)
        betaArray[ indexMapping["lam31p"] ] = (2*lam12p*lam23p + 3*g1**2*(g2**2 - lam31p) + lam31p*(-9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + lam31p + lam33)) + 32*(lam3Im**2 + lam3Re**2))/(self.pi16)
        betaArray[ indexMapping["lam31"] ] = (3*g1**4 + 9*g2**4 - 36*g2**2*lam31 - 6*g1**2*(g2**2 + 2*lam31) + 8*(2*lam12*lam23 + lam12p*lam23 + lam12*lam23p + 3*yt3**2*lam31 + 6*lam11*lam31 + 2*lam31**2 + 2*lam11*lam31p + lam31p**2 + 6*lam31*lam33 + 2*lam31p*lam33 + 4*(lam3Im**2 + lam3Re**2)))/(4.*self.pi16)
        betaArray[ indexMapping["lam33"] ] = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lam33) - 72*g2**2*lam33 + 96*lam33*(yt3**2 + 2*lam33) + 8*(-6*yt3**4 + 2*lam23**2 + 2*lam23*lam23p + lam23p**2 + 4*lam2Im**2 + 4*lam2Re**2 + 2*lam31**2 + 2*lam31*lam31p + lam31p**2 + 4*(lam3Im**2 + lam3Re**2)))/(8.*self.pi16)
        betaArray[ indexMapping["lam3Im"] ] = (-4*lam1Re*lam2Im - 4*lam1Im*lam2Re + (-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + 3*lam31p + lam33))*lam3Im)/(self.pi16)
        betaArray[ indexMapping["lam3Re"] ] = (-4*lam1Im*lam2Im + 4*lam1Re*lam2Re + (-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + 3*lam31p + lam33))*lam3Re)/(self.pi16)
        betaArray[ indexMapping["yt3"] ] = (yt3*(-17*g1**2 - 27*g2**2 - 96*g3**2 + 54*yt3**2))/(12.*self.pi16)
        betaArray[ indexMapping["mu12sqIm"] ] = (-3*g1**2*mu12sqIm - 9*g2**2*mu12sqIm + 4*(lam12 + 2*lam12p - 6*lam1Re)*mu12sqIm + 24*lam1Im*mu12sqRe)/(2.*self.pi16)
        betaArray[ indexMapping["mu12sqRe"] ] = (24*lam1Im*mu12sqIm + (-3*g1**2 - 9*g2**2 + 4*lam12 + 8*lam12p + 24*lam1Re)*mu12sqRe)/(2.*self.pi16)
        betaArray[ indexMapping["mu1sq"] ] = (-3*(g1**2 + 3*g2**2 - 8*lam11)*mu1sq + 8*lam12*mu2sq + 4*lam12p*mu2sq + 4*(2*lam31 + lam31p)*mu3sq)/(2.*self.pi16)
        betaArray[ indexMapping["mu2sq"] ] = (8*lam12*mu1sq + 4*lam12p*mu1sq - 3*(g1**2 + 3*g2**2 - 8*lam22)*mu2sq + 4*(2*lam23 + lam23p)*mu3sq)/(2.*self.pi16)
        betaArray[ indexMapping["mu3sq"] ] = (8*lam31*mu1sq + 4*lam31p*mu1sq + 4*(2*lam23 + lam23p)*mu2sq - 3*(g1**2 + 3*g2**2 - 4*yt3**2 - 8*lam33)*mu3sq)/(2.*self.pi16)
        
        return betaArray/mu

        
    def runCoupling(self, muEvaulate: float):
        
        runCoupling = {}
        ##Loop over the dict, for each key compute the spline at muEvaulate, ignore the RGScale in the dict
        for key in self._keyMapping:
            ##RHS is an array of a single number, so add item() to extract that single number
            runCoupling[key] = self.interpDict[key](muEvaulate).item()
        return runCoupling