import numpy as np
from numpy import pi as Pi

from scipy.integrate import odeint
from scipy.interpolate import CubicSpline

class BetaFunctions4D():
    """ Order of computations:
        1. Call SolveBetaFunction given a dictionary containing the coupling name and initial condition (the value of the coupling at 0T at the z mass)
        2. Call UnpackParamDict to generate an array of initial conditions (needed for odeint), also stores a mapping for coupling <-> index
           so that dictionaries/arrays can be generated with a consistent order
        3. SolveBetaFunction takes the initial condition list from UnpackParamDict gives it to odeint to solve the beta function up to some scale
        4. Interpolations of the odeint data are made and stored
        5. Call RunCoupling to evaluate the interpolations and return a dict
    """
    def __init__(self, muRange, initialParams: dict[str, float]):
        self.muRange = muRange
        self.pi16 = 16.*Pi**2

        ## Need to unpack the param dict for odeint
        self.paramsList = self.unpackParamDict(initialParams)
        
    def constructSplineDict(self):
        solution = odeint(self.hardCodeBetaFunction, self.paramsList, self.muRange)

        solution = np.transpose(solution)
        
        ## This makes it non-trivial to make this a frozen class
        self.interpDict = {}
        for i, key in enumerate(self._keyMapping):
            self.interpDict[key] =  CubicSpline(self.muRange, solution[i], extrapolate = False)
        
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

        self._indexMapping = indexMapping
        self._keyMapping = keyMapping
        return paramsList    
        
    def hardCodeBetaFunction(self, InitialConditions: np.ndarray, mu: float) -> np.ndarray:
        ## Pick params from the input array since they are appear as hardcoded symbols below
        g1 = InitialConditions[ self._indexMapping["g1"] ]
        g2 = InitialConditions[ self._indexMapping["g2"] ]
        g3 = InitialConditions[ self._indexMapping["g3"] ]
        lam11 = InitialConditions[ self._indexMapping["lam11"] ]
        lam12p = InitialConditions[ self._indexMapping["lam12p"] ]
        lam12 = InitialConditions[ self._indexMapping["lam12"] ]
        lam1Im = InitialConditions[ self._indexMapping["lam1Im"] ]
        lam1Re = InitialConditions[ self._indexMapping["lam1Re"] ]
        lam22 = InitialConditions[ self._indexMapping["lam22"] ]
        lam23p = InitialConditions[ self._indexMapping["lam23p"] ]
        lam23 = InitialConditions[ self._indexMapping["lam23"] ]
        lam2Im = InitialConditions[ self._indexMapping["lam2Im"] ]
        lam2Re = InitialConditions[ self._indexMapping["lam2Re"] ]
        lam31p = InitialConditions[ self._indexMapping["lam31p"] ]
        lam31 = InitialConditions[ self._indexMapping["lam31"] ]
        lam33 = InitialConditions[ self._indexMapping["lam33"] ]
        lam3Im = InitialConditions[ self._indexMapping["lam3Im"] ]
        lam3Re = InitialConditions[ self._indexMapping["lam3Re"] ]
        yt3 = InitialConditions[ self._indexMapping["yt3"] ]
        mu12sqIm = InitialConditions[ self._indexMapping["mu12sqIm"] ]
        mu12sqRe = InitialConditions[ self._indexMapping["mu12sqRe"] ]
        mu1sq = InitialConditions[ self._indexMapping["mu1sq"] ]
        mu2sq = InitialConditions[ self._indexMapping["mu2sq"] ]
        mu3sq = InitialConditions[ self._indexMapping["mu3sq"] ]
        
        ## Each differential equation is copy and pasted from the DRalgo file - BetaFunctions4D[]//FortranForm, except
        ## For a tiny speed up 16*pi**2 is computed in init and then called - another tiny speed up,  get rid of all the divisons
        ## Note this does lead to floating point errors
        ## Except for the gauge couplings which need to be divided by 2 g as DRalgo gives the beta function as dg^2/dmu and odeint assumes dg/dmu
        ## TODO Is it better to work with g^2?
        dg1 = (43*g1**3)/(self.pi16*6)
        dg2 = (-17*g2**3)/(self.pi16*6)
        dg3 = (-7*g3**3)/(self.pi16)
        dlam11 = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lam11) - 72*g2**2*lam11 + 8*(24*lam11**2 + 2*lam12**2 + 2*lam12*lam12p + lam12p**2 + 4*lam1Im**2 + 4*lam1Re**2 + 2*lam31**2 + 2*lam31*lam31p + lam31p**2 + 4*(lam3Im**2 + lam3Re**2)))/(8*self.pi16)
        dlam12p = (3*g1**2*(g2**2 - lam12p) - 9*g2**2*lam12p + 32*(lam1Im**2 + lam1Re**2) + 4*lam12p*(lam11 + 2*lam12 + lam12p + lam22) + 2*lam23p*lam31p)/(self.pi16)
        dlam12 = (3*g1**4 + 9*g2**4 - 36*g2**2*lam12 - 6*g1**2*(g2**2 + 2*lam12) + 8*(6*lam11*lam12 + 2*lam12**2 + 2*lam11*lam12p + lam12p**2 + 4*lam1Im**2 + 4*lam1Re**2 + 6*lam12*lam22 + 2*lam12p*lam22 + 2*lam23*lam31 + lam23p*lam31 + lam23*lam31p))/(4*self.pi16)
        dlam1Im = -0.0625*(3*g1**2*lam1Im + 9*g2**2*lam1Im - 4*lam1Im*(lam11 + 2*lam12 + 3*lam12p + lam22) + 4*lam2Re*lam3Im + 4*lam2Im*lam3Re)/Pi**2
        dlam1Re = (-3*g1**2*lam1Re - 9*g2**2*lam1Re + 4*lam1Re*(lam11 + 2*lam12 + 3*lam12p + lam22) - 4*lam2Im*lam3Im + 4*lam2Re*lam3Re)/(self.pi16)
        dlam22 = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lam22) - 72*g2**2*lam22 + 8*(2*lam12**2 + 2*lam12*lam12p + lam12p**2 + 4*lam1Im**2 + 4*lam1Re**2 + 24*lam22**2 + 2*lam23**2 + 2*lam23*lam23p + lam23p**2 + 4*(lam2Im**2 + lam2Re**2)))/(8*self.pi16)
        dlam23p = (3*g1**2*(g2**2 - lam23p) - 9*g2**2*lam23p + 6*yt3**2*lam23p + 32*(lam2Im**2 + lam2Re**2) + 2*lam12p*lam31p + 4*lam23p*(lam22 + 2*lam23 + lam23p + lam33))/(self.pi16)
        dlam23 = (3*g1**4 + 9*g2**4 - 36*g2**2*lam23 - 6*g1**2*(g2**2 + 2*lam23) + 8*(3*yt3**2*lam23 + 6*lam22*lam23 + 2*lam23**2 + 2*lam22*lam23p + lam23p**2 + 4*lam2Im**2 + 4*lam2Re**2 + 2*lam12*lam31 + lam12p*lam31 + lam12*lam31p + 6*lam23*lam33 + 2*lam23p*lam33))/(4.*self.pi16)
        dlam2Im = (lam2Im*(-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam22 + 2*lam23 + 3*lam23p + lam33)) - 4*(lam1Re*lam3Im + lam1Im*lam3Re))/(self.pi16)
        dlam2Re = (lam2Re*(-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam22 + 2*lam23 + 3*lam23p + lam33)) - 4*lam1Im*lam3Im + 4*lam1Re*lam3Re)/(self.pi16)
        dlam31p = (2*lam12p*lam23p + 3*g1**2*(g2**2 - lam31p) + lam31p*(-9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + lam31p + lam33)) + 32*(lam3Im**2 + lam3Re**2))/(self.pi16)
        dlam31 = (3*g1**4 + 9*g2**4 - 36*g2**2*lam31 - 6*g1**2*(g2**2 + 2*lam31) + 8*(2*lam12*lam23 + lam12p*lam23 + lam12*lam23p + 3*yt3**2*lam31 + 6*lam11*lam31 + 2*lam31**2 + 2*lam11*lam31p + lam31p**2 + 6*lam31*lam33 + 2*lam31p*lam33 + 4*(lam3Im**2 + lam3Re**2)))/(4.*self.pi16)
        dlam33 = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lam33) - 72*g2**2*lam33 + 96*lam33*(yt3**2 + 2*lam33) + 8*(-6*yt3**4 + 2*lam23**2 + 2*lam23*lam23p + lam23p**2 + 4*lam2Im**2 + 4*lam2Re**2 + 2*lam31**2 + 2*lam31*lam31p + lam31p**2 + 4*(lam3Im**2 + lam3Re**2)))/(8.*self.pi16)
        dlam3Im = (-4*lam1Re*lam2Im - 4*lam1Im*lam2Re + (-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + 3*lam31p + lam33))*lam3Im)/(self.pi16)
        dlam3Re = (-4*lam1Im*lam2Im + 4*lam1Re*lam2Re + (-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + 3*lam31p + lam33))*lam3Re)/(self.pi16)
        dyt3 = (yt3*(-17*g1**2 - 27*g2**2 - 96*g3**2 + 54*yt3**2))/(12.*self.pi16)
        dmu12sqIm = (-3*g1**2*mu12sqIm - 9*g2**2*mu12sqIm + 4*(lam12 + 2*lam12p - 6*lam1Re)*mu12sqIm + 24*lam1Im*mu12sqRe)/(2.*self.pi16)
        dmu12sqRe = (24*lam1Im*mu12sqIm + (-3*g1**2 - 9*g2**2 + 4*lam12 + 8*lam12p + 24*lam1Re)*mu12sqRe)/(2.*self.pi16)
        dmu1sq = (-3*(g1**2 + 3*g2**2 - 8*lam11)*mu1sq + 8*lam12*mu2sq + 4*lam12p*mu2sq + 4*(2*lam31 + lam31p)*mu3sq)/(2.*self.pi16)
        dmu2sq = (8*lam12*mu1sq + 4*lam12p*mu1sq - 3*(g1**2 + 3*g2**2 - 8*lam22)*mu2sq + 4*(2*lam23 + lam23p)*mu3sq)/(2.*self.pi16)
        dmu3sq = (8*lam31*mu1sq + 4*lam31p*mu1sq + 4*(2*lam23 + lam23p)*mu2sq - 3*(g1**2 + 3*g2**2 - 4*yt3**2 - 8*lam33)*mu3sq)/(2.*self.pi16)
        
        
        betaArray = np.zeros(len(InitialConditions))
        ##Return an array of the beta function, divided by mu as the beta function is mu dg/dmu and this seemed easier than working with log(mu)
        betaArray[ self._indexMapping["g1"] ] = dg1
        betaArray[ self._indexMapping["g2"] ] = dg2
        betaArray[ self._indexMapping["g3"] ] = dg3
        betaArray[ self._indexMapping["lam11"] ] = dlam11
        betaArray[ self._indexMapping["lam12p"] ] = dlam12p
        betaArray[ self._indexMapping["lam12"] ] = dlam12
        betaArray[ self._indexMapping["lam1Im"] ] = dlam1Im
        betaArray[ self._indexMapping["lam1Re"] ] = dlam1Re
        betaArray[ self._indexMapping["lam22"] ] = dlam22
        betaArray[ self._indexMapping["lam23p"] ] = dlam23p
        betaArray[ self._indexMapping["lam23"] ] = dlam23
        betaArray[ self._indexMapping["lam2Im"] ] = dlam2Im
        betaArray[ self._indexMapping["lam2Re"] ] = dlam2Re
        betaArray[ self._indexMapping["lam31p"] ] = dlam31p
        betaArray[ self._indexMapping["lam31"] ] = dlam31
        betaArray[ self._indexMapping["lam33"] ] = dlam33
        betaArray[ self._indexMapping["lam3Im"] ] = dlam3Im
        betaArray[ self._indexMapping["lam3Re"] ] = dlam3Re
        betaArray[ self._indexMapping["yt3"] ] = dyt3
        betaArray[ self._indexMapping["mu12sqIm"] ] = dmu12sqIm
        betaArray[ self._indexMapping["mu12sqRe"] ] = dmu12sqRe
        betaArray[ self._indexMapping["mu1sq"] ] = dmu1sq
        betaArray[ self._indexMapping["mu2sq"] ] = dmu2sq
        betaArray[ self._indexMapping["mu3sq"] ] = dmu3sq
        
        return betaArray/mu

    def runCoupling(self, muEvaulate: float):
        
        runCoupling = {}
        ##Loop over the dict, for each key compute the spline at muEvaulate, ignore the RGScale in the dict
        for key in self._keyMapping:
            ##RHS is an array of a single number, so add item() to extract that single number
            runCoupling[key] = self.interpDict[key](muEvaulate).item()
        return runCoupling
