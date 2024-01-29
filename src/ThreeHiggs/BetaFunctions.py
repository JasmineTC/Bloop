import numpy as np
import numpy.typing as npt
from typing import Tuple
import pathlib
from math import pi as Pi
import scipy.integrate
from scipy.interpolate import CubicSpline



pathToCurrentFile = pathlib.Path(__file__).parent.resolve()

BetaFunctions4DFile = str(pathToCurrentFile) + "/Data/BetaFunctions4D/BetaFunctions4D.txt"

"""Solve beta function
"""
class BetaFunctions4D():
    
    def __init__(self):
        None


    def HardCodeBetaFunction(self, InitialConditions: np.ndarray, mubar: float) -> np.ndarray:
        ## Let's run SU3 coupling too because why not

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
        
        ##Each differential equation is copy and pasted from the DRalgo file BetaFunctions4D[]//FortranForm,
        ## Except for the gauge couplings which need to be divided by 2 g as DRalgo gives the beta function as dg^2/dmu and odeint assumes dg/dmu
        ##TODO Is it better to work with g^2?
        dg1 = ((43*g1**3)/(48.*Pi**2)) / 2.
        dg2 = ((-17*g2**3)/(48.*Pi**2)) / 2.
        dg3 = ((-7*g3**3)/(8.*Pi**2)) / 2.
        dlam11 = (3*g1**2*(g2**2 - lam12p) - 9*g2**2*lam12p + 32*(lam1Im**2 + lam1Re**2) + 4*lam12p*(lam11 + 2*lam12 + lam12p + lam22) + 2*lam23p*lam31p)/(16.*Pi**2)
        dlam12p = (3*g1**2*(g2**2 - lam12p) - 9*g2**2*lam12p + 32*(lam1Im**2 + lam1Re**2) + 4*lam12p*(lam11 + 2*lam12 + lam12p + lam22) + 2*lam23p*lam31p)/(16.*Pi**2)
        dlam12 = (3*g1**4 + 9*g2**4 - 36*g2**2*lam12 - 6*g1**2*(g2**2 + 2*lam12) + 8*(6*lam11*lam12 + 2*lam12**2 + 2*lam11*lam12p + lam12p**2 + 4*lam1Im**2 + 4*lam1Re**2 + 6*lam12*lam22 + 2*lam12p*lam22 + 2*lam23*lam31 + lam23p*lam31 + lam23*lam31p))/(64.*Pi**2)
        dlam1Im = -0.0625*(3*g1**2*lam1Im + 9*g2**2*lam1Im - 4*lam1Im*(lam11 + 2*lam12 + 3*lam12p + lam22) + 4*lam2Re*lam3Im + 4*lam2Im*lam3Re)/Pi**2
        dlam1Re = (-3*g1**2*lam1Re - 9*g2**2*lam1Re + 4*lam1Re*(lam11 + 2*lam12 + 3*lam12p + lam22) - 4*lam2Im*lam3Im + 4*lam2Re*lam3Re)/(16.*Pi**2)
        dlam22 = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lam22) - 72*g2**2*lam22 + 8*(2*lam12**2 + 2*lam12*lam12p + lam12p**2 + 4*lam1Im**2 + 4*lam1Re**2 + 24*lam22**2 + 2*lam23**2 + 2*lam23*lam23p + lam23p**2 + 4*(lam2Im**2 + lam2Re**2)))/(128.*Pi**2)
        dlam23p = (3*g1**2*(g2**2 - lam23p) - 9*g2**2*lam23p + 6*yt3**2*lam23p + 32*(lam2Im**2 + lam2Re**2) + 2*lam12p*lam31p + 4*lam23p*(lam22 + 2*lam23 + lam23p + lam33))/(16.*Pi**2)
        dlam23 = (3*g1**4 + 9*g2**4 - 36*g2**2*lam23 - 6*g1**2*(g2**2 + 2*lam23) + 8*(3*yt3**2*lam23 + 6*lam22*lam23 + 2*lam23**2 + 2*lam22*lam23p + lam23p**2 + 4*lam2Im**2 + 4*lam2Re**2 + 2*lam12*lam31 + lam12p*lam31 + lam12*lam31p + 6*lam23*lam33 + 2*lam23p*lam33))/(64.*Pi**2)
        dlam2Im = (lam2Im*(-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam22 + 2*lam23 + 3*lam23p + lam33)) - 4*(lam1Re*lam3Im + lam1Im*lam3Re))/(16.*Pi**2)
        dlam2Re = (lam2Re*(-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam22 + 2*lam23 + 3*lam23p + lam33)) - 4*lam1Im*lam3Im + 4*lam1Re*lam3Re)/(16.*Pi**2)
        dlam31p = (2*lam12p*lam23p + 3*g1**2*(g2**2 - lam31p) + lam31p*(-9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + lam31p + lam33)) + 32*(lam3Im**2 + lam3Re**2))/(16.*Pi**2)
        dlam31 = (3*g1**4 + 9*g2**4 - 36*g2**2*lam31 - 6*g1**2*(g2**2 + 2*lam31) + 8*(2*lam12*lam23 + lam12p*lam23 + lam12*lam23p + 3*yt3**2*lam31 + 6*lam11*lam31 + 2*lam31**2 + 2*lam11*lam31p + lam31p**2 + 6*lam31*lam33 + 2*lam31p*lam33 + 4*(lam3Im**2 + lam3Re**2)))/(64.*Pi**2)
        dlam33 = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lam33) - 72*g2**2*lam33 + 96*lam33*(yt3**2 + 2*lam33) + 8*(-6*yt3**4 + 2*lam23**2 + 2*lam23*lam23p + lam23p**2 + 4*lam2Im**2 + 4*lam2Re**2 + 2*lam31**2 + 2*lam31*lam31p + lam31p**2 + 4*(lam3Im**2 + lam3Re**2)))/(128.*Pi**2)
        dlam3Im = (-4*lam1Re*lam2Im - 4*lam1Im*lam2Re + (-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + 3*lam31p + lam33))*lam3Im)/(16.*Pi**2)
        dlam3Re = (-4*lam1Im*lam2Im + 4*lam1Re*lam2Re + (-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lam11 + 2*lam31 + 3*lam31p + lam33))*lam3Re)/(16.*Pi**2)
        dyt3 = (yt3*(-17*g1**2 - 27*g2**2 - 96*g3**2 + 54*yt3**2))/(192.*Pi**2)
        dmu12sqIm = (-3*g1**2*mu12sqIm - 9*g2**2*mu12sqIm + 4*(lam12 + 2*lam12p - 6*lam1Re)*mu12sqIm + 24*lam1Im*mu12sqRe)/(32.*Pi**2)
        dmu12sqRe = (24*lam1Im*mu12sqIm + (-3*g1**2 - 9*g2**2 + 4*lam12 + 8*lam12p + 24*lam1Re)*mu12sqRe)/(32.*Pi**2)
        dmu1sq = (-3*(g1**2 + 3*g2**2 - 8*lam11)*mu1sq + 8*lam12*mu2sq + 4*lam12p*mu2sq + 4*(2*lam31 + lam31p)*mu3sq)/(32.*Pi**2)
        dmu2sq = (8*lam12*mu1sq + 4*lam12p*mu1sq - 3*(g1**2 + 3*g2**2 - 8*lam22)*mu2sq + 4*(2*lam23 + lam23p)*mu3sq)/(32.*Pi**2)
        dmu3sq = (8*lam31*mu1sq + 4*lam31p*mu1sq + 4*(2*lam23 + lam23p)*mu2sq - 3*(g1**2 + 3*g2**2 - 4*yt3**2 - 8*lam33)*mu3sq)/(32.*Pi**2)

        ## Probably should use the index mapping here too to avoid errors

        return np.asanyarray([dg1, dg2, dg3, dlam11, dlam12p, dlam12, dlam1Im, dlam1Re, dlam22, dlam23p, dlam23, dlam2Im, dlam2Re, dlam31p, 
                dlam31, dlam33, dlam3Im, dlam3Re, dyt3, dmu12sqIm, dmu12sqRe, dmu1sq, dmu2sq, dmu3sq])

    def SolveBetaFunction(self, initialParams: dict[str, float],  muRange: npt.ArrayLike = np.arange(91.1876, 700., 2.)) -> dict[str, float]:##what type is a scipy spline?
        ##muRange used to run the energy scale from the Z mass to around 700 in steps of 2 

        ##TODO Does anything need to be done with the RG scale in Initital conditions?
        ##taking the log so that we can work with mubar in the derivative
        muBarRange = np.log(muRange)

        ## Need to unpack the param dict for odeint
        paramsList = self._unpackParamDict(initialParams)
        
        Solution = scipy.integrate.odeint(self.HardCodeBetaFunction, paramsList, muBarRange)
        ##To make solution slightly nicer to work with we take the transpose so that each coupling is inside its own array
        Solution_trans = np.transpose(Solution)
        dg1_interp = CubicSpline(muRange, Solution_trans[0], extrapolate=False)
        dg2_interp = CubicSpline(muRange, Solution_trans[1], extrapolate=False)
        dg3_interp = CubicSpline(muRange, Solution_trans[2], extrapolate=False)
        dlam11_interp = CubicSpline(muRange, Solution_trans[3], extrapolate=False)
        dlam12p_interp = CubicSpline(muRange, Solution_trans[4], extrapolate=False)
        dlam12_interp = CubicSpline(muRange, Solution_trans[5], extrapolate=False)
        dlam1Im_interp = CubicSpline(muRange, Solution_trans[6], extrapolate=False)
        dlam1Re_interp = CubicSpline(muRange, Solution_trans[7], extrapolate=False)
        dlam22_interp = CubicSpline(muRange, Solution_trans[8], extrapolate=False)
        dlam23p_interp = CubicSpline(muRange, Solution_trans[9], extrapolate=False)
        dlam23_interp = CubicSpline(muRange, Solution_trans[10], extrapolate=False)
        dlam2Im_interp = CubicSpline(muRange, Solution_trans[11], extrapolate=False)
        dlam2Re_interp = CubicSpline(muRange, Solution_trans[12], extrapolate=False)
        dlam31p_interp = CubicSpline(muRange, Solution_trans[13], extrapolate=False)
        dlam31_interp = CubicSpline(muRange, Solution_trans[14], extrapolate=False)
        dlam33_interp = CubicSpline(muRange, Solution_trans[15], extrapolate=False)
        dlam3Im_interp = CubicSpline(muRange, Solution_trans[16], extrapolate=False)
        dlam3Re_interp = CubicSpline(muRange, Solution_trans[17], extrapolate=False)
        dyt3_interp = CubicSpline(muRange, Solution_trans[18], extrapolate=False)
        dmu12sqIm_interp = CubicSpline(muRange, Solution_trans[19], extrapolate=False)
        dmu12sqRe_interp = CubicSpline(muRange, Solution_trans[20], extrapolate=False)
        dmu1sq_interp = CubicSpline(muRange, Solution_trans[21], extrapolate=False)
        dmu2sq_interp = CubicSpline(muRange, Solution_trans[22], extrapolate=False)
        dmu3sq_interp = CubicSpline(muRange, Solution_trans[23], extrapolate=False)
        ##Hardcode a dictionary that has the beta function for each paramter, 
        ##Running_coupling_interp_dict["coupling"](mu) evaluates the interpllated beta function at the point mu

        """TODO let's not put these in dict directly. Instead, store interpolations internally
        and make a function that evaluates them at some RG scale and puts them in a dict"""
        Running_coupling_interp_dict = {
            "g1": dg1_interp,
            "g2": dg2_interp,
            "g3": dg3_interp,
            "lam11": dlam11_interp,
            "lam12p": dlam12p_interp,
            "lam12": dlam12_interp,
            "lam1Im": dlam1Im_interp,
            "lam1Re": dlam1Re_interp,
            "lam22": dlam22_interp,
            "lam23p": dlam23p_interp,
            "lam23": dlam23_interp,
            "lam2Im": dlam2Im_interp,
            "lam2Re": dlam2Re_interp,
            "lam31p": dlam31p_interp,
            "lam31": dlam31_interp,
            "lam33": dlam33_interp,
            "lam3Im": dlam3Im_interp,
            "lam3Re": dlam3Re_interp,
            "yt3": dyt3_interp,
            "mu12sqIm": dmu12sqIm_interp,
            "mu12sqRe": dmu12sqRe_interp,
            "mu1sq":dmu1sq_interp ,
            "mu2sq": dmu2sq_interp,
            "mu3sq": dmu3sq_interp,
            }
        return Running_coupling_interp_dict
    

    def _unpackParamDict(self, params: dict[str, float]) -> np.ndarray:
        """Puts a 3HDM parameter dict in array format that odeint understands.
        Also produces a name <-> index mapping for easier access to the params in beta function expressions."""

        indexMapping = {}
        paramsList = []

        idx = 0
        for key, val in params.items():
            ## Do NOT include the RG scale in this
            if (key == "RGScale"): 
                continue

            paramsList.append(val)
            indexMapping[key] = idx
            idx += 1

        paramsList = np.asanyarray(paramsList, dtype=float)

        self._indexMapping = indexMapping
        return paramsList
    



    def RunCoupling(self, muEvaulate: float):
        BetaFunctionsInterp = self.InterpolateBetaFunction()
        
        ##Easily generates the dict to return, just need to modify the value for each key
        Coupling = self.renormalizedParams
        ##When we unpacked the renormalizedParams we ignored the RGScale in the dict, this has come to bite us, 
        ##so delete it so the dict has the right length
        del Coupling["RGScale"]
        
        ##Evaluating the interp functions gives a single float in an array, so use np.ndarray to convert to float
        for i, key in enumerate (Coupling):
            Coupling[key] = np.ndarray.item(BetaFunctionsInterp[i](muEvaulate))