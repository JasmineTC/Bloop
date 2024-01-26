import numpy as np
import numpy.typing as npt
from typing import Tuple
import pathlib
from math import pi as Pi
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline



pathToCurrentFile = pathlib.Path(__file__).parent.resolve()

BetaFunctions4DFile = str(pathToCurrentFile) + "/Data/BetaFunctions4D/BetaFunctions4D.txt"

"""Solve beta function
"""
class BetaFunctions4D():
    
    def __init__(self, renormalizedParams: dict[str, float], muRange: npt.ArrayLike):
        self.renormalizedParams = renormalizedParams
        self.muRange =  muRange
        
    ##odeint needs an array to work properly, so we need to unpack the renormalized param dict
    ##This coould be done with a for loop easily, but this is safer (with for loop would be best to delete the last element (RGscale) from the dict)
    def UnpackRenormalizedParams(self):
        unpacked = np.zeros(24)
        unpacked[0] = self.renormalizedParams["lam1Re"] 
        unpacked[1]  = self.renormalizedParams["lam1Im"]
        unpacked[2]  = self.renormalizedParams["lam11"]
        unpacked[3]  = self.renormalizedParams["lam22"]
        unpacked[4]  = self.renormalizedParams["lam12"]
        unpacked[5]  = self.renormalizedParams["lam12p"]
        unpacked[6]  = self.renormalizedParams["yt3"]
        unpacked[7]  = self.renormalizedParams["g1"]
        unpacked[8]  = self.renormalizedParams["g2"]
        unpacked[9]  = self.renormalizedParams["g3"]
        unpacked[10]  = self.renormalizedParams["mu3sq"]
        unpacked[11]  = self.renormalizedParams["lam33"]
        unpacked[12]  = self.renormalizedParams["mu12sqRe"]
        unpacked[13]  = self.renormalizedParams["mu12sqIm"]
        unpacked[14]  = self.renormalizedParams["lam2Re"]
        unpacked[15]  = self.renormalizedParams["lam2Im"]
        unpacked[16]  = self.renormalizedParams["mu2sq"]
        unpacked[17]  = self.renormalizedParams["lam23"]
        unpacked[18]  = self.renormalizedParams["lam23p"]
        unpacked[19]  = self.renormalizedParams["mu1sq"]
        unpacked[20]  = self.renormalizedParams["lam3Re"]
        unpacked[21]  = self.renormalizedParams["lam3Im"]
        unpacked[22]  = self.renormalizedParams["lam31"]
        unpacked[23]  = self.renormalizedParams["lam31p"]
        return unpacked
    
    ##Takes the initial condition array, unpacks it into the parameters then computes the beta functions for that initial condition
    ##and returns the resulting array calculated points (to be feed into itself again by odeint)
    def HardCodeBetaFunction(self, initialConditions: dict[str, float], mubar: float):
        lam1Re = initialConditions[0] 
        lam1Im = initialConditions[1]
        lam11 = initialConditions[2]
        lam22 = initialConditions[3]
        lam12 = initialConditions[4]
        lam12p = initialConditions[5]
        yt3 = initialConditions[6]
        g1 = initialConditions[7]
        g2 = initialConditions[8]
        g3 = initialConditions[9]
        mu3sq = initialConditions[10]
        lam33 = initialConditions[11]
        mu12sqRe = initialConditions[12]
        mu12sqIm = initialConditions[13]
        lam2Re = initialConditions[14]
        lam2Im = initialConditions[15]
        mu2sq = initialConditions[16]
        lam23 = initialConditions[17]
        lam23p = initialConditions[18]
        mu1sq = initialConditions[19]
        lam3Re = initialConditions[20]
        lam3Im = initialConditions[21]
        lam31 = initialConditions[22]
        lam31p = initialConditions[23]
        
        ##Each differential equation is copy and pasted from the DRalgo file BetaFunctions4D[]//FortranForm,
        ## Except for the gauge couplings which need to be divided by 2 g as DRalgo gives the beta function as dg^2/dmu and odeint assumes dg/dmu
        ##TODO Is it better to work with g^2?
        dg1 = ((43*g1**4)/(48.*Pi**2))/(2*g1)
        dg2 = ((-17*g2**4)/(48.*Pi**2))/(2*g2)
        dg3 = ((-7*g3**4)/(8.*Pi**2))/(2*g3)
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
        return [dlam1Re, dlam1Im, dlam11, dlam22, dlam12, dlam12p, dyt3, dg1, dg2, dg3, dmu3sq, dlam33,
                dmu12sqRe, dmu12sqIm, dlam2Re, dlam2Im, dmu2sq, dlam23, dlam23p, dmu1sq, dlam3Re, dlam3Im, dlam31, dlam31p]  

    def SolveBetaFunction(self):
        ##Initial conditions used to store the value in renormalizedParams, len -1 to remove the renorm scale in renormalizedParams
        initialConditions = self.UnpackRenormalizedParams() 
        
        ##taking the log so that we can work with mubar in the derivative
        muBarRange = np.log(self.muRange)
        
        solution = odeint(self.HardCodeBetaFunction, initialConditions, muBarRange)
        ##To make solution slightly nicer to work with we take the transpose so that each coupling is inside its own array
        solution = np.transpose(solution)

        return solution
    
    def InterpolateBetaFunction(self):
        solution = self.SolveBetaFunction()
        
        dlam1Re = CubicSpline(self.muRange, solution[0], extrapolate=False)
        dlam1Im = CubicSpline(self.muRange, solution[1], extrapolate=False)
        dlam11 = CubicSpline(self.muRange, solution[2], extrapolate=False)
        dlam22 = CubicSpline(self.muRange, solution[3], extrapolate=False)
        dlam12 = CubicSpline(self.muRange, solution[4], extrapolate=False)
        dlam12p = CubicSpline(self.muRange, solution[5], extrapolate=False)
        dyt3 = CubicSpline(self.muRange, solution[6], extrapolate=False)
        dg1 = CubicSpline(self.muRange, solution[7], extrapolate=False)
        dg2 = CubicSpline(self.muRange, solution[8], extrapolate=False)
        dg3 = CubicSpline(self.muRange, solution[9], extrapolate=False)
        dmu3sq = CubicSpline(self.muRange, solution[10], extrapolate=False)
        dlam33 = CubicSpline(self.muRange, solution[11], extrapolate=False)
        dmu12sqRe = CubicSpline(self.muRange, solution[12], extrapolate=False)
        dmu12sqIm = CubicSpline(self.muRange, solution[13], extrapolate=False)
        dlam2Re = CubicSpline(self.muRange, solution[14], extrapolate=False)
        dlam2Im = CubicSpline(self.muRange, solution[15], extrapolate=False)
        dmu2sq = CubicSpline(self.muRange, solution[16], extrapolate=False)
        dlam23 = CubicSpline(self.muRange, solution[17], extrapolate=False)
        dlam23p = CubicSpline(self.muRange, solution[18], extrapolate=False)
        dmu1sq = CubicSpline(self.muRange, solution[19], extrapolate=False)
        dlam3Re = CubicSpline(self.muRange, solution[20], extrapolate=False)
        dlam3Im = CubicSpline(self.muRange, solution[21], extrapolate=False)
        dlam31 = CubicSpline(self.muRange, solution[22], extrapolate=False)
        dlam31p = CubicSpline(self.muRange, solution[23], extrapolate=False)
        
        return [dlam1Re, dlam1Im, dlam11, dlam22, dlam12, dlam12p, dyt3, dg1, dg2, dg3, dmu3sq, dlam33,
                dmu12sqRe, dmu12sqIm, dlam2Re, dlam2Im, dmu2sq, dlam23, dlam23p, dmu1sq, dlam3Re, dlam3Im, dlam31, dlam31p]  
        
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
        return Coupling