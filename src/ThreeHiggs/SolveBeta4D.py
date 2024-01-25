import numpy as np
import matplotlib.pylab as plt
from math import pi as Pi
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline

def beta_functions(coupling, mu):
    #Unpacking the coupling array, to reduce lines could be done with g1, g2,...Mu3sq = coupling[:]
    g1 = coupling[0] #initial[key,value]
    g2 = coupling[1]
    g3 = coupling[2]
    lambda11 = coupling[3]
    lambda12p = coupling[4]
    lambda12 = coupling[5]
    lambda1Im = coupling[6]
    lambda1Re = coupling[7]
    lambda22 = coupling[8]
    lambda23p = coupling[9]
    lambda23 = coupling[10]
    lambda2Im = coupling[11]
    lambda2Re = coupling[12]
    lambda31p = coupling[13]
    lambda31 = coupling[14]
    lambda33 = coupling[15]
    lambda3Im = coupling[16]
    lambda3Re = coupling[17]
    yt3 = coupling[18]
    Mu12sqIm = coupling[19]
    Mu12sqRe = coupling[20]
    Mu1sq = coupling[21]
    Mu2sq = coupling[22]
    Mu3sq = coupling[23]
    
    #Each differential equation is copy and pasted from the DRalgo file BetaFunctions4D[]//FortranForm 
    dg1 = ((43*g1**4)/(48.*Pi**2))/(2*g1)
    dg2 = ((-17*g2**4)/(48.*Pi**2))/(2*g2)
    dg3 = ((-7*g3**4)/(8.*Pi**2))/(2*g3)
    dlambda11 = (3*g1**2*(g2**2 - lambda12p) - 9*g2**2*lambda12p + 32*(lambda1Im**2 + lambda1Re**2) + 4*lambda12p*(lambda11 + 2*lambda12 + lambda12p + lambda22) + 2*lambda23p*lambda31p)/(16.*Pi**2)
    dlambda12p = (3*g1**2*(g2**2 - lambda12p) - 9*g2**2*lambda12p + 32*(lambda1Im**2 + lambda1Re**2) + 4*lambda12p*(lambda11 + 2*lambda12 + lambda12p + lambda22) + 2*lambda23p*lambda31p)/(16.*Pi**2)
    dlambda12 = (3*g1**4 + 9*g2**4 - 36*g2**2*lambda12 - 6*g1**2*(g2**2 + 2*lambda12) + 8*(6*lambda11*lambda12 + 2*lambda12**2 + 2*lambda11*lambda12p + lambda12p**2 + 4*lambda1Im**2 + 4*lambda1Re**2 + 6*lambda12*lambda22 + 2*lambda12p*lambda22 + 2*lambda23*lambda31 + lambda23p*lambda31 + lambda23*lambda31p))/(64.*Pi**2)
    dlambda1Im = -0.0625*(3*g1**2*lambda1Im + 9*g2**2*lambda1Im - 4*lambda1Im*(lambda11 + 2*lambda12 + 3*lambda12p + lambda22) + 4*lambda2Re*lambda3Im + 4*lambda2Im*lambda3Re)/Pi**2
    dlambda1Re = (-3*g1**2*lambda1Re - 9*g2**2*lambda1Re + 4*lambda1Re*(lambda11 + 2*lambda12 + 3*lambda12p + lambda22) - 4*lambda2Im*lambda3Im + 4*lambda2Re*lambda3Re)/(16.*Pi**2)
    dlambda22 = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lambda22) - 72*g2**2*lambda22 + 8*(2*lambda12**2 + 2*lambda12*lambda12p + lambda12p**2 + 4*lambda1Im**2 + 4*lambda1Re**2 + 24*lambda22**2 + 2*lambda23**2 + 2*lambda23*lambda23p + lambda23p**2 + 4*(lambda2Im**2 + lambda2Re**2)))/(128.*Pi**2)
    dlambda23p = (3*g1**2*(g2**2 - lambda23p) - 9*g2**2*lambda23p + 6*yt3**2*lambda23p + 32*(lambda2Im**2 + lambda2Re**2) + 2*lambda12p*lambda31p + 4*lambda23p*(lambda22 + 2*lambda23 + lambda23p + lambda33))/(16.*Pi**2)
    dlambda23 = (3*g1**4 + 9*g2**4 - 36*g2**2*lambda23 - 6*g1**2*(g2**2 + 2*lambda23) + 8*(3*yt3**2*lambda23 + 6*lambda22*lambda23 + 2*lambda23**2 + 2*lambda22*lambda23p + lambda23p**2 + 4*lambda2Im**2 + 4*lambda2Re**2 + 2*lambda12*lambda31 + lambda12p*lambda31 + lambda12*lambda31p + 6*lambda23*lambda33 + 2*lambda23p*lambda33))/(64.*Pi**2)
    dlambda2Im = (lambda2Im*(-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lambda22 + 2*lambda23 + 3*lambda23p + lambda33)) - 4*(lambda1Re*lambda3Im + lambda1Im*lambda3Re))/(16.*Pi**2)
    dlambda2Re = (lambda2Re*(-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lambda22 + 2*lambda23 + 3*lambda23p + lambda33)) - 4*lambda1Im*lambda3Im + 4*lambda1Re*lambda3Re)/(16.*Pi**2)
    dlambda31p = (2*lambda12p*lambda23p + 3*g1**2*(g2**2 - lambda31p) + lambda31p*(-9*g2**2 + 6*yt3**2 + 4*(lambda11 + 2*lambda31 + lambda31p + lambda33)) + 32*(lambda3Im**2 + lambda3Re**2))/(16.*Pi**2)
    dlambda31 = (3*g1**4 + 9*g2**4 - 36*g2**2*lambda31 - 6*g1**2*(g2**2 + 2*lambda31) + 8*(2*lambda12*lambda23 + lambda12p*lambda23 + lambda12*lambda23p + 3*yt3**2*lambda31 + 6*lambda11*lambda31 + 2*lambda31**2 + 2*lambda11*lambda31p + lambda31p**2 + 6*lambda31*lambda33 + 2*lambda31p*lambda33 + 4*(lambda3Im**2 + lambda3Re**2)))/(64.*Pi**2)
    dlambda33 = (3*g1**4 + 9*g2**4 + 6*g1**2*(g2**2 - 4*lambda33) - 72*g2**2*lambda33 + 96*lambda33*(yt3**2 + 2*lambda33) + 8*(-6*yt3**4 + 2*lambda23**2 + 2*lambda23*lambda23p + lambda23p**2 + 4*lambda2Im**2 + 4*lambda2Re**2 + 2*lambda31**2 + 2*lambda31*lambda31p + lambda31p**2 + 4*(lambda3Im**2 + lambda3Re**2)))/(128.*Pi**2)
    dlambda3Im = (-4*lambda1Re*lambda2Im - 4*lambda1Im*lambda2Re + (-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lambda11 + 2*lambda31 + 3*lambda31p + lambda33))*lambda3Im)/(16.*Pi**2)
    dlambda3Re = (-4*lambda1Im*lambda2Im + 4*lambda1Re*lambda2Re + (-3*g1**2 - 9*g2**2 + 6*yt3**2 + 4*(lambda11 + 2*lambda31 + 3*lambda31p + lambda33))*lambda3Re)/(16.*Pi**2)
    dyt3 = (yt3*(-17*g1**2 - 27*g2**2 - 96*g3**2 + 54*yt3**2))/(192.*Pi**2)
    dmu12SqIm = (-3*g1**2*Mu12sqIm - 9*g2**2*Mu12sqIm + 4*(lambda12 + 2*lambda12p - 6*lambda1Re)*Mu12sqIm + 24*lambda1Im*Mu12sqRe)/(32.*Pi**2)
    dmu12SqRe = (24*lambda1Im*Mu12sqIm + (-3*g1**2 - 9*g2**2 + 4*lambda12 + 8*lambda12p + 24*lambda1Re)*Mu12sqRe)/(32.*Pi**2)
    dmu1Sq = (-3*(g1**2 + 3*g2**2 - 8*lambda11)*Mu1sq + 8*lambda12*Mu2sq + 4*lambda12p*Mu2sq + 4*(2*lambda31 + lambda31p)*Mu3sq)/(32.*Pi**2)
    dmu2Sq = (8*lambda12*Mu1sq + 4*lambda12p*Mu1sq - 3*(g1**2 + 3*g2**2 - 8*lambda22)*Mu2sq + 4*(2*lambda23 + lambda23p)*Mu3sq)/(32.*Pi**2)
    dmu3Sq = (8*lambda31*Mu1sq + 4*lambda31p*Mu1sq + 4*(2*lambda23 + lambda23p)*Mu2sq - 3*(g1**2 + 3*g2**2 - 4*yt3**2 - 8*lambda33)*Mu3sq)/(32.*Pi**2)
    
    return [dg1,dg2, dg3, dlambda11, dlambda12p, dlambda12, dlambda1Im, dlambda1Re, dlambda22, dlambda23p, dlambda23, dlambda2Im, dlambda2Re, dlambda31p, 
            dlambda31, dlambda33, dlambda3Im, dlambda3Re, dyt3, dmu12SqIm, dmu12SqRe, dmu1Sq, dmu2Sq, dmu3Sq]


#[dg1,dg2, g3, lambda11, lambda12p, lambda12, lambda1Im, lambda1Re, lambda22, dlambda23p, dlambda23, dlambda2Im, dlambda2Re, dlambda31p, dlambda31, dlambda33, dlambda3Im, dlambda3Re, dyt3, dmu12SqIm, dmu12SqRe, dmu1Sq, dmu2Sq, dmu3Sq]
mu_initital = [np.sqrt(0.15), np.sqrt(0.4), np.sqrt(1.9), 0.1, 0.1, 0.1, 0.0, 0.1, 0.1, 1, 1.1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.61, 0.71, 1.0, 0.91, 1.01, 0.92, 0.82, 0.75]

muRange = np.linspace(90, 700, 300)
mubarRange = np.log(muRange)

solution = odeint(beta_functions, mu_initital, mubarRange)
Solution_trans = np.transpose(solution)
test1 = CubicSpline(mubarRange, Solution_trans[0], extrapolate=False)

plt.plot(mubarRange, Solution_trans[0], label ='U(1)')
plt.plot(mubarRange, test1(mubarRange), label ='U(1) spline')
#plt.plot(mubar_array, solution_trans[1], label = 'SU(2)')
#plt.plot(mubar_array, solution_trans[2], label = 'SU(3)')
#plt.plot(mubar_array, solution_trans[18], label = 'yt')
#plt.legend(loc ='best')
#plt.title("Running of the gauge couplings and top yukawa")
#plt.xlabel('$log(\mu)$')
#plt.ylabel('coupling')
plt.show()

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

Running_coupling_interp_dict = {
    "dg1": dg1_interp,
    "dg2": dg2_interp,
    "dg3": dg3_interp,
    "dlam11": dlam11_interp,
    "dlam12p": dlam12p_interp,
    "dlam12": dlam12_interp,
    "dlam1Im": dlam1Im_interp,
    "dlam1Re": dlam1Re_interp,
    "dlam22": dlam22_interp,
    "dlam23p": dlam23p_interp,
    "dlam23": dlam23_interp,
    "dlam2Im": dlam2Im_interp,
    "dlam2Re": dlam2Re_interp,
    "dlam31p": dlam31p_interp,
    "dlam31": dlam31_interp,
    "dlam33": dlam33_interp,
    "dlam3Im": dlam3Im_interp,
    "dlam3Re": dlam3Re_interp,
    "dyt3": dyt3_interp,
    "dmu12sqIm": dmu12sqIm_interp,
    "dmu12sqRe": dmu12sqRe_interp,
    "dmu1Sq":dmu1sq_interp ,
    "dmu2sq": dmu2sq_interp,
    "dmu3sq": dmu3sq_interp,
    }

print (Running_coupling_interp_dict["dg1"](150))