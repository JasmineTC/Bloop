import numpy as np
import matplotlib.pylab as plt
from math import pi as Pi
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline

renormalizedParams = {'lam1Re': 0.1, 'lam1Im': 0.0, 'lam11': 0.1, 'lam22': 0.1, 'lam12': 0.1, 'lam12p': 0.1, 'yt3': 0.9922814354462509, 
                      'g1': 0.3498206230347181, 'g2': 0.6528876614409878, 'g3': 1.2192627459570353, 
                      'mu3sq': 7812.5, 'lam33': 0.12886749199352251, 'mu12sqRe': 117.5, 'mu12sqIm': 0.0, 'lam2Re': -0.0005734361980235576, 
                      'lam2Im': 0.00099322062987593, 'mu2sq': -4710.528856347395, 'lam23': 0.30007679706315876, 'lam23p': -0.29827978978801506, 
                      'mu1sq': -4710.528856347395, 'lam3Re': -0.0005734361980235576, 'lam3Im': 0.00099322062987593, 'lam31': 0.30007679706315876, 
                      'lam31p': -0.29827978978801506, 'RGScale': 91.1876}
VenusBenchMark = np.zeros(24)
del renormalizedParams['RGScale']
for i, value in enumerate(renormalizedParams.values()):
    VenusBenchMark[i] = value

def beta_functions(coupling, mu):
    #Unpacking the coupling array, to reduce lines could be done with g1, g2,...Mu3sq = coupling[:]
    lam1Re = coupling[0] #initial[key,value]
    lam1Im = coupling[1]
    lam11 = coupling[2]
    lam22 = coupling[3]
    lam12 = coupling[4]
    lam12p = coupling[5]
    yt3 = coupling[6]
    g1 = coupling[7]
    g2 = coupling[8]
    g3 = coupling[9]
    mu3sq = coupling[10]
    lam33 = coupling[11]
    mu12sqRe = coupling[12]
    mu12sqIm = coupling[13]
    lam2Re = coupling[14]
    lam2Im = coupling[15]
    mu2sq = coupling[16]
    lam23 = coupling[17]
    lam23p = coupling[18]
    mu1sq = coupling[19]
    lam3Re = coupling[20]
    lam3Im = coupling[21]
    lam31 = coupling[22]
    lam31p = coupling[23]
       
    #Each differential equation is copy and pasted from the DRalgo file BetaFunctions4D[]//FortranForm 
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
    
    return  [dlam1Re, dlam1Im, dlam11, dlam22, dlam12, dlam12p, dyt3, dg1, dg2, dg3, dmu3sq, dlam33,
     dmu12sqRe, dmu12sqIm, dlam2Re, dlam2Im, dmu2sq, dlam23, dlam23p, dmu1sq, dlam3Re, dlam3Im, dlam31, dlam31p]   


#[dg1,dg2, g3, lambda11, lambda12p, lambda12, lambda1Im, lambda1Re, lambda22, dlambda23p, dlambda23, dlambda2Im, dlambda2Re, dlambda31p, dlambda31, dlambda33, dlambda3Im, dlambda3Re, dyt3, dmu12SqIm, dmu12SqRe, dmu1Sq, dmu2Sq, dmu3Sq]
InitialConditions = [np.sqrt(0.15), np.sqrt(0.4), np.sqrt(1.9), 0.1, 0.1, 0.1, 0.0, 0.1, 0.1, 1, 1.1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.61, 0.71, 1.0, 0.91, 1.01, 0.92, 0.82, 0.75]

muRange = np.linspace(90, 700, 300)
mubarRange = np.log(muRange)

solution = odeint(beta_functions, VenusBenchMark, mubarRange)
solution = np.transpose(solution)
#test1 = CubicSpline(mubarRange, Solution_trans[0], extrapolate=False)

plt.plot(mubarRange, solution[6], label ='yt')
#plt.plot(mubarRange, test1(mubarRange), label ='U(1) spline')
plt.plot(mubarRange, solution[7], label = 'U(1)')
plt.plot(mubarRange, solution[8], label = 'SU(2)')
plt.plot(mubarRange, solution[9], label = 'SU(3)')
plt.legend(loc ='best')
plt.title("Running of the gauge couplings and top yukawa")
plt.xlabel('$log(\mu)$')
plt.ylabel('coupling')
plt.show()

def InterpolateBetaFunction():     
    dlam1Re = CubicSpline(muRange, solution[0], extrapolate=False)
    dlam1Im = CubicSpline(muRange, solution[1], extrapolate=False)
    dlam11 = CubicSpline(muRange, solution[2], extrapolate=False)
    dlam22 = CubicSpline(muRange, solution[3], extrapolate=False)
    dlam12 = CubicSpline(muRange, solution[4], extrapolate=False)
    dlam12p = CubicSpline(muRange, solution[5], extrapolate=False)
    dyt3 = CubicSpline(muRange, solution[6], extrapolate=False)
    dg1 = CubicSpline(muRange, solution[7], extrapolate=False)
    dg2 = CubicSpline(muRange, solution[8], extrapolate=False)
    dg3 = CubicSpline(muRange, solution[9], extrapolate=False)
    dmu3sq = CubicSpline(muRange, solution[10], extrapolate=False)
    dlam33 = CubicSpline(muRange, solution[11], extrapolate=False)
    dmu12sqRe = CubicSpline(muRange, solution[12], extrapolate=False)
    dmu12sqIm = CubicSpline(muRange, solution[13], extrapolate=False)
    dlam2Re = CubicSpline(muRange, solution[14], extrapolate=False)
    dlam2Im = CubicSpline(muRange, solution[15], extrapolate=False)
    dmu2sq = CubicSpline(muRange, solution[16], extrapolate=False)
    dlam23 = CubicSpline(muRange, solution[17], extrapolate=False)
    dlam23p = CubicSpline(muRange, solution[18], extrapolate=False)
    dmu1sq = CubicSpline(muRange, solution[19], extrapolate=False)
    dlam3Re = CubicSpline(muRange, solution[20], extrapolate=False)
    dlam3Im = CubicSpline(muRange, solution[21], extrapolate=False)
    dlam31 = CubicSpline(muRange, solution[22], extrapolate=False)
    dlam31p = CubicSpline(muRange, solution[23], extrapolate=False)
        
    return [dlam1Re, dlam1Im, dlam11, dlam22, dlam12, dlam12p, dyt3, dg1, dg2, dg3, dmu3sq, dlam33,
                dmu12sqRe, dmu12sqIm, dlam2Re, dlam2Im, dmu2sq, dlam23, dlam23p, dmu1sq, dlam3Re, dlam3Im, dlam31, dlam31p]  
array = np.zeros(24)
interp = InterpolateBetaFunction()

# for i, (key, value) in enumerate (renormalizedParams.items()):
#     array[i] = dlam1Re(200)
#     #EvalulateCoupling[key] = BetaFunctionsInterp[i](muEvaulate)
# print (array)
#print (interp[0](150))
print (   type( np.ndarray.item(interp[0](150))  )   )
