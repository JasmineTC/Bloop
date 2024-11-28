'''Generate benchmarks at 0T satisfying some constraints. 
If the functions return True then the benchmark point is allowed by those constraints
Variables that are bools are denoted with a b prefix'''

import numpy as np
import math as m

def bPositiveMassStates(mu2sq, mu12sq, lam23, lam23p, lambdaMinus, lambdaPlus, vsq) -> bool:
    return -mu2sq - mu12sq + lam23*vsq/2 >0 and \
           -mu2sq + mu12sq + lam23*vsq/2 >0 and \
           -mu2sq + (lam23 + lam23p)*vsq/2 - lambdaMinus >0 and \
           -mu2sq + (lam23 + lam23p)*vsq/2 + lambdaMinus >0 and \
           -mu2sq + (lam23 + lam23p)*vsq/2 - lambdaPlus >0 and \
           -mu2sq + (lam23 + lam23p)*vsq/2 + lambdaPlus>0

def bIsBounded(lam11, lam22, lam33, lamx, lamy, lamz) -> bool:
    return lam11 > 0 and \
           lam22 > 0 and \
           lam33 > 0 and \
           lamx > -2*m.sqrt(lam11*lam22) and \
           lamy > -2*m.sqrt(lam11*lam33) and \
           lamz > -2*m.sqrt(lam22*lam33) and \
           (m.sqrt(lam33)*lamx + m.sqrt(lam11)*lamz + m.sqrt(lam22)*lamy >=0 or \
           lam33*lamx**2 + lam11*lamz**2 + lam22*lamy**2 -lam11*lam22*lam33 - 2*lamx*lamy*lamz < 0)

def bNoLightCharged(mSpm1, mSpm2) -> bool:
    return mSpm1 >= 90 and \
           mSpm2 >= 90

def lagranianParamGen(mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy, bmNumber):
    vsq = 246.22**2
    ## Some 'dark' sector params we keep fixed to keep the scan managable
    lam1Re = 0.1
    lam1Im = 0.0
    lam11 = 0.11
    lam22 = 0.12
    lam12 = 0.13
    lam12p = 0.14
    ## SM params
    mu3sq = 125**2/2 ## Higgs mass*2/2
    lam33 = ((125.00)**2 / (2.*vsq)) ## Higgs mass^2 / 2*vev**2
    ## Helper param
    lamx = lam12 + min(0, lam12p - 2*m.sqrt(lam1Re**2 + lam1Im**2) )
    

    mS2 = delta12 + mS1
    mSpm1 = delta1c + mS1
    mSpm2 = deltac + mSpm1
    mu12sq = (mSpm2**2 - mSpm1**2)/2
    
    sinTheta, cosTheta = m.sin(thetaCPV), m.cos(thetaCPV)
    lam2absInsideSqR = (2.*mu12sq*cosTheta)**2 + (mS2**2 - mS1**2)**2 - (mSpm2**2 - mSpm1**2)**2
    if lam2absInsideSqR < 0:
        return False
    
    lam2Abs = ( mu12sq*cosTheta + m.sqrt( lam2absInsideSqR )/4 ) /vsq
    lam2Re = lam2Abs*cosTheta
    lam2Im = lam2Abs*sinTheta
    
    lambdaMinus = m.sqrt( mu12sq**2 + vsq**2*lam2Abs**2 - 2.*vsq*mu12sq*lam2Abs*cosTheta)
    lambdaPlus = m.sqrt( mu12sq**2 + vsq**2*lam2Abs**2 + 2.*vsq*mu12sq*lam2Abs*cosTheta)
    alpha = (-mu12sq + vsq*lam2Abs*cosTheta - lambdaMinus) / ( (vsq*lam2Abs*sinTheta) +1e-100 )
    mu2sq = vsq/2. * ghDM - vsq / (alpha**2 + 1.) * (2.*alpha* sinTheta + (alpha**2 - 1.)*cosTheta) * lam2Abs - (mS2**2 + mS1**2)/2
    lam23 = (2.*mu2sq + mSpm2**2 + mSpm1**2)/vsq
    lam23p = (mS2**2 + mS1**2 - mSpm2**2 - mSpm1**2)/vsq
    
    mu1sq = darkHieracy*mu2sq
    lam3Re = darkHieracy*lam2Re
    lam3Im = darkHieracy*lam2Im
    lam31 = darkHieracy*lam23
    lam31p = darkHieracy*lam23p
    
    lamy = lam31 + min(0, lam31p - 2*m.sqrt(lam3Re**2 + lam3Im**2) )
    lamz = lam23 + min(0, lam23p - 2*m.sqrt(lam2Re**2 + lam2Im**2) )

    if bIsBounded(lam11, lam22, lam33, lamx, lamy, lamz) \
        and bNoLightCharged(mSpm1, mSpm2) \
        and bPositiveMassStates(mu2sq, mu12sq, lam23, lam23p, lambdaMinus, lambdaPlus, vsq):
        return{
            "bmNumber": bmNumber,
            "RGScale": 91.1876,
            "bPreCalculated": True,
            
            "bmInput": {"thetaCPV" : thetaCPV,
                         "ghDM" : ghDM,
                         "mS1" : mS1,
                         "delta12" : delta12,
                         "delta1c" : delta1c,
                         "deltac" : deltac,
                         "darkHieracy": darkHieracy},
            
            "massTerms": {"mu12sqRe" : mu12sq, 
                          "mu12sqIm" : 0, 
                          "mu2sq": mu2sq, 
                          "mu3sq":mu3sq,
                          "mu1sq": mu1sq,},
                                
            "couplingValues":{"lam1Re" : lam1Re, 
                              "lam1Im" : lam1Im, 
                              "lam2Re" : lam2Abs*cosTheta, 
                              "lam2Im" : lam2Abs*sinTheta, 
                              "lam11" : lam11, 
                              "lam22" : lam22, 
                              "lam12" : lam12, 
                              "lam12p" : lam12p, 
                              "lam23": lam23,
                              "lam23p": lam23p,
                              "lam3Re": lam3Re,
                              "lam3Im": lam3Im,
                              "lam31": lam31,
                              "lam31p": lam31p,
                              "lam33": lam33}
            
        }
    return False

def randomBmParam():
    bmdictList = []
    while len(bmdictList) < 1e6:
        mS1 = np.random.uniform(63, 100)
        delta12 = np.random.uniform(5, 100)
        delta1c = np.random.uniform(5, 100)
        deltac = np.random.uniform(5, 100)
        ghDM = np.random.uniform(0, 1)
        thetaCPV = np.random.uniform(np.pi/2, 3*np.pi/2)
        darkHieracy = 1
        bmDict = lagranianParamGen(mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy, len(bmdictList))
        if bmDict:
            bmdictList.append(bmDict)
    return bmdictList

def notRandombmParam():
    bmdictList = []
    bmInputList = [[300, 0, 0, 0, 0.0, 0.0, 1],
                   [67, 4.0, 50.0, 1.0, 0.0,  2.*np.pi/3, 1],
                   [57, 8.0, 50.0, 1.0, 0.0,  2.*np.pi/3, 1],
                   [70, 12.0, 50.0, 1.0, 0.0,  2.*np.pi/3, 1],
                   [48, 20.0, 50.0, 1.0, 0.0, 5.*np.pi/6., 1],
                   [75, 55.0, 50.0, 1.0, 0.0,  2.*np.pi/3, 1],
                   [74, 55.0, 50.0, 15.0, 0.0,  2.*np.pi/3, 1],
                   [90, 5.0, 1.0, 1.0, 0.0,  2.*np.pi/3, 1],
                   [90, 55.0, 1.0, 1.0, 0.0,  2.*np.pi/3, 1],
                   [90, 55.0, 1.0, 22.0, 0.0,  2.*np.pi/3, 1],
                   [50, 55, 70, 25, 3/9, np.pi/2., 1]]
    for bmInput in bmInputList:
        mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy = bmInput
        bmDict = lagranianParamGen(mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy, len(bmdictList))
        if bmDict:
            bmdictList.append(bmDict)
    return bmdictList


from json import dump                       
with open("Test.json", "w") as outfile: 
    dump(notRandombmParam(), outfile, indent = 4)
