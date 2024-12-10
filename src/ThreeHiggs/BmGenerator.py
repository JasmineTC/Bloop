'''Generate benchmarks at 0T satisfying some constraints. 
If the functions return True then the benchmark point is allowed by those constraints
Variables that are bools are denoted with a b prefix'''

import math as m
def bIsBounded(param : dict[str, float]) -> bool:
    ## Taking equations 26-31 from the draft that ensure the potential is bounded from below.
    if not param["lam11"] > 0:
        return False
    if not param["lam22"] > 0:
        return False
    if not param["lam33"] > 0:
        return False
    
    lamx = param["lam12"] + min(0, param["lam12p"] - 2*m.sqrt(param["lam1Re"]**2 + param["lam1Im"]**2) )
    lamy = param["lam31"] + min(0, param["lam31p"] - 2*m.sqrt(param["lam3Re"]**2 + param["lam3Im"]**2) )
    lamz = param["lam23"] + min(0, param["lam23p"] - 2*m.sqrt(param["lam2Re"]**2 + param["lam2Im"]**2) )
    if not lamx > -2*m.sqrt(param["lam11"]*param["lam22"]):
        return False
    if not lamy > -2*m.sqrt(param["lam11"]*param["lam33"]):
        return False
    if not lamz > -2*m.sqrt(param["lam22"]*param["lam33"]):
        return False
    if not (m.sqrt(param["lam33"])*lamx + m.sqrt(param["lam11"])*lamz + m.sqrt(param["lam22"])*lamy >= 0 or \
            param["lam33"]*lamx**2 + param["lam11"]*lamz**2 + param["lam22"]*lamy**2 -param["lam11"]*param["lam22"]*param["lam33"] - 2*lamx*lamy*lamz < 0):
        return False
    return True
    
def _bPositiveMassStates(mu2sq, mu12sq, lam23, lam23p, lambdaMinus, lambdaPlus, vsq) -> bool:
    return -mu2sq - mu12sq + lam23*vsq/2 >0 and \
           -mu2sq + mu12sq + lam23*vsq/2 >0 and \
           -mu2sq + (lam23 + lam23p)*vsq/2 - lambdaMinus >0 and \
           -mu2sq + (lam23 + lam23p)*vsq/2 + lambdaMinus >0 and \
           -mu2sq + (lam23 + lam23p)*vsq/2 - lambdaPlus >0 and \
           -mu2sq + (lam23 + lam23p)*vsq/2 + lambdaPlus>0

def _bNoLightCharged(mSpm1, mSpm2) -> bool:
    return mSpm1 >= 90 and \
           mSpm2 >= 90
           
from math import sin, cos
def _lagranianParamGen(mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy, bmNumber):
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

    mS2 = delta12 + mS1
    mSpm1 = delta1c + mS1
    mSpm2 = deltac + mSpm1
    mu12sq = (mSpm2**2 - mSpm1**2)/2
    
    sinTheta, cosTheta = sin(thetaCPV), cos(thetaCPV)
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
    
    if not _bNoLightCharged(mSpm1, mSpm2):
        return False
    if not _bPositiveMassStates(mu2sq, mu12sq, lam23, lam23p, lambdaMinus, lambdaPlus, vsq):
        return False
    paramDict = {"bmNumber": bmNumber,
            "RGScale": 91.1876,
            
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
                              "lam33": lam33}}
    if not bIsBounded(paramDict["massTerms"] | paramDict["couplingValues"]):
        return False
    return paramDict


import numpy as np
def _randomBmParam(num):
    bmdictList = []
    while len(bmdictList) < num:
        mS1 = np.random.uniform(63, 100)
        delta12 = np.random.uniform(5, 100)
        delta1c = np.random.uniform(5, 100)
        deltac = np.random.uniform(5, 100)
        ghDM = np.random.uniform(0, 1)
        thetaCPV = np.random.uniform(np.pi/2, 3*np.pi/2)
        darkHieracy = 1
        bmDict = _lagranianParamGen(mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy, len(bmdictList))
        if bmDict:
            bmdictList.append(bmDict)
    return bmdictList

def _notRandomBmParam():
    bmdictList = []
    bmInputList = [[300, 0, 0, 0, 0.0, 0.0, 1],
                   [67, 4.0, 50.0, 1.0, 0.0, 2.*np.pi/3, 1],
                   [57, 8.0, 50.0, 1.0, 0.0, 2.*np.pi/3, 1],
                   [70, 12.0, 50.0, 1.0, 0.0, 2.*np.pi/3, 1],
                   [48, 20.0, 50.0, 1.0, 0.0, 5.*np.pi/6., 1],
                   [75, 55.0, 50.0, 1.0, 0.0, 2.*np.pi/3, 1],
                   [74, 55.0, 50.0, 15.0, 0.0, 2.*np.pi/3, 1],
                   [90, 5.0, 1.0, 1.0, 0.0, 2.*np.pi/3, 1],
                   [90, 55.0, 1.0, 1.0, 0.0, 2.*np.pi/3, 1],
                   [90, 55.0, 1.0, 22.0, 0.0, 2.*np.pi/3, 1],
                   [50, 55, 70, 25, 3/9, np.pi/2., 1]]
    for bmInput in bmInputList:
        mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy = bmInput
        bmDict = _lagranianParamGen(mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy, len(bmdictList))
        if bmDict:
            bmdictList.append(bmDict)
    return bmdictList

def _strongSubSet():
    from json import load
    dictList = load(open("Benchmarks/StrongBmList.json", "r"))
    bmDict = []
    for ele in dictList:
        bmDict.append(_lagranianParamGen(ele["mS1"],
                    ele["delta12"],
                    ele["delta1c"], 
                    ele["deltac"], 
                    ele["ghDM"], 
                    ele["thetaCPV"],
                    ele["darkHieracy"],
                    ele["bmNumber"]))
    return bmDict
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", action = "store",  dest = "mode", default = "handPicked",
                    choices = ["handPicked", "random", "randomSSS"],
                    help = "Str: Specify the mode to generate bm with.")
    parser.add_argument("--randNum",type = int, action = "store",  dest = "randNum", default = 1_000_000,
                    help = "Int: Specify how many random bm to generate.")
    

    args = parser.parse_args()
    from json import dump                       
    if args.mode == "handPicked":
        dump(_notRandomBmParam(), open("Benchmarks/handPicked.json", "w"), indent = 4)
    elif args.mode == "random":
        print("no")
        #dump(_randomBmParam(args.randNum), open("Benchmarks/randomScan.json", "w"), indent = 4)
    else:
        dump(_strongSubSet(), open("Benchmarks/randomScanSSS.json", "w"), indent = 4)



    
from unittest import TestCase
class BmGeneratorUnitTests(TestCase):
    def test_bIsBoundedFalse(self):
        source = {'mu12sqRe': 12724.595168103033, 'mu12sqIm': 0, 'mu2sq': -20604.175986862854, 'mu3sq': 7812.5, 'mu1sq': -20604.175986862854, 'lam1Re': 0.1, 'lam1Im': 0.0, 'lam2Re': 0.08368703688875163, 'lam2Im': 0.05054893896685599, 'lam11': 0.11, 'lam22': 0.12, 'lam12': 0.13, 'lam12p': 0.14, 'lam23': 0.13184965469208287, 'lam23p': -0.10179114069822925, 'lam3Re': 0.08368703688875163, 'lam3Im': 0.05054893896685599, 'lam31': 0.13184965469208287, 'lam31p': -0.10179114069822925, 'lam33': 0.12886749199352251}
        self.assertEqual(False, bIsBounded(source))
        
    def test_bIsBoundedTrue(self):
        source = {'mu12sqRe': 4799.941333804141, 'mu12sqIm': 0, 'mu2sq': -11505.594825493996, 'mu3sq': 7812.5, 'mu1sq': -11505.594825493996, 'lam1Re': 0.1, 'lam1Im': 0.0, 'lam2Re': 0.010823997158779158, 'lam2Im': 0.017584425457946057, 'lam11': 0.11, 'lam22': 0.12, 'lam12': 0.13, 'lam12p': 0.14, 'lam23': 0.33218995023667736, 'lam23p': -0.29472720912675215, 'lam3Re': 0.010823997158779158, 'lam3Im': 0.017584425457946057, 'lam31': 0.33218995023667736, 'lam31p': -0.29472720912675215, 'lam33': 0.12886749199352251}
        self.assertEqual(True, bIsBounded(source))
    
    def test___bPositiveMassStatesFalse(self):
        source = (-7446.958243389069, 14066.835009399923, 1.0335363761664145, -0.8185339583108823, 6510.844284793283, 23075.04842411651, 60624.2884)
        self.assertEqual(False, _bPositiveMassStates(*source))
        
    def test___bPositiveMassStatesTrue(self):
        source = (5142.163347176796, 2301.350882335473, 0.6120747865923423, -0.12441080416604636, 2931.315410805045, 2448.0647739410247, 60624.2884)
        self.assertEqual(True, _bPositiveMassStates(*source))
        
    def test_lagranianParamGen(self):
        reference = {'bmNumber': 0, 'RGScale': 91.1876, 'bmInput': {'thetaCPV': 3.11308902835221, 'ghDM': 0.15520161865427817, 'mS1': 89.15641588128479, 'delta12': 87.17952518246265, 'delta1c': 14.020273320699415, 'deltac': 5.129099092707543, 'darkHieracy': 1}, 'massTerms': {'mu12sqRe': 542.3572917258725, 'mu12sqIm': 0, 'mu2sq': -9572.907910345595, 'mu3sq': 7812.5, 'mu1sq': -9572.907910345595}, 'couplingValues': {'lam1Re': 0.1, 'lam1Im': 0.0, 'lam2Re': -0.08646867756101999, 'lam2Im': 0.0024653384763691356, 'lam11': 0.11, 'lam22': 0.12, 'lam12': 0.13, 'lam12p': 0.14, 'lam23': 0.05327497010465927, 'lam23p': 0.274933662245124, 'lam3Re': -0.08646867756101999, 'lam3Im': 0.0024653384763691356, 'lam31': 0.05327497010465927, 'lam31p': 0.274933662245124, 'lam33': 0.12886749199352251}}
        source = (89.15641588128479, 87.17952518246265, 14.020273320699415, 5.129099092707543, 0.15520161865427817, 3.11308902835221, 1, 0)
        self.assertEqual( reference, _lagranianParamGen(*source) )
