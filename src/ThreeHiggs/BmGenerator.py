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


def potential(params, field):
    lam11 = params["lam11"]
    lam22 = params["lam22"]
    lam33 = params["lam33"]
    lam12 = params["lam12"]
    lam12p = params["lam12p"]
    lam23 = params["lam23"]
    lam23p = params["lam23p"]
    lam31 = params["lam31"]
    lam31p = params["lam31p"]
    lam1Re = params["lam1Re"]
    lam2Re = params["lam2Re"]
    lam3Re = params["lam3Re"]
    mu1sq = params["mu1sq"]
    mu12sqRe = params["mu12sqRe"]
    mu2sq = params["mu2sq"]
    mu3sq = params["mu3sq"]

    return float((field[0]**4*lam11 + field[1]**4*lam22 + field[2]**4*lam33 - 4*field[0]*field[1]*mu12sqRe + field[0]**2*(field[1]**2*(lam12 + lam12p + 2*lam1Re) + field[2]**2*(lam31 + lam31p + 2*lam3Re) - 2*mu1sq) + field[1]**2*(field[2]**2*(lam23 + lam23p + 2*lam2Re) - 2*mu2sq) - 2*field[2]**2*mu3sq)/4)
    #return float((field**4*lam11 + field**4*lam22 + field**4*lam33 - 4*field*field*mu12sqRe + field**2*(field**2*(lam12 + lam12p + 2*lam1Re) + field**2*(lam31 + lam31p + 2*lam3Re) - 2*mu1sq) + field**2*(field**2*(lam23 + lam23p + 2*lam2Re) - 2*mu2sq) - 2*field**2*mu3sq)/4)

import nlopt
from functools import partial
## TODO import the nlopt stuff from EP and use the eval to compute potential at LO
## Overkill for threeHiggs since venus calculated the equations for tree level min
## But this is model indepdent
def bPhysicalMinimum(params):
    minimumInitialGuesses = [[0,0,0],
                             [0,0,246],
                             [100,100,100],
                             [-100,100,100],
                             [50,50,50],
                             [299,299,299],
                             [-299,299,299]]
                             
    potentialPartial = partial(potential, params)
    minFunc = lambda x, _: potentialPartial(x)
    opt = nlopt.opt(nlopt.GN_DIRECT_NOSCAL, 3)
    opt.set_min_objective(minFunc)
    opt.set_lower_bounds([-300,0,0])
    opt.set_upper_bounds([300,300,300])
    opt.set_xtol_abs(0.5)
    opt.set_xtol_rel(0.5)
    minLocation, minValue = opt.optimize([0,0,0]), opt.last_optimum_value()
    opt = nlopt.opt(nlopt.LN_BOBYQA, 3)
    opt.set_min_objective(minFunc)
    opt.set_lower_bounds([-300,0,0])
    opt.set_upper_bounds([300,300,300])
    opt.set_xtol_abs(0.5)
    opt.set_xtol_rel(0.5)
    for guess in minimumInitialGuesses:
        minLocationTemp, minValueTemp = opt.optimize(guess), opt.last_optimum_value()
        if minValueTemp < minValue:
            minLocation, minValue =  minLocationTemp, minValueTemp
    return np.all(np.isclose(minLocation, [0,0, 246], atol=1))

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
    params = paramDict["massTerms"] | paramDict["couplingValues"]
    if not bIsBounded(params):
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

def _handPickedBm():
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

def _strongSubSet(prevResultDir):
    from json import load
    from os.path import join
    from glob import glob
    bmDict = []
    for fileName in glob(join(prevResultDir, '*.json')):
        resultDic = load(open(fileName, "r"))
        if resultDic["strong" ] > 0.6:
            bmDict.append(_lagranianParamGen(resultDic["bmInput"]["mS1"],
                                             resultDic["bmInput"]["delta12"],
                                             resultDic["bmInput"]["delta1c"], 
                                             resultDic["bmInput"]["deltac"], 
                                             resultDic["bmInput"]["ghDM"], 
                                             resultDic["bmInput"]["thetaCPV"],
                                             resultDic["bmInput"]["darkHieracy"],
                                             resultDic["bmNumber"]))
    return bmDict

def generateBenchmarks(benchmarkOutput: str, mode: str, randomNum: int = None, prevResultDir: str = None)-> None:
    from pathlib import Path
    (output_file := Path(benchmarkOutput)).parent.mkdir(exist_ok=True, parents=True)   
    
    from json import dump  
    
    if mode == "handPicked":
        dump(_handPickedBm(), open(output_file, "w"), indent = 4)
        return
    elif mode == "random":
        dump(_randomBmParam(randomNum), open(output_file, "w"), indent = 4)
        return
    elif mode == "randomSSS":
        ## THIS NEEDS UPDATING BUT NOT IN THIS COMMIT
        dump(_strongSubSet(prevResultDir), open(output_file, "w"), indent = 4)
        return 
    return


    
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
