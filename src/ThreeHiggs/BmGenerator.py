'''Generate benchmarks at 0T satisfying some constraints. 
If the functions return True then the benchmark point is allowed by those constraints
Variables that are bools are denoted with a b prefix'''

import math as m
from json import load, dump
from ThreeHiggs.ParsedExpression import ParsedExpression
from ThreeHiggs.EffectivePotential import cNlopt
from pathlib import Path
import numpy as np
from os.path import join
from glob import glob

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

def bPhysicalMinimum(nloptInst,
                     potentialFunc,
                     params):
    ## Move these to user args or somewhere out of the way
    minimumInitialGuesses = [[0,0,0],
                             [0,0,246],
                             [100,100,100],
                             [-100,100,100],
                             [50,50,50],
                             [299,299,299],
                             [-299,299,299]]
    
    ## I don't like having the function wrapped here, 
    ## should be defined once at a higher level, not every func call but params
    ## makes that non-trivial 
    potentialFuncWrapped = lambda fields, _: potentialFunc(fields, params)
    
    minValue = nloptInst.nloptGlobal(potentialFuncWrapped, minimumInitialGuesses[0])[1]
    
    for guess in minimumInitialGuesses:
        minLocationTemp, minValueTemp = nloptInst.nloptLocal(potentialFuncWrapped, guess)
        
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
           
def _lagranianParamGen(mS1,
                       delta12, 
                       delta1c, 
                       deltac, 
                       ghDM, 
                       thetaCPV, 
                       darkHieracy, 
                       bmNumber, 
                       nloptInst, 
                       potentialFunc):
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
    if not bIsBounded(params) or not bPhysicalMinimum(nloptInst, potentialFunc, params):
        return False
    return paramDict

def _randomBmParam(randomNum, nloptInst, potentialFunc):
    bmdictList = []
    ## TODO put in some upper limit for this while loop
    darkHieracy = 1
    while len(bmdictList) < randomNum:
        mS1 = np.random.uniform(63, 100)
        delta12 = np.random.uniform(5, 100)
        delta1c = np.random.uniform(5, 100)
        deltac = np.random.uniform(5, 100)
        ghDM = np.random.uniform(0, 1)
        thetaCPV = np.random.uniform(np.pi/2, 3*np.pi/2)
        bmDict = _lagranianParamGen(mS1, delta12, delta1c, deltac, ghDM, thetaCPV, darkHieracy, len(bmdictList), nloptInst, potentialFunc)
        
        if bmDict:
            bmdictList.append(bmDict)
            
    return bmdictList

def _handPickedBm(nloptInst, potentialFunc):
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
        bmDict = _lagranianParamGen(*bmInput, len(bmdictList), nloptInst, potentialFunc)
        if bmDict:
            bmdictList.append(bmDict)
    return bmdictList

def _strongSubSet(prevResultDir, nloptInst, potentialFunc):
    bmdictList = []
    
    for fileName in glob(join(prevResultDir, '*.json')):
        resultDic = load(open(fileName, "r"))
        if resultDic["strong" ] > 0.6:
            bmdictList.append(_lagranianParamGen(resultDic["bmInput"]["mS1"],
                                             resultDic["bmInput"]["delta12"],
                                             resultDic["bmInput"]["delta1c"], 
                                             resultDic["bmInput"]["deltac"], 
                                             resultDic["bmInput"]["ghDM"], 
                                             resultDic["bmInput"]["thetaCPV"],
                                             resultDic["bmInput"]["darkHieracy"],
                                             resultDic["bmNumber"],
                                             nloptInst,
                                             potentialFunc))
    return bmdictList

def generateBenchmarks(args)-> None:
    
    (output_file := Path(args.benchmarkFile)).parent.mkdir(exist_ok=True, 
                                                             parents=True)   
    
    ## Need to take these from args
    nloptInst = cNlopt(config = {"nbrVars": 3, 
                                 "absGlobalTol" : 0.5,
                                 "relGlobalTol" :0.5, 
                                 "absLocalTol" : 0.5, 
                                 "relLocalTol" : 0.5,
                                 "varLowerBounds" : [-300, 0, 0],
                                 "varUpperBounds" : [300, 300, 300]})
    
    parsedExpressions = load(open("../parsedExpressions.json", "r"))
    ## Take the pythonised tree level potential we've generated 
    potential = ParsedExpression(parsedExpressions["veff"]["expressions"][0], None)
    
    ## Feels weird to have nested function but not sure how else to go until can
    ## use arrays with nlopt
    def potentialFunc(fields, params):
        params["v1"] = fields[0]
        params["v2"] = fields[1]
        params["v3"] = fields[2]
        return potential.evaluate(params)
    
    
    ## subset has already been checked to be physical so
    ## shouldn't need nloptInst and potential
    if args.benchmarkType == "randomSSS":
        dump(_strongSubSet(args.prevResultDir, nloptInst, potentialFunc), open(output_file, "w"), indent = 4)
        return

    elif args.benchmarkType == "handPicked":
        dump(_handPickedBm(nloptInst, potentialFunc), open(output_file, "w"), indent = 4)
        return
    
    elif args.benchmarkType == "random":
        dump(_randomBmParam(args.randomNum, nloptInst, potentialFunc), open(output_file, "w"), indent = 4)
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
        
        nloptInst = cNlopt(config = {"nbrVars": 3, 
                                     "absGlobalTol" : 0.5,
                                     "relGlobalTol" :0.5, 
                                     "absLocalTol" : 0.5, 
                                     "relLocalTol" : 0.5,
                                     "varLowerBounds" : [-300, 0, 0],
                                     "varUpperBounds" : [300, 300, 300]})
        ## Gives resource warning for unknown reason
        parsedExpressions = load(open("parsedExpressions.json", "r"))
        
        potential = ParsedExpression(parsedExpressions["veff"]["expressions"][0], None)
        
        ## Feels weird to have nested function but not sure how else to go until can
        ## use arrays with nlopt
        def potentialFunc(fields, params):
            params["v1"] = fields[0]
            params["v2"] = fields[1]
            params["v3"] = fields[2]
            return potential.evaluate(params)
        
        self.assertEqual( reference, _lagranianParamGen(*source, nloptInst, potentialFunc, ) )
