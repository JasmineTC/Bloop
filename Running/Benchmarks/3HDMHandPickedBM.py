import numpy as np
__MZ = 91.1876

##These are benchmarks taken from the draft in table 1
bmList = []
## Tries to be close to the SM limit, full decoupling of extra Higgses is not possible though because of the gauge sector.
bm0 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,
    "bmNumber": 0,
    ## "Physical" masses of the scalar sector at tree level
    "mS1" : 300,
    "bMassSplittingInput": True, ## Use mass splittings as Venus does or directly put in masses
    "preCalculated": False, ## Are all non SM couplings calculated before hand or does generic model need to do it?
    "delta12" : 0, # mS2 - mS1
    "delta1c" : 0, # mSpm1 - mS1 => set 0 to make lam23p = 0
    "deltac" : 0, # mSpm2 - mSpm1, splitting between charged scalars => if degenerate, mu12^2 = 0 

    "thetaCPV" : 0.0, ##CP violating angle
    "ghDM" : 0.0, ## Higgs coupling to the DM candidate i.e. H S1 S1
    "darkHierarchy": 1, ## Multiplicative value between dark sector params
    ## Dark sector parameters
    "lam1Re" : 1e-12,
    "lam1Im" : 1e-12,
    "lam11" : 1e-12,
    "lam22" : 1e-12,
    "lam12" : 1e-12,
    "lam12p" : 1e-12
}
bmList.append(bm0)

bm1 = {
    "bmNumber": 1,
    "RGScale" : __MZ,
    "mS1" : 67,
    "bMassSplittingInput": True,
    "bPreCalculated": False,
    "delta12" : 4.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,
    "darkHierarchy": 1,
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm1)

bm2 = {
    "RGScale" : __MZ,
    "bmNumber": 2,
    "mS1" : 57,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 8.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm2)

bm3 = {
    "RGScale" : __MZ,
    "bmNumber": 3,
    "mS1" : 70,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 12.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm3)

bm4 = {
    "RGScale" : __MZ,
    "bmNumber": 4,
    "mS1" : 48,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 20.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 5.*np.pi/6.,
    "ghDM" : 0.0,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm4)

bm5 = {
    "RGScale" : __MZ,
    "bmNumber": 5,
    "mS1" : 75,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm5)

bm6 = {
    "RGScale" : __MZ,
    "bmNumber": 6,
    "mS1" : 74,
    "bMassSplittingInput": True, 
    "preCalculated": False,
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 15.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm6)

bm7 = {
    "RGScale" : __MZ,
    "bmNumber": 7,
    "mS1" : 90,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 5.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    "lam1Re"    : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm7)

bm8 = {
    "RGScale" : __MZ,
    "bmNumber": 8,
    "mS1" : 90,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm8)

bm9 = {
    "RGScale" : __MZ,
    "bmNumber": 9,
    "mS1" : 90,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 22.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm9)

##Fake benchmark that is quick to run to keep naming convention the same as draft
bm10 = {
    "RGScale" : __MZ,
    "bmNumber": 10,
    "mS1" : 300,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 0,
    "delta1c" : 0,
    "deltac" : 0, 

    "thetaCPV" : 0.0,
    "ghDM" : 0.0,

    "lam1Re" : 1e-12,
    "lam1Im" : 1e-12,
    "lam11" : 1e-12,
    "lam22" : 1e-12,
    "lam12" : 1e-12,
    "lam12p" : 1e-12
}
bmList.append(bm10)
##Same as benchmark as 10 before it but with higgs dark matter coupling turned on
bm11 = {
    "RGScale" : __MZ,
    "bmNumber": 11,
    "mS1" : 67,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 4.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm11)

bm12 = {
    "RGScale" : __MZ,
    "bmNumber": 12,
    "mS1" : 57,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 8.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm12)

bm13 = {
    "RGScale" : __MZ,
    "bmNumber": 13,
    "mS1" : 70,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 12.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm13)

bm14 = {
    "RGScale" : __MZ,
    "bmNumber": 14,
    "mS1" : 48,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 20.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 5.*np.pi/6.,
    "ghDM" : 0.01,
    
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm14)

bm15 = {
    "RGScale" : __MZ,
    "bmNumber": 15,
    "mS1" : 75,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm15)

bm16 = {
    "RGScale" : __MZ,
    "bmNumber": 16,
    "mS1" : 74,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 15.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm16)

bm17 = {
    "RGScale" : __MZ,
    "bmNumber": 18,
    "mS1" : 90,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 5.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm17)

bm18 = {
    "RGScale" : __MZ,
    "bmNumber": 18,
    "mS1" : 90,
    "bMassSplittingInput": True, 
    "preCalculated": False,
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm18)

bm19 = {
    "RGScale" : __MZ,
    "bmNumber": 19,
    "mS1" : 90,
    "bMassSplittingInput": True,
    "preCalculated": False,
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 22.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.11,
    "lam22" : 0.12,
    "lam12" : 0.13,
    "lam12p" : 0.14
}
bmList.append(bm19)

import json
with open('Running/Benchmarks/Benchmarks_3HDM.json', 'w') as file:
    json.dump(bmList, file, indent=4)