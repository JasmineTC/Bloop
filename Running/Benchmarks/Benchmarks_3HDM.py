import numpy as np

## This just defines some fixed benchmark points to use when testing

__MZ = 91.1876

##These are benchmarks taken from the draft in table 1
bmList = []

## Tries to be close to the SM limit, full decoupling of extra Higgses is not possible though because of the gauge sector.
bm0 = {

    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 300,

    "delta12" : 0, # mS2 - mS1
    "delta1c" : 0, # mSpm1 - mS1 => set 0 to make lam23p = 0
    "deltac" : 0, # mSpm2 - mSpm1, splitting between charged scalars => if degenerate, mu12^2 = 0 

    "thetaCPV" : 0.0,
    "ghDM" : 0.0,

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
    
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 67,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 4.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm1)

bm2 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 57,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 8.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm2)

bm3 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 70,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 12.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm3)

bm4 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 48,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 20.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 5.*np.pi/6.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm4)

bm5 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 75,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm5)

bm6 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 74,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 15.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm6)

bm7 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 5.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm7)

bm8 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm8)

bm9 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 22.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm9)

##Fake benchmark that is quick to run to keep naming convention the same as draft
bm10 = {

    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 300,

    "delta12" : 0, # mS2 - mS1
    "delta1c" : 0, # mSpm1 - mS1 => set 0 to make lam23p = 0
    "deltac" : 0, # mSpm2 - mSpm1, splitting between charged scalars => if degenerate, mu12^2 = 0 

    "thetaCPV" : 0.0,
    "ghDM" : 0.0,

    ## Dark sector parameters
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
    
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 67,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 4.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm11)

bm12 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 57,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 8.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm12)

bm13 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 70,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 12.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm13)

bm14 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 48,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 20.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 5.*np.pi/6.,
    "ghDM" : 0.01,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm14)

bm15 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 75,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm15)

bm16 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 74,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 15.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm16)

bm17 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 5.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm17)

bm18 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm18)

bm19 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 22.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0.02,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm19)
##TODO BUG FIXING NOT REAL BENCHMARKS
bm20 = {
    
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 67,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 4.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : 0,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm20)
bugfix = 1
bm21 = {
    
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 67,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 4.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm21)

bm22 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 57,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 8.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm22)

bm23 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 70,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 12.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm23)

bm24 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 48,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 20.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 5.*np.pi/6.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm24)

bm25 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 75,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm25)

bm26 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 74,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 50.,
    "deltac" : 15.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm26)

bm27 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 5.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm27)

bm28 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 1.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm28)

bm29 = {
    ## Renormalization scale is also an input since we're inputting some action parameters directly
    "RGScale" : __MZ,

    ## "Physical" input in the scalar sector
    "mS1" : 90,
    ## Mass splittings. If you want to input masses directly instead of these deltas, set model.bMassSplittingInput = False
    "delta12" : 55.,
    "delta1c" : 1.,
    "deltac" : 22.,

    "thetaCPV" : 2.*np.pi/3.,
    "ghDM" : bugfix,

    ## Input these dark sector parameters directly. Set to 0.1 because otherwise things become untractable
    "lam1Re" : 0.1,
    "lam1Im" : 0.0,
    "lam11" : 0.1,
    "lam22" : 0.1,
    "lam12" : 0.1,
    "lam12p" : 0.1
}
bmList.append(bm29)