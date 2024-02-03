import numpy as np

## This just defines some fixed benchmark points to use when testing

__MZ = 91.1876

##These are benchmarks taken from the draft in table 1
BM1 = {
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

BM2 = {
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

BM3 = {
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

BM4 = {
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

BM5 = {
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

BM6 = {
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

BM7 = {
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

BM8 = {
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

BM9 = {
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


## Tries to be close to the SM limit, full decoupling of extra Higgses is not possible though because of the gauge sector.
BM_SM_like = {

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

