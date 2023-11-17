import numpy as np
import numpy.typing as npt

class EffectivePotential:

    def __init__(self, modelParameters: dict):
        
        self.setModelParameters(modelParameters)


    def setModelParameters(self, modelParameters: dict) -> None:
        self.mu1sq = modelParameters["mu1sq"]
        self.mu2sq = modelParameters["mu2sq"]
        self.mu3sq = modelParameters["mu3sq"]

        self.lam11 = modelParameters["lam11"]
        self.lam22 = modelParameters["lam22"]
        self.lam33 = modelParameters["lam33"]
        self.lam12 = modelParameters["lam12"]
        self.lam23 = modelParameters["lam23"]
        self.lam31 = modelParameters["lam31"]
        self.lam12Prime = modelParameters["lam12Prime"]
        self.lam23Prime = modelParameters["lam23Prime"]
        self.lam31Prime = modelParameters["lam31Prime"]

        self.mu12sqRe = modelParameters["mu12sqRe"]
        self.mu12sqIm = modelParameters["mu12sqIm"]
        self.lam1Re = modelParameters["lam1Re"]
        self.lam1Im = modelParameters["lam1Im"]
        self.lam2Re = modelParameters["lam2Re"]
        self.lam2Im = modelParameters["lam2Im"]
        self.lam3Re = modelParameters["lam3Re"]
        self.lam3Im = modelParameters["lam3Im"]

    ## Evaluate Veff(fields, T) using current set of model parameters
    def evaluate(self, temperature: npt.ArrayLike, fields: npt.ArrayLike) -> complex:
        
        VTotal = self.V0(fields)
        return VTotal

    # just calls self.evaluate
    def __call__(self, temperature: npt.ArrayLike, fields: npt.ArrayLike) -> complex:

        self.evaluate(temperature, fields)


    ## Tree level potential
    def V0(self, fields: npt.ArrayLike) -> float:
        v3, _ = fields 

        # this assumes just one background field (in the phi3 doublet)

        res = -0.5*self.mu3sq * v3**2 + 0.25*self.lam33*v3**4
        return res