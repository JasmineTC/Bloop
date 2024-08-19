class GenericModel():
    def __init__(self, effectivePotential, dimensionalReduction):
        self.inputParams = {}
        self.effectivePotential = effectivePotential
        self.dimensionalReduction = dimensionalReduction

    def setInputParams(self, inputParams):
        self.inputParams = inputParams

