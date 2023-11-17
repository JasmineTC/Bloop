from ParameterMatching import ParameterMatching


softScaleMatching = ParameterMatching()

inputParams = {
    'T' : 100,
    'Lb' : 1,
    'Lf' : 2,
    'g1' : 0.3,
    'g2' : 0.6,
    'g3' : 1.5
}


matchedParams = softScaleMatching.getMatchedParams(inputParams)

print(matchedParams)



