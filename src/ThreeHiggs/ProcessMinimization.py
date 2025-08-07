import numpy as np
from math import sqrt

def makeFieldDimensionless(temp: list[float], field: list[float]) -> list[float]:
    return field/np.sqrt(temp)

def jumpFinder(array: np.ndarray[float])-> np.ndarray[int]:
    return np.nonzero(np.abs(array) > 0.3)[0]

def interpretData(result: dict, index: int, bmInput: dict[str, float]):
    interpResult = {"badReason": result["badReason"],
                    "bIsPerturbative": result["bIsPerturbative"],
                    "jumpsv1": [],
                    "jumpsv2": [],
                    "jumpsv3": [],
                    "strong": False,
                    "bmNumber": index,
                    "bmInput": bmInput}

    v1Dimless = makeFieldDimensionless(result["T"], result["minimumLocation"][0])
    v2Dimless = makeFieldDimensionless(result["T"], result["minimumLocation"][1])
    v3Dimless = makeFieldDimensionless(result["T"], result["minimumLocation"][2])
    
    v1DimlessDiff = np.diff( v1Dimless )
    v2DimlessDiff = np.diff( v2Dimless )
    v3DimlessDiff = np.diff( v3Dimless )
    
    jumpv1 = jumpFinder( v1DimlessDiff )
    jumpv2 = jumpFinder( v2DimlessDiff )
    jumpv3 = jumpFinder( v3DimlessDiff )
    
    if len(jumpv1) > 0: 
        for val in jumpv1:
            interpResult["jumpsv1"].append(( v1DimlessDiff[val], 
                                             result["T"][val] ))
    if len(jumpv2) > 0:
        for val in jumpv2:
            interpResult["jumpsv2"].append(( v2DimlessDiff[val], 
                                             result["T"][val] ))
    if len(jumpv3) > 0:
        for val in jumpv3:
            interpResult["jumpsv3"].append(( v3DimlessDiff[val], 
                                             result["T"][val] ))
    strength = 0
    if max(abs(v3DimlessDiff[jumpv3]), default = 0) > 0.3:
        phaseJumpIdx = jumpv3[ np.argmax( abs( v3DimlessDiff[jumpv3] ) ) ]
        strength = sqrt( v1Dimless[phaseJumpIdx]**2 +\
                                       v2Dimless[phaseJumpIdx]**2 + \
                                       v3Dimless[phaseJumpIdx]**2) - \
                                 sqrt( v1Dimless[phaseJumpIdx+1]**2 +\
                                       v2Dimless[phaseJumpIdx+1]**2 + \
                                       v3Dimless[phaseJumpIdx+1]**2)
    
    interpResult["strong"] = strength if strength > 0.3 else False
    interpResult["step"] = 2 if (len(jumpv1)+len(jumpv2)) >0  else 1
    return interpResult
    
