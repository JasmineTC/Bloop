import numpy as np

def makeFieldDimensionless(temp: list[float], field: list[float]) -> list[float]:
    return field/np.sqrt(temp)

def jumpFinder(array: np.ndarray[float])-> np.ndarray[int]:
    return np.nonzero(np.abs(array) > 0.2)[0]

def interpretData(result: dict, index: int, bmInput: dict[str, float]):
    interpResult = {"bmNumber": index,
                    "failureReason": result["failureReason"],
                    "bIsPerturbative": result["bIsPerturbative"], 
                    "jumpsv1": [],
                    "jumpsv2": [],
                    "jumpsv3": [],
                    "UltraSoftTemp": result["UltraSoftTemp"],
                    "bmInput": bmInput}

    if result["failureReason"]:
        return interpResult 
    
    v1ListRenormDiff = np.diff( makeFieldDimensionless(result["T"], result["minimumLocation"][0]) )
    v2ListRenormDiff = np.diff( makeFieldDimensionless(result["T"], result["minimumLocation"][1]) )
    v3ListRenormDiff = np.diff( makeFieldDimensionless(result["T"], result["minimumLocation"][2]) )
    
    jumpv1 = jumpFinder(v1ListRenormDiff)
    jumpv2 = jumpFinder(v2ListRenormDiff)
    jumpv3 = jumpFinder(v3ListRenormDiff)
    
    if len(jumpv1) > 0:
        for val in jumpv1:
            interpResult["jumpsv1"].append(( v1ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
    
    if len(jumpv2) > 0:
        for val in jumpv2:
            interpResult["jumpsv2"].append(( v2ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
    if len(jumpv3) > 0:
        for val in jumpv3:
            interpResult["jumpsv3"].append(( v3ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
    return interpResult
    
