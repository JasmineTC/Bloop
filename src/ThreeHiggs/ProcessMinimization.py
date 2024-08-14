import numpy as np

def makeFieldDimensionless(temp: list[float], field: list[float]) -> list[float]:
    return field/np.sqrt(temp)

def jumpFinder(array: np.ndarray[float])-> np.ndarray[int]:
    return np.nonzero(np.abs(array) > 0.2)[0]

def interpretData(result,index: int, bmInput: dict[str, float]):
    interpResult = {"bmNumber": index,
                    "bBoundFromBelow": result["bBoundFromBelow"],
                    "bIsPerturbative": result["bIsPerturbative"],
                    "bad": False, 
                    "jumpsv1": [],
                    "jumpsv2": [],
                    "jumpsv3": [],
                    "UltraSoftTemp": result["UltraSoftTemp"],
                    "bmInput": bmInput}
    if not result["bBoundFromBelow"]:
        return interpResult ##If unbounded from below no point doing these calcs
    
    v1ListRenormDiff = np.diff( makeFieldDimensionless(result["T"], result["minimumLocation"][0]) )
    v2ListRenormDiff = np.diff( makeFieldDimensionless(result["T"], result["minimumLocation"][1]) )
    v3ListRenormDiff = np.diff( makeFieldDimensionless(result["T"], result["minimumLocation"][2]) )
    
    jumpv1 = jumpFinder(v1ListRenormDiff)
    jumpv2 = jumpFinder(v2ListRenormDiff)
    jumpv3 = jumpFinder(v3ListRenormDiff)

    if len(jumpv1) == 0 and len(jumpv2) == 0 and len(jumpv3) == 1: ##Standard one step case
        interpResult["jumpsv3"].append(( v3ListRenormDiff[int(jumpv3)], 
                                         result["T"][int(jumpv3)] ))
    
    elif len(jumpv1) >0 and len(jumpv2) >0 and len(jumpv3) >0: ## 2 step case
    ## If v1 and v2 are large by our lowest temp (i.e. the 0th element T~50GeV) then v3 might not be global min or transition might happen too late
        if result["minimumLocation"][0][0] + result["minimumLocation"][1][0] > 0.1:
            result["bad"] = True
        for val in jumpv1:
            interpResult["jumpsv1"].append(( v1ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
        for val in jumpv2:
            interpResult["jumpsv2"].append(( v2ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
        for val in jumpv3:
            interpResult["jumpsv3"].append(( v3ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
    else: ## Something has probably gone wrong if this gets hit
        result["bad"] = True
        for val in jumpv1:
            interpResult["jumpsv1"].append(( v1ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
        for val in jumpv2:
            interpResult["jumpsv2"].append(( v2ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
        for val in jumpv3:
            interpResult["jumpsv3"].append(( v3ListRenormDiff[int(val)], 
                                             result["T"][int(val)] ))
    return interpResult
    
