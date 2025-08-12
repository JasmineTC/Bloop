import numpy as np

def PTStrength(
    idx, 
    fields
): 
    store = np.zeros(2)
    for i in range(2):
        for field in fields:
            store[i] += field[idx+i]**2
    return np.sqrt(store[0] - store[1])

def interpretData(
    result, 
    bmNumber, 
    bmInput,
    fieldNames
):
    interpResult = {"bmNumber": bmNumber,
                    "bmInput": bmInput,
                    "strong": False,
                    "results":{}}
    PTTemps = set()
    allFieldValues = result["vevLocation"]/np.sqrt(result["T"])
    for idx, fieldValues in enumerate(allFieldValues):
        ## Find the indices where a field (dimentionless) changes by more than 0.3
        PTindices = np.nonzero(np.abs(np.diff( fieldValues)) > 0.3)[0]

        if len(PTindices) > 0:
            
            strengthResults = []
            for PTindex in PTindices:
                
                strength = float(PTStrength(PTindex, allFieldValues))                
                T = float(result["T"][PTindex])
                
                strengthResults.append([strength, T])
                PTTemps.add(T)
                
                if not interpResult["strong"]:
                    interpResult["strong"] = strength if strength > 0.6 else False
                if  interpResult["strong"]:
                    interpResult["strong"] = strength if strength > interpResult["strong"] else interpResult["strong"]
                    
            interpResult["results"][f"{fieldNames[idx]}"] = strengthResults
            
    interpResult["steps"] = len(PTTemps)
    interpResult["bIsPerturbative"] = bool(np.all(result["bIsPerturbative"]))
    return interpResult 

    
