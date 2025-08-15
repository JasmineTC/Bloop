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
    processedResult = {"bmNumber": bmNumber,
                    "bmInput": bmInput,
                    "strong": False,
                    "complex": False,
                    "results":{}}
    
    if result["failureReason"]:
        return processedResult | {"failureReason": result["failureReason"]}
    
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
                
                if not processedResult["strong"]:
                    processedResult["strong"] = strength if strength > 0.6 else False
                if  processedResult["strong"]:
                    processedResult["strong"] = strength if strength > processedResult["strong"] else processedResult["strong"]
                    
            processedResult["results"][f"{fieldNames[idx]}"] = strengthResults
            
    processedResult["steps"] = len(PTTemps)
    processedResult["bIsPerturbative"] = bool(np.all(result["bIsPerturbative"]))
    imag2RealRatio = abs(np.array(result["vevDepthImag"])/np.array(result["vevDepthReal"]))
    processedResult["complex"] = bool(np.any(imag2RealRatio > 1e-8))
    return processedResult 

