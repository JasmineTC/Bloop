'''---prints the number of SFOPT in a dir ---
Loads all the file names in the results directory 
Strips the Results/ part of the file name
Loads the benchmark information from the file
If the benchmark is strong save the file to a subsetresult dir'''

from os.path import join
from glob import glob 
from json import load
from collections import defaultdict

totalCount = 0
strongCount = 0
mutliV3Jump = 0
multiStepCount = 0
complexCount = 0
failDict = defaultdict(int)
strongestBM = (0,0)
TcMin = 1e100

directory = "Results/2LoopResults/Combined01SSS"
allFiles = glob(join(directory,'*.json'))

if len(allFiles) == 0:
    print("Empty directory, exiting")
    exit()
    
for filePointer in allFiles:
    fileName = filePointer.split('/')[1].lstrip().split(' ')[0]
    with open(filePointer, "r") as f:
        resultDic = load(f)
        totalCount +=1
        if resultDic["failureReason"]:
            failDict[resultDic["failureReason"]] +=1
        else:
            if len(resultDic["jumpsv3"])>1:
                mutliV3Jump += 1
                
            if resultDic["strong"] >0.6:
                strongestBM = (resultDic["strong"], resultDic["bmNumber"]) if resultDic["strong"] > strongestBM[0] else strongestBM
                strongCount += 1
                TcMin = resultDic["jumpsv3"][0][1] if resultDic["jumpsv3"][0][1] < TcMin else TcMin
                
            if resultDic["step"] > 1:
                multiStepCount += 1
                
            if resultDic["complexMin"]:
                complexCount += 1

print(f"Summary of the results in the directory '{directory}':")
print(f"The lowest Tc min is: {TcMin}")
print(f"The strongest BM is: {strongestBM}")
print(f"The number of benchmarks is: {totalCount}")
print(f"The number of strong benchmarks is: {strongCount}")
print(f"The number of mutli step phase transitions is: {multiStepCount}")
print(f"The number of failed benchmarks is: {failDict.items()}")
print(f"The mutli v3 jump count is: {mutliV3Jump}")
print(f"The number of benchmarks with a complex min is: {complexCount}")
