'''---prints the number of SFOPT in a dir ---
Loads all the file names in the results directory 
Strips the Results/ part of the file name
Loads the benchmark information from the file
If the benchmark is strong save the file to a subsetresult dir'''
from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump
directory = "1LoopResults/MaximalCPV01"
totalCount = 0
strongCount = 0
mutliV3Jump = 0
multiStepCount = 0
complexCount = 0
from collections import defaultdict
failDict = defaultdict(int)

for filePointer in glob(join(directory, '*.json')):
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
                strongCount += 1
            if resultDic["step"] > 1:
                multiStepCount += 1
            if resultDic["complexMin"]:
                complexCount += 1
print(f"Summary of the results in the directory '{directory}':")
print(f"The number of benchmarks is: {totalCount}")
print(f"The number of strong benchmarks is: {strongCount}")
print(f"The number of mutli step phase transitions is: {multiStepCount}")
print(f"The number of failed benchmarks is: {failDict.items()}")
print(f"The mutli v3 jump count is: {mutliV3Jump}")
print(f"The number of benchmarks with a complex min is: {complexCount}")
