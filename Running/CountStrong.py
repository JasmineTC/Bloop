'''---prints the number of SFOPT in a dir ---
Loads all the file names in the results directory 
Strips the Results/ part of the file name
Loads the benchmark information from the file
If the benchmark is strong save the file to a subsetresult dir'''
from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump
directory = "2LoopResults/Filler01"
totalCount = 0
strongCount = 0
mutliV3Jump = 0
multiStepCount = 0
complexCount = 0
from collections import defaultdict
failDict = defaultdict(int)

TcMin = 1e100
for filePointer in glob(join(directory, '*.json')):
    fileName = filePointer.split('/')[1].lstrip().split(' ')[0]
    with open(filePointer, "r") as f:
        resultDic = load(f)
        totalCount +=1
        if resultDic["failureReason"]:
            failDict[resultDic["failureReason"]] +=1
        else:
            if len(resultDic["jumpsv3"])>1:
                print(resultDic)
                mutliV3Jump += 1
            if resultDic["strong"] >0.6:
                strongCount += 1
                TcMin = resultDic["jumpsv3"][0][1] if resultDic["jumpsv3"][0][1] < TcMin else TcMin
            if resultDic["step"] > 1:
                multiStepCount += 1
                if resultDic["strong"]>0.6:
                    print(resultDic["strong"])
            if resultDic["complexMin"]:
                complexCount += 1
print(f"Summary of the results in the directory '{directory}':")
print(TcMin)
print(f"The number of benchmarks is: {totalCount}")
print(f"The number of strong benchmarks is: {strongCount}")
print(f"The number of mutli step phase transitions is: {multiStepCount}")
print(f"The number of failed benchmarks is: {failDict.items()}")
print(f"The mutli v3 jump count is: {mutliV3Jump}")
print(f"The number of benchmarks with a complex min is: {complexCount}")
