'''---prints the number of SFOPT in a dir ---
Loads all the file names in the results directory 
Strips the Results/ part of the file name
Loads the benchmark information from the file
If the benchmark is strong save the file to a subsetresult dir'''
from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump
directory = "1LoopResults/Results01GeV/"
strongCount = 0
errorCount = 0
multiStepCount = 0
complexCount = 0
failCount = 0

for filePointer in glob(join(directory, '*.json')):
    fileName = filePointer.split('/')[1].lstrip().split(' ')[0]
    with open(filePointer, "r") as f:
        resultDic = load(f)
        if len(resultDic["jumpsv3"])>1:
            errorCount += 1
        if resultDic["failureReason"]:
            failCount +=1 
        if resultDic["strong"]:
            strongCount += 1
        if resultDic["step"] > 1:
            multiStepCount += 1
            #print(resultDic)
        if resultDic["complexMin"]:
            complexCount += 1

print(f"The number of strong benchmarks in the {directory} directory is: {strongCount}")
print(f"The error count is: {errorCount}")
print(f"The number of mutli step phase transitions is: {multiStepCount}")
print(f"The numer of benchmarks with a complex min is: {complexCount}")
print(f"The numer of failed benchmarks is: {failCount}")
