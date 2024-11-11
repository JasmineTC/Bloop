'''Creates a txt file containing all the bm numbers that are strong, 
to passed to Benchmarks/ScanSubSet.py '''

from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump

strongBmList = []

for fileName in glob(join('Results', '*.json')):
    with open(fileName, "r") as f:
        resultDic = load(f)
        if resultDic["strong"]:
            strongBmList.append(resultDic["bmNumber"])

print(f"The number of benchmark files tagged as strong is: {len(strongBmList)}")

with open("Benchmarks/StrongBMList.txt", "w") as fp:
    dump(strongBmList, fp)
