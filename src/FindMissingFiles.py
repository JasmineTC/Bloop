from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load
import numpy as np

bmNumList = []
for filename in glob(join('Results', '*.json')):
   bmNumList.append(int(filename.split("_")[1])) ##The fastest way to extract the BM number from a file is to read it from the file name (saves loading the file)
   ## BM number is seperated by the rest of the string by the _ charact so spilt based on that character

bmNumArraySort = np.sort(bmNumList)
bmNumArrayDiff = np.diff( bmNumArraySort, append = 2 ) ##If there's a missing BM point then the sorted BM list will be [0,1,(missing), 3] taking the difference between each element
## and checking where the difference is >1 will find the missing elements
bMissingTag = bmNumArrayDiff >1

print(bmNumArraySort[bMissingTag]+1) ## Add 1 to get the actual mising bm
