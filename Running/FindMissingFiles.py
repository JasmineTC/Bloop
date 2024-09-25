from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load
import numpy as np

bmNumList = []
for filename in glob(join('Results', '*.json')):
   bmNumList.append(int(filename.split("_")[1]))

bmNumArraySort = np.sort(bmNumList)
bmNumArrayDiff = np.diff( bmNumArraySort, append = 2 )
bMissingTag = bmNumArrayDiff >1

print(bmNumArraySort[bMissingTag])

