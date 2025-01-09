from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump
from matplotlib import pylab as plt
from numpy import transpose, unique, asarray


strengthList = []
bmInputList = []

for fileName in glob('ResultsSSS/*.json'):
    resultDic = load( open(fileName, "r") )
    if resultDic["strong"]:
        strengthList.append((resultDic['strong']))
        bmInputList.append((list(resultDic['bmInput'].values())))
#Take the transpose so each list is one type of bm input
bmInputList = transpose(bmInputList)
#Sort the bmInputs by order of strength, this is so the colour of the scatter plot is set by the strong PT at that point
strength, theta, gHDM, ms1, delta12, delta1c, deltac,_ = zip(*sorted(zip(strengthList, *bmInputList)))
#Finds each unique combination of 2d list, the first index each unique element shows up and the multiplicity of the element
strongest = max(strengthList)
for idx, val in enumerate(ms1):
    plt.scatter(val, theta[idx],s = 4.2**2,c = strength[idx], marker = shape(strength[idx]),alpha = alphaNum(strength[idx]) , vmin=0.6, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$\\theta_{CPV}$")
plt.show()

for idx, val in enumerate(ms1):
    plt.scatter(val, gHDM[idx],  s = 4.2**2,c = strength[idx], marker = shape(strength[idx]) , vmin=0.6, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$gHDM$")
plt.show()

for idx, val in enumerate(ms1):
    plt.scatter(val, delta12 [idx], s = 4.2**2, c = strength[idx], marker = shape(strength[idx]) , vmin=0.6, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{12}$")
plt.show()


for idx, val in enumerate(ms1):
    plt.scatter(val, delta1c [idx], s = 4.2**2, c = strength[idx], marker = shape(strength[idx]) , vmin=0.6, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{1c}$")
plt.show()

for idx, val in enumerate(ms1):
    plt.scatter(val, deltac [idx], s = 4.2**2, c = strength[idx], marker = shape(strength[idx]), vmin=0.6, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{c}$")
plt.show()
