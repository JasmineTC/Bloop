from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump
from matplotlib import pylab as plt
from numpy import transpose, unique, asarray

def shape(x:float) -> str:
    if x < 0.7:
        return 'o'
    elif x < 0.8:
        return '^'
    elif x < 0.9:
        return 'D'
    elif x < 1:
        return 's'
    return 'h'

strengthList = []
bmInputList = [] 
for fileName in glob(join('SubSetResults', '*.json')):
    with open(fileName, "r") as f:
        resultDic = load(f)
        if resultDic["strong"]:
            strengthList.append((resultDic['strong']))
            bmInputList.append((list(resultDic['bmInput'].values())))
#Take the transpose so each list is one type of bm input
bmInputList = transpose(bmInputList)
#Sort the bmInputs by order of strength, this is so the colour of the scatter plot is set by the strong PT at that point
strength, theta, gHDM, ms1, delta12, delta1c, deltac = zip(*sorted(zip(strengthList, *bmInputList)))

zipped = tuple(zip(ms1, theta))
#Finds each unique combination of 2d list, the first index each unique element shows up and the multiplicity of the element
uniqueElement, firstUniqueIndex, multiplicity = unique(zipped, return_index = True, return_counts= True, axis = 0)
lastUniqueIndex = [len(zipped) - zipped[::-1].index(tuple(element)) -1 for element in uniqueElement]  
#Plots the last unqiue occurance of a index (this is to get the strongest point) and annonates the point with its multiplicity
for idx, value in enumerate(lastUniqueIndex):
    plt.scatter(ms1[value], theta[value],  c = strength[value], marker= shape(strength[value]), vmin=0.6, vmax=1.1)
    plt.annotate(multiplicity[idx], (ms1[value]-2, theta[value]+0.1)) # --Note the off set for the label
#Plots the first unqiue occurance of a index (gets the weakest point) --Note the off set in the x axis
for idx in firstUniqueIndex:
    plt.scatter(ms1[idx]-1.2, theta[idx],  c = strength[idx], marker= shape(strength[idx]), vmin=0.6, vmax=1.1)
plt.colorbar()
plt.xlabel("$ms_1$")
plt.ylabel("$\\theta_{CPV}$")
plt.show()

zipped = tuple(zip(ms1, gHDM))
uniqueElement, firstUniqueIndex, multiplicity = unique(zipped, return_index = True, return_counts= True, axis = 0)
lastUniqueIndex = [len(zipped) - zipped[::-1].index(tuple(element)) -1 for element in uniqueElement]  
for idx, value in enumerate(lastUniqueIndex):
    plt.scatter(ms1[value], gHDM[value],  c = strength[value], marker= shape(strength[value]), vmin=0.6, vmax=1.1)
    plt.annotate(multiplicity[idx], (ms1[value]-2, gHDM[value]+0.01))
for idx in firstUniqueIndex:
    plt.scatter(ms1[idx]-1.2, gHDM[idx],  c = strength[idx], marker= shape(strength[idx]), vmin=0.6, vmax=1.1)
plt.colorbar()
plt.xlabel("$ms_1$")
plt.ylabel("$gHDM$")
plt.show()

zipped = tuple(zip(ms1, delta12))
uniqueElement, firstUniqueIndex, multiplicity = unique(zipped, return_index = True, return_counts= True, axis = 0)
lastUniqueIndex = [len(zipped) - zipped[::-1].index(tuple(element)) -1 for element in uniqueElement]  
for idx, value in enumerate(lastUniqueIndex):
    plt.scatter(ms1[value], delta12[value],  c = strength[value], marker= shape(strength[value]), vmin=0.6, vmax=1.1)
    plt.annotate(multiplicity[idx], (ms1[value]-2, delta12[value]+0.1))
for idx in firstUniqueIndex:
    plt.scatter(ms1[idx]-1.2, delta12[idx],  c = strength[idx], marker= shape(strength[idx]), vmin=0.6, vmax=1.1)
plt.colorbar()
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{12}$")
plt.show()


zipped = tuple(zip(ms1, delta1c))
uniqueElement, firstUniqueIndex, multiplicity = unique(zipped, return_index = True, return_counts= True, axis = 0)
lastUniqueIndex = [len(zipped) - zipped[::-1].index(tuple(element)) -1 for element in uniqueElement]  
for idx, value in enumerate(lastUniqueIndex):
    plt.scatter(ms1[value], delta1c[value],  c = strength[value], marker= shape(strength[value]), vmin=0.6, vmax=1.1)
    plt.annotate(multiplicity[idx], (ms1[value]-2, delta1c[value]+0.1))
for idx in firstUniqueIndex:
    plt.scatter(ms1[idx]-1.2, delta1c[idx],  c = strength[idx], marker= shape(strength[idx]), vmin=0.6, vmax=1.1)
plt.colorbar()
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{1c}$")
plt.show()

zipped = tuple(zip(ms1, deltac))
uniqueElement, firstUniqueIndex, multiplicity = unique(zipped, return_index = True, return_counts= True, axis = 0)
lastUniqueIndex = [len(zipped) - zipped[::-1].index(tuple(element)) -1 for element in uniqueElement]  
for idx, value in enumerate(lastUniqueIndex):
    plt.scatter(ms1[value], deltac[value],  c = strength[value], marker= shape(strength[value]), vmin=0.6, vmax=1.1)
    plt.annotate(multiplicity[idx], (ms1[value]-2, deltac[value]+0.1))
for idx in firstUniqueIndex:
    plt.scatter(ms1[idx]-1.2, deltac[idx],  c = strength[idx], marker= shape(strength[idx]), vmin=0.6, vmax=1.1)
plt.colorbar()
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{c}$")
plt.show()

