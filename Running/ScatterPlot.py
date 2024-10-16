from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump
from matplotlib import pylab as plt
from numpy import transpose, unique

strengthList = []
bmInputList = [] 
count = 0
for fileName in glob(join('Results', '*.json')):
    with open(fileName, "r") as f:
        resultDic = load(f)
        if resultDic["strong"]:
            strengthList.append((resultDic['strong']))
            bmInputList.append((list(resultDic['bmInput'].values())))
#Take the transpose so each list is one type of bm input
bmInputList = transpose(bmInputList)

#Sort the bmInputs by order of strength, this is so the colour of the scatter plot is set by the strong PT at that point
strengthSorted, theta, gHDM, ms1, delta12, delta1c, deltac = zip(*sorted(zip(strengthList, *bmInputList)))
#Star producing figures
# plt.scatter(x, y, variable to set colour)
plt.colorbar(plt.scatter(ms1, theta, c=strengthSorted), label ="$\\frac{v_c}{T_c}$")
plt.xlabel("$ms_1$")
plt.ylabel("$\\theta_{CPV}$")
ms1ThetaTuple = tuple(zip(ms1, theta))
#Finds each unique combination of 2d list, the first index each unique element shows up and the multiplicity of the element
uniqueElement, uniqueIndex, multiplicity = unique(ms1ThetaTuple, return_index = True, return_counts= True, axis = 0)
#Labels each point with its multiplicity
for idx, element in enumerate(uniqueElement):
    plt.annotate(multiplicity[idx], (element[0]-1, element[1]+0.05))
#Adds the weakest point to the plot, just to the left of the orginal
for idx in uniqueIndex:
    plt.scatter(ms1ThetaTuple[idx][0] -1, ms1ThetaTuple[idx][1], c=strengthSorted[idx]) 
plt.show()

plt.colorbar(plt.scatter(ms1, gHDM, c=strengthSorted), label ="$\\frac{v_c}{T_c}$")
plt.xlabel("$ms_1$")
plt.ylabel("$gHDM$")
ms1gHDMTuple = tuple(zip(ms1, gHDM))
uniqueElement, uniqueIndex, multiplicity = unique(ms1gHDMTuple, return_index = True, return_counts= True, axis = 0)
for idx, element in enumerate(uniqueElement):
    plt.annotate(multiplicity[idx], (element[0]-1, element[1]+0.005))
for idx in uniqueIndex:
    plt.scatter(ms1gHDMTuple[idx][0] -1, ms1gHDMTuple[idx][1], c=strengthSorted[idx]) 
plt.show()


plt.colorbar(plt.scatter(ms1, delta12, c=strengthSorted), label ="$\\frac{v_c}{T_c}$")
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{12}$")
ms1delta12Tuple = tuple(zip(ms1, delta12))
uniqueElement, uniqueIndex, multiplicity = unique(ms1delta12Tuple, return_index = True, return_counts= True, axis = 0)
for idx, element in enumerate(uniqueElement):
    plt.annotate(multiplicity[idx], (element[0]-1, element[1]+0.05))
for idx in uniqueIndex:
    plt.scatter(ms1delta12Tuple[idx][0] -1, ms1delta12Tuple[idx][1], c=strengthSorted[idx]) 
plt.show()


plt.colorbar(plt.scatter(ms1, delta1c, c=strengthSorted), label ="$\\frac{v_c}{T_c}$")
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{1c}$")
ms1delta1cTuple = tuple(zip(ms1, delta1c))
uniqueElement, uniqueIndex, multiplicity = unique(ms1delta1cTuple, return_index = True, return_counts= True, axis = 0)
for idx, element in enumerate(uniqueElement):
    plt.annotate(multiplicity[idx], (element[0]-1, element[1]+0.05))
for idx in uniqueIndex:
    plt.scatter(ms1delta1cTuple[idx][0] -1, ms1delta1cTuple[idx][1], c=strengthSorted[idx]) 
plt.show()


plt.colorbar(plt.scatter(ms1, deltac, c=strengthSorted), label ="$\\frac{v_c}{T_c}$")
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{c}$")
ms1deltacTuple = tuple(zip(ms1, deltac))
uniqueElement, uniqueIndex, multiplicity = unique(ms1deltacTuple, return_index = True, return_counts= True, axis = 0)
for idx, element in enumerate(uniqueElement):
    plt.annotate(multiplicity[idx], (element[0]-1, element[1]+0.05))
for idx in uniqueIndex:
    plt.scatter(ms1deltacTuple[idx][0] -1, ms1deltacTuple[idx][1], c=strengthSorted[idx]) 
plt.show()

