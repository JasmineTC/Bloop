from glob import glob  ## I think glob lets you do the * thingy
from json import load
from matplotlib import pylab as plt
from numpy import transpose


strengthList = []
bmInputList = []

for fileName in glob("ResultsSSS/*.json"):
    resultDic = load(open(fileName, "r"))
    if resultDic["strong"]:
        strengthList.append((resultDic["strong"]))
        bmInputList.append((list(resultDic["bmInput"].values())))
# Take the transpose so each list is one type of bm input
bmInputList = transpose(bmInputList)
# Sort the bmInputs by order of strength, this is so the colour of the scatter plot is set by the strong PT at that point
strength, theta, gHDM, ms1, delta12, delta1c, deltac, _ = zip(
    *sorted(zip(strengthList, *bmInputList))
)
# Finds each unique combination of 2d list, the first index each unique element shows up and the multiplicity of the element
strongest = max(strengthList)
weakest = min(strengthList)
plt.scatter(ms1, theta, s=4.2**2, c=strength, vmin=weakest, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$\\theta_{CPV}$")
plt.show()

plt.scatter(ms1, gHDM, s=4.2**2, c=strength, vmin=weakest, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$gHDM$")
plt.show()

plt.scatter(ms1, delta12, s=4.2**2, c=strength, vmin=weakest, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{12}$")
plt.show()


plt.scatter(ms1, delta1c, s=4.2**2, c=strength, vmin=weakest, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{1c}$")
plt.show()

plt.scatter(ms1, deltac, s=4.2**2, c=strength, vmin=weakest, vmax=strongest)
plt.colorbar().solids.set(alpha=1)
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{c}$")
plt.show()
