from os.path import join
from glob import glob  ## I think glob lets you do the * thingy
from json import load, dump
from matplotlib import pylab as plt
from numpy import transpose, unique, asarray

passedBm = []
v3NotGlobalMinBm = []
unBoundedBm = []
minFailedBm = []
complexBm = []
complexCount = 0
nonPertCount = 0

for fileName in glob("Results/*.json"):
    resultDic = load(open(fileName, "r"))
    if not resultDic["bIsPerturbative"]:
        nonPertCount += 1
    if resultDic["complexMin"]:
        complexCount += 1
    if resultDic["failureReason"] == None:
        passedBm.append((list(resultDic["bmInput"].values())))

    elif resultDic["failureReason"] == "v3NotGlobalMin":
        v3NotGlobalMinBm.append((list(resultDic["bmInput"].values())))

    elif resultDic["failureReason"] == "Unbounded":
        unBoundedBm.append((list(resultDic["bmInput"].values())))

    elif resultDic["failureReason"] == "v3NotGlobalMin":
        v3NotGlobalMinBm.append((list(resultDic["bmInput"].values())))

    elif resultDic["failureReason"] == "MinimisationFailed":
        minFailed.append((list(resultDic["bmInput"].values())))

    else:
        print(resultDic["failureReason"])
        exit()


print(f"The number of benchmark points that failed minimisation is: {len(minFailedBm)}")
print(
    f"The number of benchmark points that don't have v3 global min at 50GeV is: {len(v3NotGlobalMinBm)}"
)
print(
    f"The number of benchmark points that are unbounded at finite T is: {len(unBoundedBm)}"
)
print(f"The number of benchmark points that contain a complex  min is: {complexCount}")
print(f"The number of benchmark points are non pert is: {nonPertCount}")
# Take the transpose so each list is one type of bm input
passedBm = transpose(passedBm)
v3NotGlobalMinBm = transpose(v3NotGlobalMinBm)
unBoundedBm = transpose(unBoundedBm)
minFailedBm = transpose(minFailedBm)
complexBm = transpose(complexBm)
# theta, gHDM, ms1, delta12, delta1c, deltac,_ = passedBm

plt.scatter(passedBm[2], passedBm[0], s=4.2**2, c="green", label="okay")
plt.xlabel("$ms_1$")
plt.ylabel("$\\theta_{CPV}$")
plt.show()

plt.scatter(passedBm[2], passedBm[0], s=4.2**2, c="green", label="okay")
plt.scatter(
    v3NotGlobalMinBm[2], v3NotGlobalMinBm[0], s=4.2**2, c="blue", label="v3NotGlobal"
)
plt.scatter(unBoundedBm[2], unBoundedBm[0], s=4.2**2, c="black", label="unbounded")
plt.legend(loc=(0.03, 0.8))
plt.xlabel("$ms_1$")
plt.ylabel("$\\theta_{CPV}$")
plt.show()

plt.scatter(passedBm[2], passedBm[1], s=4.2**2, c="green", label="okay")
plt.xlabel("$ms_1$")
plt.ylabel("$gHDM$")
plt.show()

plt.scatter(passedBm[2], passedBm[1], s=4.2**2, c="green", label="okay")
plt.scatter(
    v3NotGlobalMinBm[2], v3NotGlobalMinBm[1], s=4.2**2, c="blue", label="v3NotGlobal"
)
plt.scatter(unBoundedBm[2], unBoundedBm[1], s=4.2**2, c="black", label="unbounded")
plt.legend(loc=(0.03, 0.8))
plt.xlabel("$ms_1$")
plt.ylabel("$gHDM$")
plt.show()

plt.scatter(passedBm[2], passedBm[3], s=4.2**2, c="green", label="okay")
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{12}$")
plt.show()


plt.scatter(passedBm[2], passedBm[3], s=4.2**2, c="green", label="okay")
plt.scatter(
    v3NotGlobalMinBm[2], v3NotGlobalMinBm[3], s=4.2**2, c="blue", label="v3NotGlobal"
)
plt.scatter(unBoundedBm[2], unBoundedBm[3], s=4.2**2, c="black", label="unbounded")
plt.legend(loc=(0.03, 0.8))
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{12}$")
plt.show()

plt.scatter(passedBm[2], passedBm[4], s=4.2**2, c="green", label="okay")
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{1c}$")
plt.show()


plt.scatter(passedBm[2], passedBm[4], s=4.2**2, c="green", label="okay")
plt.scatter(
    v3NotGlobalMinBm[2], v3NotGlobalMinBm[4], s=4.2**2, c="blue", label="v3NotGlobal"
)
plt.scatter(unBoundedBm[2], unBoundedBm[4], s=4.2**2, c="black", label="unbounded")
plt.legend(loc=(0.03, 0.8))
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{1c}$")
plt.show()

plt.scatter(passedBm[2], passedBm[5], s=4.2**2, c="green", label="okay")
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{c}$")
plt.show()


plt.scatter(passedBm[2], passedBm[5], s=4.2**2, c="green", label="okay")
plt.scatter(
    v3NotGlobalMinBm[2], v3NotGlobalMinBm[5], s=4.2**2, c="blue", label="v3NotGlobal"
)
plt.scatter(unBoundedBm[2], unBoundedBm[5], s=4.2**2, c="black", label="unbounded")
plt.legend(loc=(0.03, 0.8))
plt.xlabel("$ms_1$")
plt.ylabel("$\\delta_{c}$")
plt.show()
