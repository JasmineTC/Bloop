"""---prints the number of SFOPT in a dir ---
Loads all the file names in the results directory
Strips the Results/ part of the file name
Loads the benchmark information from the file
If the benchmark is strong save the file to a subsetresult dir"""

from os.path import join
from glob import glob
from json import load
from statistics import mean

directory1 = "1LoopResults/Combined01SSS/"
directory2 = "2LoopResults/Combined01SSS/"
becomeWeakCount = 0
pdDifP = []
pdDifM = []
stronger = []
Tc1Loop = []
Tc2Loop = []
for filePointer in glob(join(directory2, "*.json")):
    fileName = filePointer.split("/")[2]
    Loop2 = load(open(filePointer, "r"))
    Loop1 = load(open(directory1 + fileName, "r"))
    strength2Loop = Loop2["strong"]
    strength1Loop = Loop1["strong"]
    Tc1Loop.append(Loop1["jumpsv3"][0][1])
    Tc2Loop.append(Loop2["jumpsv3"][0][1])
    if strength2Loop > strength1Loop:
        stronger.append(Loop2["bmNumber"])

    if not strength2Loop:
        becomeWeakCount += 1
        continue
    pd = 100 * (strength2Loop - strength1Loop) / strength2Loop
    if pd >= 0:
        pdDifP.append(pd)
    if pd < 0:
        pdDifM.append(pd)
print("Comparing 2 loop scan to 1 loop scan (0.1 GeV step)")
print(f"Total BM at 2 loop:  {becomeWeakCount + len(pdDifP) + len(pdDifM)}")
print(
    f"Number of benchmarks that were strong at one loop but aren't at 2 loop: {becomeWeakCount}"
)
print(f"Number of benchmarks that become stronger at 2 loop: {len(pdDifP)}")
print(
    f"Number of benchmarks that become weaker at 1 loop: {len(pdDifM) + becomeWeakCount}"
)
print(f"Mean percentage difference of benchmarks that become stronger: {mean(pdDifP)}")
print(f"--- THE BELOW IS AN UNDER ESTIMATE AS I DON'T TRACK STRENGTHS BELOW 0.3 ---")
print(f"Mean percentage difference of benchmarks that become weaker: {mean(pdDifM)}")
print(f"Mean percentage difference overall: {mean(pdDifM + pdDifP)}")
meanTc1Loop = mean(Tc1Loop)
meanTc2Loop = mean(Tc2Loop)
print(f"Mean Tc at 1 loop: {meanTc1Loop}")
print(f"Mean Tc at 2 loop: {meanTc2Loop}")
print(f"Percentage difference: {100 * (meanTc2Loop - meanTc1Loop) / meanTc2Loop}")
