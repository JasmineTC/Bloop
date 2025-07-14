from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump
from matplotlib import pylab as plt
from numpy import transpose, unique, asarray
from matplotlib.cm import ScalarMappable

def loadData(directoryList):
    strength, bmInput, Tc, bmNumber = [], [], [], []
    for directory in directoryList:
        allFileNames = glob(directory+"/*.json")
        if len(allFileNames) == 0:
            print(f"{directory} is empty, exiting")
            exit()
        for fileName in allFileNames:
            resultDic = load( open(fileName, "r") )
            if resultDic["strong"]>0.6:
                strength.append((resultDic['strong']))
                bmInput.append((list(resultDic['bmInput'].values())))
                Tc.append(resultDic["jumpsv3"][0][1])
                bmNumber.append(resultDic["bmNumber"])
    return strength, bmInput, Tc, bmNumber

strength2Loop, bmInput2Loop, Tc2Loop, bmNumber2Loop = loadData(["2LoopResults/Combined01SSS"] )
strength1Loop, bmInput1Loop, Tc1Loop, bmNumber1Loop = loadData(["1LoopResults/Combined01SSS"] )

#Sort the bmInputs by order of strength, this is so the colour of the scatter plot is set by the strong PT at that point, transpose taken so each row is just on variable type
strength2Loop, theta2Loop, gHDM12Loop, ms12Loop, delta122Loop, delta1c2Loop, deltac2Loop , _, Tc2Loop, bmNumber2Loop = zip(*sorted(zip(strength2Loop, *transpose(bmInput2Loop), Tc2Loop, bmNumber2Loop)))
strength1Loop, theta1Loop, gHDM11Loop, ms11Loop, delta121Loop, delta1c1Loop, deltac1Loop , _, Tc1Loop, bmNumber1Loop = zip(*sorted(zip(strength1Loop, *transpose(bmInput1Loop), Tc1Loop, bmNumber1Loop)))
#Finds each unique combination of 2d list, the first index each unique element shows up and the multiplicity of the element
strongest = max(strength1Loop + strength2Loop)
weakest = 0.6
norm=plt.Normalize(weakest,strongest)
from numpy import pi

#for index, theta in enumerate(theta2Loop):
#    if theta < 1.08*pi/2 and strength2Loop[index] > 0.8:
#        print(2*theta/pi)
#        print(strength2Loop[index])
#        print(bmNumber2Loop[index])
#        print()
    
inputData = ((Tc2Loop, ms12Loop, theta2Loop, gHDM12Loop, delta122Loop, delta1c2Loop, deltac2Loop, strength2Loop), 
             (Tc1Loop, ms11Loop, theta1Loop, gHDM11Loop, delta121Loop, delta1c1Loop, deltac1Loop, strength1Loop))

from pathlib import Path
ylabels = ["$m_{s1}$ (GeV)", "$\\theta_{\\text{CPV}}$","$g_{\\text{hDM}}$", "$\\delta_{12} \\  (GeV)$", "$\\delta_{1c} \\  (GeV)$","$\\delta_{c} \\  (GeV)$"]
fileNames = ["ms101","Theta01","GHDM01","Delta1201","Delta1c01","Deltac01"]
loopList = ["2Loop", "1Loop"]
for idxL, dataRow in enumerate(inputData):
    #Makes plots of ms1 vs other stuff
    for idxI, data in enumerate(dataRow[2:-1]):
       #+1 to skip over Tc 
       idxI+=1
       plt.scatter(dataRow[1], data, s=4.2**2, c = dataRow[-1], marker = "o", norm=norm)
       plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
       plt.ylabel(ylabels[idxI], labelpad = 5, fontsize = 12)
       plt.colorbar(ticks = (0.60, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
       fileDir = f"TestFigures/{loopList[idxL]}"
       Path(fileDir).mkdir(parents=True, exist_ok = True)
       plt.savefig(fileDir+f"/{fileNames[idxI]}")
       plt.close()
    #Makes plots of Tc vs other stuff
    for idxI, data in enumerate(dataRow[1:-1]):
        plt.scatter(data, dataRow[0], s=4.2**2, c = dataRow[-1], marker = "o", norm=norm)
        plt.xlabel(ylabels[idxI], labelpad = 5, fontsize = 12 )
        plt.ylabel("$T_c$ (GeV)", labelpad = 5, fontsize = 12)
        plt.colorbar(ticks = (0.60, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
        fileDir = f"TestFigures/{loopList[idxL]}"
        Path(fileDir).mkdir(parents=True, exist_ok = True)
        plt.savefig(fileDir+f"/Tc{fileNames[idxI]}")
        plt.close()
    
exit()

items = 50
plt.scatter(ms12Loop[-items:], theta2Loop[-items:],s=4.6**2, c = strength2Loop[-items:], marker = "o", norm=norm)
for bm in bmNumber2Loop[-items:]:
    index = bmNumber1Loop.index(bm)
    plt.scatter(ms11Loop[index]*1.005, theta1Loop[index], s=4.6**2, c = strength1Loop[index], marker = "*", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\theta_{\\text{CPV}}$ ", labelpad = 5, fontsize = 12)
plt.colorbar(ticks = (0.60, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/TEST")
plt.close()

##Puts 1 and 2 loop into one figure
#f, axes = plt.subplots(nrows = 1, ncols = 2)
#norm=plt.Normalize(weakest,strongest)
#sc = axes[0].scatter(ms11Loop, theta1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
#axes[1].scatter(ms12Loop, theta2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
#for i in range(len(axes)):
#    axes[i].set_xlabel("$m_{s1}$ (GeV)", labelpad = 5)
#    axes[i].set_ylabel("$\\theta_{\\text{cpv}}$ ", labelpad = 5)
#f.set_figwidth(15)
#cbar_ax = f.add_axes([.92, .15, 0.02, 0.7])
#
#f.colorbar(sc, cax=cbar_ax, ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05), label = "strength")
#
#plt.show()

