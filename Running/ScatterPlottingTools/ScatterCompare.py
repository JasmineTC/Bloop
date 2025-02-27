from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump
from matplotlib import pylab as plt
from numpy import transpose, unique, asarray
from matplotlib.cm import ScalarMappable
strengthList2Loop = []
bmInputList2Loop = []
TcList2Loop = []
for fileName in glob("2LoopResults/2Loop1GeV/*.json"):
    resultDic = load( open(fileName, "r") )
    if resultDic["strong"]:
        strengthList2Loop.append((resultDic['strong']))
        bmInputList2Loop.append((list(resultDic['bmInput'].values())))
        TcList2Loop.append(resultDic["jumpsv3"][0][1])


strengthList1Loop = []
bmInputList1Loop = []
TcList1Loop = []
for fileName in glob("1LoopResults/Results01GeV/*.json"):
    resultDic = load( open(fileName, "r") )
    if resultDic["strong"]:
        strengthList1Loop.append((resultDic['strong']))
        bmInputList1Loop.append((list(resultDic['bmInput'].values())))
        TcList1Loop.append(resultDic["jumpsv3"][0][1])

#Take the transpose so each list is one type of bm input
bmInputList1Loop = transpose(bmInputList1Loop)
bmInputList2Loop = transpose(bmInputList2Loop)
#Sort the bmInputs by order of strength, this is so the colour of the scatter plot is set by the strong PT at that point
strength1Loop, theta1Loop, gHDM1Loop, ms11Loop, delta121Loop, delta1c1Loop, deltac1Loop,_, Tc1Loop = zip(*sorted(zip(strengthList1Loop, *bmInputList1Loop, TcList1Loop)))
strength2Loop, theta2Loop, gHDM2Loop, ms12Loop, delta122Loop, delta1c2Loop, deltac2Loop,_, Tc2Loop = zip(*sorted(zip(strengthList2Loop, *bmInputList2Loop, TcList2Loop)))
#Finds each unique combination of 2d list, the first index each unique element shows up and the multiplicity of the element
strongest = max(strength1Loop + strength2Loop)
weakest = min(strength1Loop + strength2Loop)

norm=plt.Normalize(weakest,strongest)
plt.scatter(ms11Loop, theta1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\theta_{\\text{CPV}}$ ", labelpad = 5, fontsize = 12)
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopTheta_01")
plt.close()

plt.scatter(ms12Loop, theta2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\theta_{\\text{CPV}}$ ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopTheta_1")
plt.close()

plt.scatter(ms11Loop, gHDM1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$g_{\\text{hDM}}$ ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopGHDM_01")
plt.close()


plt.scatter(ms12Loop, gHDM2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$g_{\\text{hDM}}$ ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopGHDM_1")
plt.close()

plt.scatter(ms11Loop, delta121Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\delta_{12} \\  (GeV)$", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopDelta12_01")
plt.close()

plt.scatter(ms12Loop, delta122Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\delta_{12}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopDelta12_1")
plt.close()

plt.scatter(ms11Loop, delta1c1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\delta_{1c}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopDelta1c_01")
plt.close()

plt.scatter(ms12Loop, delta1c2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\delta_{1c}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopDelta1c_1")
plt.close()

plt.scatter(ms11Loop, deltac1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\delta_{c}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopDeltac_01")
plt.close()

plt.scatter(ms12Loop, deltac2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$\\delta_{c}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopDeltac_1")
plt.close()

plt.scatter(ms11Loop, Tc1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c$ (GeV)", labelpad = 5, fontsize = 12)
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopTcMs1_01")
plt.close()

plt.scatter(ms12Loop, Tc2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$m_{s1}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c $ (GeV) ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopTcMs1_1")
plt.close()

plt.scatter(gHDM1Loop, Tc1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$g_{\\text{hDM}}$ ", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c$ (GeV)", labelpad = 5, fontsize = 12)
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopTcgHDM_01")
plt.close()

plt.scatter(gHDM2Loop, Tc2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$g_{\\text{hDM}}$ ", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c $ (GeV) ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopTcgHDM_1")
plt.close()

plt.scatter(theta1Loop, Tc1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$\\theta_{\\text{CPV}}$ ", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c$ (GeV)", labelpad = 5, fontsize = 12)
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopTcTheta_01")
plt.close()

plt.scatter(theta2Loop, Tc2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$\\theta_{\\text{CPV}}$ ", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c $ (GeV) ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopTcTheta_1")
plt.close()

plt.scatter(delta121Loop, Tc1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$\\delta_{12}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c$ (GeV)", labelpad = 5, fontsize = 12)
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopTcDelta12_01")
plt.close()

plt.scatter(delta122Loop, Tc2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$\\delta_{12}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c $ (GeV) ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopTcDelta12_1")
plt.close()

plt.scatter(delta1c1Loop, Tc1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$\\delta_{1c}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c$ (GeV)", labelpad = 5, fontsize = 12)
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopTcDelta1c_01")
plt.close()

plt.scatter(delta1c2Loop, Tc2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$\\delta_{1c}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c $ (GeV) ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopTcDelta1c_1")
plt.close()

plt.scatter(deltac1Loop, Tc1Loop,s=4.2**2, c = strength1Loop, marker = "o", norm=norm)
plt.xlabel("$\\delta_{c}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c$ (GeV)", labelpad = 5, fontsize = 12)
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/1LoopTcDeltac_01")
plt.close()

plt.scatter(deltac2Loop, Tc2Loop,s=4.2**2, c = strength2Loop, marker = "o", norm=norm)
plt.xlabel("$\\delta_{c}$ (GeV)", labelpad = 5, fontsize = 12 )
plt.ylabel("$T_c $ (GeV) ", labelpad = 5, fontsize = 12 )
plt.colorbar(ticks = (0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.15), label = "strength")
plt.savefig("ScanFigures/2LoopTcDeltac_1")
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


