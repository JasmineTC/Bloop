import matplotlib.pylab as plt
import numpy as np

class PlotResult():

    def __init__():
        None

    @staticmethod
    def RenormTemp(tempList : list, minimumList : list):# -> list
        for i, value in enumerate(minimumList):
            minimumList[i] =  minimumList[i]/np.sqrt(tempList[i])
        return minimumList

    @staticmethod
    def PlotData(result : list, index : int): #->png

        result = result.transpose()

        tempList = result[0]
        minimumList =  result[4]

        minimumList = PlotResult.RenormTemp(tempList, minimumList)

        plt.plot(tempList, minimumList, '.')
        plt.title(f"Results of bench mark {index}")
        plt.ylabel("Minimum")
        plt.xlabel("T")
        plt.savefig(f"Results/bm{index}.png")

        result.transpose()