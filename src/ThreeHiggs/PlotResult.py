import matplotlib.pyplot as plt
import numpy as np


def makeFieldDimensionless(temp: list[float], field: list[float]) -> list[float]:
    return field / np.sqrt(temp)


def plotData(result: dict, index: int, loopOrder: int, filename: str):  # ->png
    if result["failureReason"]:
        return
    tempList = result["T"]

    v1List = makeFieldDimensionless(tempList, result["minimumLocation"][0])
    v2List = makeFieldDimensionless(tempList, result["minimumLocation"][1])
    v3List = makeFieldDimensionless(tempList, result["minimumLocation"][2])

    plt.plot(tempList, v1List, "v", label="$v_1$", color="orange")
    plt.plot(tempList, v2List, "^", label="$v_2$", color="black")
    plt.plot(tempList, v3List, ".", label="$v_3$", color="blue")
    plt.legend(loc="best")
    # plt.title(f"Field values at global minimum for benchmark {index} at loop order {loopOrder}")
    plt.ylabel(r"$\dfrac{v}{\sqrt{T}}$", rotation=0, labelpad=10)
    plt.xlabel("T (GeV)")
    plt.savefig(f"{filename}.png")
    plt.show()
