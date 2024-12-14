from numba import njit
import numpy as np
@njit
def diagonalizeNumba(matrices, matrixNumber, matrixSize, T):
    ##Gives complex cast warning
    subEigenValues = np.empty( (matrixNumber, matrixSize) ) 
    subRotationMatrix = np.empty( (matrixNumber, matrixSize, matrixSize) )
    for idx, matrix in enumerate(matrices):
         subEigenValues[idx], subRotationMatrix[idx] = np.linalg.eigh(matrix)
    return subEigenValues*T**2, subRotationMatrix

