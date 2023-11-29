import numpy as np
import numpy.typing as npt

#------------------ 3D integrals 


""" UV-finite integral for 1-loop effective potential.
    J3(m^2) = 1/2 * int_p ln(p^2 + m^2)
    Result is proportional to m^3, so can be complex...
"""
def J3(msq: npt.ArrayLike) -> npt.ArrayLike:
    return -1. / (12.*np.pi) * msq * np.sqrt(msq + 0j)