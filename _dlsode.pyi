"""
Python wrapper to the DLSODE solver from the ODEPACK family of solvers, originally
written in Fortran.
"""

from typing import Callable, Tuple

import numpy as np
from numpy.typing import ArrayLike

def dlsode(
    f: Callable,
    y: ArrayLike[np.float64],
    t: float,
    tout: float,
    rtol: float | ArrayLike[np.float64],
    atol: ArrayLike[np.float64],
    itask: int,
    istate: int,
    rwork: ArrayLike[np.float64],
    iwork: ArrayLike[np.int64],
    jac: Callable,
    mf: int,
    f_extra_args: Tuple = (),
    overwrite_y: int = 0,
) -> Tuple[ArrayLike, float, int]: ...
