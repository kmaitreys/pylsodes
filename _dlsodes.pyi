"""
Python wrapper to the DLSODES solver from the ODEPACK family of solvers, originally
written in Fortran.
"""

from typing import Callable, Tuple

import numpy as np
from numpy.typing import ArrayLike

def dlsodes(
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
    jac_extra_args: Tuple = (),
) -> Tuple[ArrayLike, float, int]:
    """
    Description
    -----------
    Solves the initial value problem for ordinary differential equations using 
    the `DLSODES` solver from the `ODEPACK` family of ODE solver.

    `DLSODES` is a variant of the `DLSODE` package, and is intended for
    problems in which the Jacobian matrix `df/dy` has an arbitrary
    sparse structure (when the problem is stiff).

    Parameters
    ----------
    f : callable
        Right-hand side of the system of equations. The calling signature is `f(t, y, *f_extra_args)`.
    
    y : array_like
        Initial conditions on the dependent variables.
    
    t : float
        Initial value of the independent variable.
    
    tout : float
        Final value of the independent variable.
    
    rtol : float or array_like
        Relative tolerance(s) for the solver. If a scalar, the same tolerance is applied to all elements of `y`.
    
    atol : array_like
        Absolute tolerance(s) for the solver. Must have the same shape as `y`.
    
    itask : int
        Task indicator. The solver integrates from `t` to `tout` and returns the solution at `tout`.
        Default is `1`.
    
    istate : int
        State indicator. If `istate = 1`, the solver performs initialization steps.
        If `istate = 2`, the solver continues the integration. Default is `1`.
    
    rwork : array_like
        Real workspace array. Must have a length of at least `lrw`.
    
    iwork : array_like
        Integer workspace array. Must have a length of at least `liw`.
    
    jac : callable
        Jacobian matrix of the system. The calling signature is `jac(t, y, j, *jac_extra_args)`
        where `j` is the `j`-th column of the Jacobian matrix.
    
    mf : int
        Method flag, `mf = 10` for nonstiff (Adams) method, no Jacobian used,
        `mf = 121` for stiff (BDF) method, user-supplied sparse Jacobian,
        `mf = 222` for stiff method, internally generated sparse Jacobian

    Returns
    -------
    yout : array_like
        Solution at `tout`.
    
    tout : float
        Final value of the independent variable.
    
    istate : int
        State indicator. If `istate = 2`, the integration was successful.


    Authors
    -------

    -   Alan C. Hindmarsh \\
        Center for Applied Scientific Computing, L-561 \\
        Lawrence Livermore National Laboratory \\
        Livermore, CA 945


    -   Andrew H. Sherman \\
        J. S. Nolen and Associates \\
        Houston, TX 77084

    References
    ----------
     1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
         Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
         North-Holland, Amsterdam, 1983, pp. 55-64.

     2.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
         Yale Sparse Matrix Package: I. The Symmetric Codes,
         Int. J. Num. Meth. Eng., 18 (1982), pp. 1145-1151.

     3.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
         Yale Sparse Matrix Package: II. The Nonsymmetric Codes,
         Research Report No. 114, Dept. of Computer Sciences, Yale
         University, 1977.

    
    """
    ...
