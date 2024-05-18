"""
Python wrapper to the DLSODES solver from the ODEPACK family of solvers, originally
written in Fortran.
"""

def dlsodes(
    y: list[float],
    t: float,
    tout: float,
    itol: int,
    rtol: float,
    atol: float,
    itask: int,
    istate: int,
    iopt: int,
    rwork: list[float],
    lrw: int,
    iwork: list[int],
    liw: int,
    jac: callable,
) -> None:
    """
    Solves the initial value problem for ordinary differential equations using the DLSODES solver.

    Parameters:
    y (list[float]): Array of initial conditions and output for solution vector.
    t (float): The initial value of the independent variable.
    tout (float): The next value of t at which a solution is desired.
    itol (int): An indicator for the type of error control.
    rtol (float): Relative error tolerance.
    atol (float): Absolute error tolerance.
    itask (int): An index specifying the task to be performed.
    istate (int): An index used for input and output to specify the state of the calculation.
    iopt (int): An indicator for whether optional inputs are being used.
    rwork (list[float]): Work array of length lrw.
    lrw (int): The length of rwork.
    iwork (list[int]): Work array of length liw.
    liw (int): The length of iwork.
    jac (callable): A function to compute the Jacobian matrix or None.

    Returns:
    None
    """
    ...
