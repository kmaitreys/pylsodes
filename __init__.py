"""
Python wrapper for DLSODES solver from Fortran ODEPACK family of
differential equation solvers.
"""

from ._dlsode import dlsode
from ._dlsodes import dlsodes


class Solver:
    def __init__(self, f, jac):
        self.f = f
        self.jac = jac
        self.f_params = ()
        self.jac_params = ()
        self._y = []
        self.issparse = False

    @property
    def y(self):
        return self._y

    def integrate(self, t, tout):
        if self.issparse:
            dlsodes()
        else:
            dlsode()

    def update(self): ...
