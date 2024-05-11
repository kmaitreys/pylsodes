"""
This module acts as a wrapper for the D/LSODES solver from the ODEPACK library.
Adapted from `odespy` by Hans Petter Langtangen originally written in Python 2.
"""

import numpy as np


class LSODES:
    _description = """
    LSODES is a Python wrapper for the D/LSODES solver from the ODEPACK library.
    D/LSODES is a solver for stiff differential equations, with sparse Jacobians.
    """

    _optional_parameters = [
        "order",
        "moss",
        "seth",
        "jac_column",
        "ia",
        "ja",
        "jac_column_f77",
    ]

    def set_iter_method(self):
        with_jac_column = hasattr(self, "jac_column") or hasattr(self, "jac_column_f77")
        with_ia_ja = hasattr(self, "ia") and hasattr(self, "ja")
        if not hasattr(self, "moss"):
            if with_ia_ja:
                self.moss = 0
            elif with_jac_column:
                self.moss = 1
            else:
                self.moss = 2
        if not hasattr(self, "corrector_iter_method"):
            self.iter_method = int(with_jac_column)

    def set_iwork_rwork(self):
        """
        Initialization of work arrays with caculated length and optional inputs.
        In ODEPACK, "iwork" & "rwork" should be initialized with the specific
        optional parameters in all the solvers.
        "liw" & "lrw" represented the length requirement of work arrays.
        Specially, in Dlsodes, ia & ja should be attached to iwork_in.
        """

        # initialize integer work array (iwork)
        self.iwork_in = [0] * 30
        if (self.moss == 0) and hasattr(self, "ia"):
            self.iwork_in += list(self.ia) + list(self.ja)
        nnz = (
            len(getattr(self, "ja", [])) if hasattr(self, "ja") else self.neq**2 / 2
        )  # default value
        self.liw_min = len(self.iwork_in)
        for index in self._iwork_index:
            self.iwork_in[index] = getattr(self, self._iwork_index[index], 0)

        # calculate the minimum length of float work array (rwork)
        maxord = 5 if self.adams_or_bdf == "bdf" else 12
        self.lrw_min = 20 + self.neq * (maxord + 4)
        lrw_arg = [(0, 0, 0, 0), (0, 2, 2, 9), (0, 2, 2, 10), (self.neq + 2, 0, 0, 0)][
            self.iter_method
        ]
        self.lrw_min += (
            lrw_arg[0]
            + nnz * lrw_arg[1]
            + self.neq * lrw_arg[2]
            + ((nnz + lrw_arg[3] * self.neq) / 2)
        )

        # Initializereal input work arrays
        lrw_in = max(self._rwork_index.keys()) + 1
        self.rwork_in = np.zeros(lrw_in, float)
        for index in self._rwork_index:
            self.rwork_in[index] = getattr(self, self._rwork_index[index], 0.0)

    def set_iopt(self):
        self.iopt = int(self.iwork_in[4:7] > 0 or any(self.rwork_in[4:7] > 0))

    def initialize_for_solve(self):
        self.set_iter_method()
        self.func_wrapper()
        self.mf = (
            self.iter_method
            + (1 + (self.adams_or_bdf == "bdf")) * 10
            + getattr(self, "moss", 0) * 100
        )
        self.set_iwork_rwork()
        self.set_iopt()
        self.set_ydoti()
        self.itol = (
            (not isinstance(self.rtol, (int, float))) * 2
            + (not isinstance(self.atol, (int, float)))
            + 1
        )
        

    def set_ydoti(self):
        """
        ``ydoti`` is an array used in linearly implicit solvers.
        It has to be extended if its length is smaller than neq.
        """
        ydoti = getattr(self, "ydoti", [])
        # flag to indicate whether ydoti is supplied
        self.ydoti_flag = int(ydoti != [])
        # extend the length of 'ydoti' to neq
        self.ydoti = list(ydoti) + [0.0] * (self.neq - len(ydoti))
