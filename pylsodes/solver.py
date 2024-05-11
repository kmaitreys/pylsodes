"""
This module acts as a wrapper for the D/LSODES solver from the ODEPACK library.
Adapted from `odespy` by Hans Petter Langtangen originally written in Python 2.
"""

import pprint

import numpy as np

_parameters = dict(
    f=dict(help="Right-hand side ``f(u,t)`` defining the ODE.", type=callable),
    f_args=dict(
        help="Extra positional arguments to f: ``f(u, t, *f_args, **f_kwargs).``",
        type=(tuple, list, np.ndarray),
        default=(),
    ),
    f_kwargs=dict(
        help="Extra keyword arguments to f: ``f(u, t, *f_args, **f_kwargs)``.",
        type=dict,
        default={},
    ),
    complex_valued=dict(help="True if f is complex valued.", default=False, type=bool),
    f_is_linear=dict(help="True if f(u,t) is linear in u.", type=bool, default=False),
    jac=dict(
        help="Jacobian of right-hand side function f (df/du).",
        default=None,
        type=callable,
    ),
    jac_args=dict(
        help="Extra positional arguments to jac: ``jac(u, t, *jac_args,"
        "**jac_kwargs)``.",
        type=(tuple, list),
        default=(),
    ),
    jac_kwargs=dict(
        help="Extra keyword arguments to jac: ``jac(u, t, *jac_args,"
        "**jac_kwargs)``.",
        type=dict,
        default={},
    ),
    h_in_fd_jac=dict(
        help="h in finite difference approximation of the Jacobian.",
        default=1e-4,
        type=float,
    ),
    verbose=dict(
        help="Integer reflecting output of intermediate quantities.",
        default=0,
        type=int,
    ),
    disk_storage=dict(
        help="Indicates whether u is stored in memory or in file. "
        'If string, it is the filename; if False or "", u is '
        "kept in memory; if True, a default filename tmp_odspy.dat "
        "is used.",
        default=False,
        type=(str, bool),
    ),
    u_exact=dict(
        help="Function of t returning exact solution.", default=None, type=callable
    ),
    start_method=dict(
        help="Method for the first steps in multi-step solvers.",
        default="RK2",
        type=str,
    ),
    nonlinear_solver=dict(
        help="Standard Newton or Picard nonlinear solver, or Picard2 (u*f(u_,t)/u_ implicit trick).",
        default="Picard",
        type=str,
        range=("Newton", "Picard", "Picard2"),
    ),
    eps_iter=dict(
        help="Max error tolerance in nonlinear solver", default=1e-4, type=float
    ),
    max_iter=dict(
        help="Max no of iterations in nonlinear solver", default=50, type=int
    ),
    g=dict(
        help="Constraint function of (u, t) in differential-algebraic systems.",
        type=callable,
    ),
    ng=dict(help="No of components in constraint function g.", type=int),
    theta=dict(
        help='Weight in [0,1] used for "theta-rule" finite difference approx.',
        default=0.5,
        type=(int, float),
        range=[0, 1],
    ),
    # Parameters for adaptive methods
    atol=dict(
        help="Absolute tolerance for solution.",
        type=(float, list, tuple, np.ndarray),
        default=1e-8,
    ),
    rtol=dict(
        help="Relative tolerance for solution.",
        type=(list, tuple, np.ndarray, float),
        default=1e-6,
    ),
    min_step=dict(help="Minimum step size for an adaptive algorithm.", type=float),
    max_step=dict(help="Maximum step size for an adaptive algorithm.", type=float),
    first_step=dict(
        help="Suggested first time step size for an adaptive algorithm.", type=float
    ),
    solver=dict(
        help="Name of solver class in solvers that need an extra solver "
        "(e.g., AdaptiveResidual).",
        default="RK4",
        type=str,
    ),
    butcher_tableau=dict(
        help="2d-array which contains the butcher table for user-supplied "
        "Runge-Kutta method. (n,n) array for 1-level Runge-Kutta "
        "methods.(n+1,n) array for 2-level Runge-Kutta methods.",
        type=np.ndarray,
    ),
    # vode parameters
    adams_or_bdf=dict(
        help='Method in vode or solvers in odepack: "adams" or "bdf".',
        type=str,
        default="adams",
        range=["adams", "bdf"],
    ),
    order=dict(
        help="Maximum order used by the integrator "
        '(<= 12 for "adams", <= 5 for "bdf").',
        type=int,
        default=4,
    ),
    nsteps=dict(
        help="Max no of internal solver steps per time step.", type=int, default=1000
    ),
    method_order=dict(
        help="Method order for user-defined method if known."
        "A integer for 1-level methods, or a pair of   "
        "integer for 2-levels methods.",
        type=(int, tuple, list, np.ndarray),
    ),
    # beta, ifactor and dfactor are intended for adaptive Dormand&Prince
    # methods like dopri5 or dop853 in scipy
    beta=dict(
        help="Beta argument for stabilized step size control in "
        "Dormand&Prince methods from scipy",
        type=float,
    ),
    ifactor=dict(
        help="Maximum factor for increasing the step size", type=float, default=2
    ),
    dfactor=dict(
        help="Maximum factor for decreasing the step size", type=float, default=0.5
    ),
    safety=dict(help="Safety factor on new step selection", default=0.9, type=float),
    # odelab parameters
    odelab_solver=dict(
        help="Name of Solver class in odelab", default="RungeKutta4", type=str
    ),
    # Vode_PyDS parameters
    init_step=dict(help="Fixed step size for time mesh.", type=float),
    strictdt=dict(
        help="Uniform time mesh vs exact dt spacings", type=bool, default=True
    ),
    stiff=dict(help="Boolean flag to indicate stiffness.", type=bool),
    use_special=dict(help="Switch for using special times", type=bool),
    specialtimes=dict(
        help="List of special times to use during iteration",
        type=lambda float_seq: np.asarray(
            map(lambda x: isinstance(x, float), float_seq)
        ).all(),
    ),
    ode_method=dict(
        help='solver type: "adams" or "bdf"',
        alias="method",
        type=str,
        default="adams",
        range=("adams", "bdf"),
    ),
    relaxation=dict(
        help="relaxation argument (r): new_solution = r*solution + "
        "(1-r)*old_solution",
        default=1.0,
        type=float,
    ),
    # parameters for Jacobian
    jac_banded=dict(
        help="""\
Banded Jacobian matrix: ``jac_banded(u, t, ml, mu)``.
``ml`` and ``mu`` are the number of lower and upper
diagonals. The returned rectangular array should have
shape ``neq, ml+mu+1``.""",
        type=callable,
    ),
    jac_constant=dict(
        help="Flag to show whether Jacobian is constant, 0 (false) or 1 (true)",
        default=0,
        type=int,
    ),
    # parameters for linearly implicit ODE solvers: Lsodi, Lsoibt, Lsodis
    res=dict(
        help="User-supplied function to calculate the residual vector,"
        "defined by ``r = g(t,y) - A(t,y) * s``."
        "Used in Lsodi, Lsoibt, Lsodis",
        type=callable,
    ),
    ydoti=dict(
        help="Real array for the initial value of dy/dt.",
        type=(list, tuple, np.ndarray),
        extra_check=lambda float_seq: np.asarray(
            map(lambda x: isinstance(x, float), float_seq)
        ).all(),
        default=[],
    ),
    # ja, ia, jc & ic are used to describe the sparse structure
    # of matrices
    ja=dict(
        help="Integer array containing the row indices of nonzero entries "
        "in a sparse matrix (CRS storage scheme).",
        type=(list, tuple, np.ndarray),
    ),
    ia=dict(
        help="Integer array containing info where the different rows "
        "of a sparse matrix start (CRS storage scheme).",
        type=(list, tuple, np.ndarray),
    ),
    # ml, mu describe banded Jacobian matrix.
    ml=dict(help="Number of lower non-zero diagonals in a banded Jacobian.", type=int),
    mu=dict(help="Number of upper non-zero diagonals in a banded Jacobian.", type=int),
    # mb, nb describe the block-tridiagonal form of matrix.
    # Used in Lsoibt.
    mb=dict(
        help="Block size,  mb>=1, mb*nb = neq (number of equations).",
        type=int,
        extra_check=lambda x: x >= 1,
    ),
    nb=dict(
        help="Number of blocks in the main diagonal. nb>=4.",
        type=int,
        extra_check=lambda x: x >= 4,
    ),
    # Fortran versions of f, jac, g (can be used when solver is in Fortran)
    f_f77=dict(help="Fortran subroutine for f.", type=callable),
    g_f77=dict(help="Fortran subroutine for g.", type=callable),
    jac_f77=dict(help="Fortran subroutine for jac.", type=callable),
    myadvance=dict(
        help="User supplied function to advance current solution"
        " one step forward. See documents of class MySolver.",
        type=callable,
    ),
)


class Solver:
    _required_parameters = ["f"]
    _optional_parameters = [
        "f_args",
        "f_kwargs",
        "complex_valued",
        "disk_storage",
        "verbose",
        "u_exact",
    ]

    def __init__(self, f, **kwargs):
        self._parameters = dict(
            (key, value.copy())
            for key, value in _parameters.items()
            if key in self._optional_parameters or key in self._required_parameters
        )

        self.adjust_parameters()

        for name in self._parameters:
            if "default" in self._parameters[name]:
                setattr(self, name, self._parameters[name]["default"])

        nones = [name for name in kwargs if kwargs[name] is None]
        for name in nones:
            del kwargs[name]

        self.set(**kwargs)

        if f is not None:
            self.users_f = f
            if not callable(f):
                raise ValueError("f must be a callable function.")
            if "f_args" in self._optional_parameters:
                self.f = lambda u, t: np.asarray(f(u, t, *self.f_args, **self.f_kwargs))
            else:
                self.f = lambda u, t: np.asarray(f(u, t))

        self.initialize()

    def name(self):
        return self.__class__.__name__

    def initialize(self):
        return None

    def set(self, strict=False, **kwargs):
        kwargs_copy = kwargs.copy()
        for name in kwargs_copy:
            if name not in self._parameters:
                if strict:
                    raise ValueError(f"Unknown parameter: {name}")
                del kwargs[name]
            elif kwargs[name] is not None:
                del kwargs[name]
                if hasattr(self, name):
                    del self.__dict__[name]

        self.check_input_types(**kwargs)
        self.check_input_range(**kwargs)

        self.check_extra(**kwargs)

        for name in kwargs:
            setattr(self, name, kwargs[name])

        self.check_conditional_parameters()

    def check_input_types(self, **kwargs):
        parameters = self._parameters

        arg_type_list = [
            (name, parameters[name]["type"], kwargs[name])
            for name in parameters
            if name in kwargs and "type" in parameters[name]
        ]
        # name in parameters -->  valid inputs in current class
        # name in kwargs       -->  existing inputs in current instance
        # 'type' in parameters[name]   --> type is specified to be checked

        for name, types, value in arg_type_list:
            # (Ex: types = (callable,int)
            if not isinstance(types, (list, tuple)):
                types = [types]  # make a type list
            ok_type = False
            for tp in types:
                if tp == callable:
                    if callable(value):
                        # value should be a callable object
                        ok_type = True
                else:
                    if isinstance(value, tp):
                        ok_type = True
            if not ok_type:
                # make types more read-friendly in error message if
                # it contains callable as a type
                include_callable = callable in types
                if include_callable:
                    types = [tp for tp in types if tp != callable]
                    types.append("callable function")
                raise TypeError(
                    f"{name}={value} is illegal - type={types}", f"Must be from {types}"
                )
        return True

    def check_input_range(self, **kwargs):
        """Check whether all existing inputs are in right specified range."""

        parameters = self._parameters
        arg_type_list = [
            (name, parameters[name]["range"], kwargs[name])
            for name in parameters
            if name in kwargs and "range" in parameters[name]
        ]
        # name in parameters -->  valid inputs in current class
        # name in kwargs       -->  existing inputs in current instance
        # 'range' in parameters[name]   --> range is specified to be checked

        for name, ranges, value in arg_type_list:
            if isinstance(value, (float, int, complex)):
                # value is a comargble number
                if len(ranges) == 2:  # ranges is an interval
                    low, high = ranges
                    if not ((low <= value <= high) or (low >= value >= high)):
                        raise ValueError(f"{name}={value} is illegal - range={ranges}")
                else:  # range is a list of valid values
                    if value not in ranges:
                        raise ValueError(f"{name}={value} is illegal - range={ranges}")
        return True

    def check_extra(self, **kwargs):
        p = self._parameters
        prm_type_list = [
            (name, p[name]["extra_check"], kwargs[name])
            for name in p
            if name in kwargs and "extra_check" in p[name]
        ]
        # name in parameters -->  valid inputs in current class
        # name in kwargs       -->  existing inputs in current instance
        # 'extra_check' in parameters[name]
        #           --> extra functions is specified to check the value

        for name, check_funcs, value in prm_type_list:
            try:
                if not check_funcs(value):
                    raise ValueError(f"{name}={value} is illegal")
            except Exception as e:
                raise ValueError(f"{name}={value} is illegal - {e}")
        return True

    def get(self, parameter_name=None, print_info=False):
        if parameter_name is None:
            all_args = dict(
                [
                    (name, getattr(self, name, None))
                    for name in self._parameters
                    if hasattr(self, name)
                ]
            )
            del all_args["f"]
            all_args["name of f"] = self.users_f.func_name
            if "jac" in all_args:
                del all_args["jac"]
                all_args["name of jac"] = self.users_jac.func_name

            if print_info:
                print(pprint.pformat(all_args))
            return all_args
        else:
            if hasattr(self, parameter_name):
                value = getattr(self, parameter_name)
                if print_info:
                    print(f"{parameter_name} = {value}")
                return value
            else:
                raise AttributeError(f"{parameter_name} is not set")

    def get_parameter_info(self, print_info=False):
        if print_info:
            print(
                f"Legal parameters for {self.name()}:", f"are {self._parameters.keys()}"
            )
            return None
        return self._parameters

    def set_initial_conditions(self, u0):
        try:
            self.neq = len(u0)
            u0 = np.asarray(u0)
        except TypeError:
            self.neq = 1
            if isinstance(u0, int):
                u0 = float(u0)

        self.u0 = u0

    def solve(self, time_points, terminate=None):
        if terminate is None:

            def terminate(u, t, step_no):
                return False

        self.t = np.asarray(time_points)
        self.n = 0  # time step counter
        self.initialize_for_solve()
        self.validate_data()

        N = self.t.size - 1
        for n in range(N):
            self.n = n
            self.u[n + 1] = self.advance()

            if self.verbose > 2:
                print(f"{self.__class__.__name__}: step {n+1}, t={self.t[n+1]}")
            if terminate(self.u, self.t, n + 1):
                print(f"{self.__class__.__name__}: terminated at t = {self.t[n+1]}")
            if not self.disk_storage:
                self.u, self.t = self.u[: n + 2], self.t[: n + 2]
            break

        if self.disk_storage:
            self.u.flush()

        return self.u, self.t

    def advance(self):
        raise NotImplementedError

    def initialize_for_solve(self):
        if not hasattr(self, "u0"):
            raise AttributeError("Initial conditions not set")

        self._allocate_u(self.t)
        self.u[0] = self.u0

        return None

    def _allocate_u(self, t_array):
        if isinstance(self.disk_storage, bool) and self.disk_storage:
            self.disk_storage = "tmp_lsodes.dat"
        N = t_array.size - 1
        if self.neq == 1:
            if self.disk_storage:
                self.u = np.memmap(
                    self.disk_storage,
                    dtype=float,
                    mode="w+",
                    shape=(N + 1),
                )
            else:
                self.u = np.zeros(N + 1, float)
        else:
            if self.disk_storage:
                self.u = np.memmap(
                    self.disk_storage,
                    dtype=float,
                    mode="w+",
                    shape=(N + 1, self.neq),
                )
            else:
                self.u = np.zeros((N + 1, self.neq), float)

    def validate_data(self):
        if (not isinstance(self.t, (list, tuple, np.ndarray))) or (
            not np.asarray(
                # all items in self.t should be numbers
                [isinstance(t, (int, float)) for t in self.t]
            ).all()
        ):
            raise TypeError(f"Time_points={self.t} is illegal - type={type(self.t)}")

        # self.t should be supplied in an asscending/descending order
        t_sorted = sorted(self.t, reverse=self.t[0] > self.t[-1])
        if list(self.t) != list(t_sorted):
            raise ValueError(f"Time_points={self.t} is illegal - not sorted")

        # Test whether all required parameters are provided
        for arg in self._required_parameters:
            if not hasattr(self, arg):
                raise ValueError(f"Required parameter {arg} is not set")
        return True

    def check_conditional_parameters(self):
        parameters = self._parameters
        # Parameters with condition-list settings
        with_condition_args = [
            (name, parameters[name]["condition-list"], str(getattr(self, name)))
            for name in parameters
            if name in self.__dict__ and "condition-list" in parameters[name]
        ]

        # name in parameters   -->  valid inputs for current class
        # name in self.__dict__ -->  existing inputs for curremt instance
        # 'condition-list' in parameters[name] -->
        #                       'condition-list' is specified to check

        for name, conditions, value in with_condition_args:
            # Ex: name = 'iter_method'
            #     conditions = {'1':(('jac','jac_f77'),), '4':.., '5':..})
            #     value = '1'
            if value in conditions:
                # There is conditional requirements for current value
                condition_args = conditions[value]
                # list/tuple for conditional parameters
                for arg in condition_args:
                    if not isinstance(arg, str):
                        # arg is a list for alternative parameters
                        # e.g. ('jac', 'jac_f77')
                        # Either 'jac' or 'jac_f77' should be supplied
                        found = bool(
                            [p_name for p_name in arg if hasattr(self, p_name)]
                        )
                        arg_print = "One of %s" % str(arg)
                    else:  # arg is a single parameter
                        found = hasattr(self, arg)
                        arg_print = arg
                    if not found:
                        raise ValueError(f"{name}={value} requires {arg_print}")


class ODEPACK(Solver):
    def adjust_parameters(self): ...


class LSODES(ODEPACK):
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

    def adjust_parameters(self):
        # If jac_column is input in form of jac(u,t,j),
        # wrap jac_column to jac_column_f77(t,u,j-1) for Fortran code.
        self._parameters["jac_column"]["paralist_old"] = "u,t,j-1,ia,ja"
        self._parameters["jac_column"]["paralist_new"] = "t,u,j,ia,ja"
        self._parameters["jac_column"]["name_wrapped"] = "jac_column_f77"
        # If jac_column is input in form of jac(t,u,j),
        # wrap it to the general form jac_column(u,t,j) for switch_to().
        self._parameters["jac_column_f77"]["paralist_old"] = "t,u,j+1"
        self._parameters["jac_column_f77"]["paralist_new"] = "u,t,j"
        self._parameters["jac_column_f77"]["name_wrapped"] = "jac_column"

        self._parameters["moss"]["range"] = range(4)
        self._parameters["moss"]["condition-list"] = {
            1: [
                ("jac_column", "jac_column_f77"),
            ],
            0: ["ia", "ja"],
        }

        self._parameters["moss"]["help"] = (
            """\
Method choice to obtain sparse structure with 3 possible
values:

  0. The user has supplied IA, JA.
  1  The user has supplied JAC_COLUMN and the
     sparse structure will be obtained from NEQ initial
     calls to JAC_COLUMN.
  2. The sparse structure will be obtained from
      NEQ+1 initial calls to F.""",
        )

        self._parameters["corrector_iter_method"]["condition-list"] = {
            "1": [
                ("jac_column", "jac_column_f77"),
            ],
        }
        self._parameters["corrector_iter_method"]["help"] = (
            """\
Corrector iteration method choice with 4
possible values:

 0. Functional iteration without any Jacobian
    matrix.
 1. Chord iteration with user-supplied sparse
    Jacobian.
 2. Chord iteration with internally generated
    sparse Jacobian matrix.
 3. Chord iteration with internally generated
    diagonal Jacobian matrix.""",
        )
        ODEPACK.adjust_parameters(self)

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
