import ctypes
import sys
from ctypes import pythonapi

sys.path.append("/Users/kmaitreys/Documents/software/pylsodes/build")

import numba as nb
import numpy as np
from _dlsodes import dlsodes
from numba import types
from numba.core import cgutils
from numba.extending import intrinsic


# Function to extract the function pointer
def capsule_name(capsule):
    pythonapi.PyCapsule_GetName.restype = ctypes.c_char_p
    pythonapi.PyCapsule_GetName.argtypes = [ctypes.py_object]
    return pythonapi.PyCapsule_GetName(capsule)


def get_f2py_function_address(capsule):
    name = capsule_name(capsule)
    pythonapi.PyCapsule_GetPointer.restype = ctypes.c_void_p
    pythonapi.PyCapsule_GetPointer.argtypes = [ctypes.py_object, ctypes.c_char_p]
    return pythonapi.PyCapsule_GetPointer(capsule, name)


# Access function pointer for `dlsodes`
func_ptr = get_f2py_function_address(dlsodes._cpointer)


# Signature of the Fortran function
"""
```python
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
) -> Tuple[ArrayLike, float, int]: ...
"""

# Define its signature
func_type = ctypes.CFUNCTYPE(
    # Result type :: Returns updated `y`, `t`, and `istate`
    ctypes.c_void_p,
    # Argument types
    # f: Callable -> must be a function pointer
    ctypes.c_void_p,
    # y: ArrayLike[np.float64] -> must be a pointer to a float64 array
    ctypes.POINTER(ctypes.c_double),
    # t: float -> must be a float
    ctypes.c_double,
    # tout: float -> must be a float
    ctypes.c_double,
    # rtol: float | ArrayLike[np.float64] -> must be a pointer to a float64 array
    ctypes.POINTER(ctypes.c_double),
    # atol: ArrayLike[np.float64] -> must be a pointer to a float64 array
    ctypes.POINTER(ctypes.c_double),
    # itask: int -> must be an integer
    ctypes.c_int,
    # istate: int -> must be an integer
    ctypes.c_int,
    # rwork: ArrayLike[np.float64] -> must be a pointer to a float64 array
    ctypes.POINTER(ctypes.c_double),
    # iwork: ArrayLike[np.int64] -> must be a pointer to a int64 array
    ctypes.POINTER(ctypes.c_int),
    # jac: Callable -> must be a function pointer
    ctypes.c_void_p,
    # mf: int -> must be an integer
    ctypes.c_int,
    # f_extra_args: Tuple = () -> not used in the Fortran function
    ctypes.c_void_p,
    # overwrite_y: int = 0 -> not used in the Fortran function
    ctypes.c_int,
    # jac_extra_args: Tuple = () -> not used in the Fortran function
    ctypes.c_void_p,
)
dlsodes_ctypes = func_type(func_ptr)


# Intrinsics for scalar-to-pointer and pointer-to-scalar
@intrinsic
def val_to_ptr(typingctx, data):
    def codegen(context, builder, signature, args):
        ptr = cgutils.alloca_once_value(builder, args[0])
        return ptr

    sig = types.CPointer(data)(data)
    return sig, codegen


@intrinsic
def ptr_to_val(typingctx, ptr):
    def codegen(context, builder, signature, args):
        val = builder.load(args[0])
        return val

    sig = ptr.dtype(ptr)
    return sig, codegen


# Wrap the Fortran function in Numba
@nb.njit(
    "void(float64(float64, float64[:]), float64[:], float64[:], float64[:], float64[:], float64[:], int64[:], int64[:], float64[:], int64[:], float64(float64, float64[:]), int64[:])"
)
def wrapped_dlsodes(f, y, t, tout, rtol, atol, itask, istate, rwork, iwork, jac, mf):
    f_ptr = val_to_ptr(f)
    y_ptr = val_to_ptr(y)
    rtol_ptr = val_to_ptr(rtol)
    atol_ptr = val_to_ptr(atol)
    rwork_ptr = val_to_ptr(rwork)
    iwork_ptr = val_to_ptr(iwork)
    jac_ptr = val_to_ptr(jac)
    t_ptr = val_to_ptr(t)
    tout_ptr = val_to_ptr(tout)
    itask_ptr = val_to_ptr(itask)
    istate_ptr = val_to_ptr(istate)
    mf_ptr = val_to_ptr(mf)

    dlsodes_ctypes(
        f_ptr,
        y_ptr,
        t_ptr,
        tout_ptr,
        rtol_ptr,
        atol_ptr,
        itask_ptr,
        istate_ptr,
        rwork_ptr,
        iwork_ptr,
        jac_ptr,
        mf_ptr,
        ctypes.c_void_p(),
        0,
        ctypes.c_void_p(),
    )

    # Extract the updated values
    y[:] = ptr_to_val(y_ptr)
    t[:] = ptr_to_val(t_ptr)
    istate[:] = ptr_to_val(istate_ptr)


# Test the wrapped function
def f(t, y):
    return -y


def jac(t, y):
    return -1


y = np.array([1.0])
t = np.array([0.0])
tout = np.array([1.0])
rtol = np.array([1e-8])
atol = np.array([1e-8])
itask = np.array([1])
istate = np.array([1])
rwork = np.zeros(100)
iwork = np.zeros(100)
mf = np.array([22])

wrapped_dlsodes(f, y, t, tout, rtol, atol, itask, istate, rwork, iwork, jac, mf)

print(y, t, istate)
