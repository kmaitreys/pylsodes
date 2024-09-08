# pylsodes
Python wrapper for `DLSODE` and `DLSODES` solvers from Fortran ODEPACK family of differential equation solvers.

## Prerequisites
A Fortran compiler (I built with `gfortran`) and Python < 3.11 with `numpy`. That's about it.

## Installation
Acquire the latest release. Download the zip and extract into a folder of your liking, and then:
```
cd pylsodes
pip install .
```

Tested with 64-bit Linux and Apple silicon.

## Usage
There is no Python frontend right now (might be there in the future), so you have to write your logic to populate work arrays as per your use case. The solver can easily be called by importing directly from the module and using as is used in Fortran.

```python
from pylsodes import dlsodes # or dlsode
```

Signature of the `rhs` and `jac` functions and description of arguments to the solver is detailed in the typing stub file.
