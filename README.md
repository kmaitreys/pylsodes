# pylsodes
Python wrapper for `DLSODES` solver from Fortran ODEPACK family of differential equation solvers

## Installation
Acquire the latest release. Download the zip and extract into a folder of your liking, and then:
```
cd pylsodes
pip install .
```

Only 64-bit Linux support for now.

## Usage
There is no Python frontend right now (might be there in the future), so you have to write your logic to populate work arrays as per your use case. The solver can easily be called by importing directly from the module and using as is used in Fortran.

```
from pylsodes import dlsodes
```

Signature of the `rhs` and `jac` functions and description of arguments to the solver is detailed in the typing stub file.
