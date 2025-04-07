# CafeMol

CafeMol is a general-purpose coarse-grained (CG) biomolecular modeling and simulation software.
It can simulate proteins, nucleic acids, and their mixture with various CG models.

## How to use

1. Edit your input file (e.g. cafemol\_go\_1chain.inp)
1. Run cafemol with this command: `./bin/cafemol <path to your input file>`

## How to build

```sh
$ mkdir build && cd build
$ cmake -DCMAKE_Fortran_COMPILER=<your favorite compiler> [-DCMAKE_BUILD_TYPE=<Release or Debug>] ..
$ make -j8
```

- To compile the code with OpenMP support, add '-DOMP\_PAR=ON' to the cmake command.

- To compile the code with a different precision for floats, use '-DPREC=[REAL32, REAL64, REAL128]'

- To test the compiled binary, run 'make test' inside the build directory.

## References

* Kenzaki, H., Koga, N., Hori, N., Kanada, R., Li, WF., Okazaki, K., Yao, XQ., Takada, S.
CafeMol: A Coarse-Grained Biomolecular Simulator for Simulating Proteins at Work
Journal of Chemical Theory and Computation (2011) 7(6), pp 1979-1989
DOI: [10.1021/ct2001045](http://dx.doi.org/10.1021/ct2001045)
# Cafemol_jigakuyou
