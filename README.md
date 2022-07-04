# Peacemaker setup

Peacemaker program for Quantum Cluster Equilibrium calculations
## Compiling Peacemaker
Peacemaker is a modern FORTRAN code and thus requires a modern FORTRAN compiler. We recommend a recent version of gfortran which is used for active development. Peacemaker can be built by running

```$ make release```

which should produce a run time optimized binary called `peacemaker`. In case of errors, adjust the makefile to your compiler. We recommend the following compiler flags or your compiler’s equivalents:

* `-O3 highest` optimization level that guarantees standard compliance
* `-fopenmp OpenMP` parallelization
* `-flto` link-time optimization

A version suitable for development and debugging can be built by running

```$ make debug```

Note: Older versions of gfortran are subject to a bug which prevents OpenMP parallelization. If you receive the error message “Attempting to allocate already allocated variable ‘ib’ ”, compile without OpenMP support, or upgrade to a newer compiler version.

## Running Peacemaker

Peacemaker is run by

```$ peacemaker [input] [clusterset]```

where `[input]` is the location of the input file and `[clusterset]` is the location of the clusterset file. The structure of both files is explained in Section 4 of the [manual](manual/manual.pdf). If Peacemaker was compiled with OpenMP parallelization, it can be run in parallel by

```$ OMP_NUM_THREADS=[N] peacemaker [input] [clusterset]```

In this case, `[N]` specifies the number of threads to run with.
