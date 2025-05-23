# LinAlg - Linear algebra functions and examples

The aim of this repository is not to provide an exhaustive set of
hyper-optimised linear algebra functions, but to act as a reference for commonly
used algorithms (defined by whether I need to use them regularly) with
well-documented implementations. Often functions will be designed to be copied
and adapted to a specific use case rather than be a perfect generalised
implementation. For example, as written, the block-solver functions don't
provide any benefit over calling LAPACK's `dgesv`, but are meant as a framework
to which other sparse solvers can be inserted.


## Testing
Although this is meant as a reference to be chopped, changed, and copied
into other projects, you can compile the code and run tests for optimised
versions:

```bash
make
make check
```

and debug versions:

```bash
make DEBUG=1 rebuild # to ensure that non-debug code is recompiled
make check
```

Note that the Makefile is currently set up to use Homebrew's `clang`, so you
might need to change the compiler to work on your machine, e.g.

```bash
make CC=gcc
```

## Notes

* Whenever a function takes a matrix as input, it is assumed to be in
  flattened, **row-major** order
* All functions are written using `double`s

## Contents
The various functions are grouped into the following categories:

### Helpers

* [Matrix memory management](/src/alloc.h)
* [Matrix IO](/src/io.h)

### Solvers

* [General LU solvers](/src/lu_solve.h)
* [Block-decomposed solvers](/src/block_solve.h)
* [Pentadiagonal solvers](/src/pent_solve.h)
