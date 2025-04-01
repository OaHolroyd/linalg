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
make clean # to ensure that non-debug code is recompiled
make debug
make check
```

## Contents
The various functions are grouped into the following categories:
* none yet, I haven't written any!
