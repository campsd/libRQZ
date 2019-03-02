# libRQZ
Fortran 90 library implementing the multishift,
multipole rational QZ method with aggressive early deflation.

## Preprints

This implementation is based on the articles:

* Daan Camps, Karl Meerbergen, and Raf Vandebril. [A multishift, multipole rational QZ method with aggressive early deflation](https://arxiv.org/abs/1902.10954). eprint arXiv:1902.10954.
* Daan Camps, Karl Meerbergen, and Raf Vandebril. [A rational QZ method](https://arxiv.org/abs/1802.04094). eprint arXiv:1802.04094.

## Webpage

More online resources are available on our [webpage](http://numa.cs.kuleuven.be/software/rqz/).

## Installation and usage

We have tested the code with **gfortran** on Ubuntu, other compilers and operating systems have *not* been tested.

An example makefile is included. You need to fill out the installation directory in `make.inc`. In that same file the compiler and compilation flags can be set and the references to the LAPACK and BLAS libraries can be provided.

If everything is set up correctly, the library can be installed via:

```
make install
```

*Note:* it is possible that the installation fails the first time, but works the second time once all object files have been built in the correct directories.

Deinstallation is done via:

```
make uninstall
```
### Usage

The main computational routine for real-double problems is `d_RQZm` located in the module `dRQZtp`;
and for complex-double problems it is `z_RQZm` located in the module `zRQZtp`.

## Version and open issues

The  current version of the library is v0.1.0.

We are aware of some open issues with the software and working on solving them in a future update.
If you happen to find a bug, feel free to raise it on the [issues page](https://github.com/campsd/libRQZ/issues)

