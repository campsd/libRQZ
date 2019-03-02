# libRQZ
# Description:
#   Set compiler and install directories
# ___________________________________________________________________

# Project directory
INSTALLDIR := *FILL OUT INSTALLATION DIRECTORY*

# fortran compiler and flags
FC := gfortran
FFLAGS := -O3

# libraries required by compiler (LAPACK AND BLAS)
LIBS := -llapack -lblas

# Dependent directories
SRCDIR := $(INSTALLDIR)/src
OBJDIR := $(INSTALLDIR)/obj
INCDIR := $(INSTALLDIR)/include
BINDIR := $(INSTALLDIR)/bin
