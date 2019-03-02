# libRQZ
# Description:
#   Master make file
# ___________________________________________________________________

include make.inc .master.inc

SRCS := $(wildcard ./src/*/*.f90)
OBJS := $(wildcard obj/*.o)

all: lib$(LIBNAME).$(SLIB).$(VERSION)

install: lib$(LIBNAME).$(SLIB).$(VERSION)
	@mkdir -p $(INSTALLDIR)/$(LIBNAME)/lib
	@cp ./lib$(LIBNAME).$(SLIB).$(VERSION) $(INSTALLDIR)/$(LIBNAME)/lib

lib$(LIBNAME).$(SLIB).$(VERSION): srcs
	$(FC) $(FFLAGS) -shared -o lib$(LIBNAME).$(SLIB).$(VERSION) $(OBJS)

srcs:
	@$(MAKE) $@ -C ./src

uninstall: clean
	@rm -f $(INSTALLDIR)/$(LIBNAME)/lib/lib$(LIBNAME).$(SLIB).$(VERSION)
	@rm -f $(INCDIR)/*.mod

clean:
	@$(MAKE) clean -C ./src
	@rm -f lib$(LIBNAME).$(SLIB).$(VERSION)
