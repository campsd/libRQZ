include ../make.inc

DOBJS := d_memorymgmt.o d_ctransformations.o d_setpoles.o d_swappoles12.o \
         d_swappoles22.o d_swappolesmk.o d_computepoles.o d_deflations.o \
				 d_RQZ.o d_aed.o d_RQZtp.o
ZOBJS := z_memorymgmt.o z_ctransformations.o z_setpoles.o z_swappoles.o \
         z_computepoles.o z_deflations.o z_RQZ.o z_aed.o z_RQZtp.o
UOBJS := u_parameters.o u_activeparts.o

all: srcs

srcs: $(UOBJS) $(ZOBJS) $(DOBJS)

# utilities
UDIR := utilities
u_parameters.o:
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(UDIR)/u_parameters.f90 \
	-I../include -o $(OBJDIR)/$@

u_activeparts.o:
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(UDIR)/u_activeparts.f90 \
	-I../include -o $(OBJDIR)/$@

# d (double implementation)
DDIR := d

d_memorymgmt.o: u_parameters.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_memorymgmt.f90 \
	-I../include -o $(OBJDIR)/$@

d_ctransformations.o: u_parameters.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_ctransformations.f90 \
	-I../include -o $(OBJDIR)/$@

d_computepoles.o: u_parameters.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_computepoles.f90 \
	-I../include -o $(OBJDIR)/$@

d_setpoles.o: u_parameters.o d_ctransformations.o d_computepoles.o \
	            d_swappoles12.o d_swappoles22.o d_deflations.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_setpoles.f90 \
	-I../include -o $(OBJDIR)/$@

d_deflations.o: u_activeparts.o u_parameters.o d_ctransformations.o \
	              d_memorymgmt.o d_computepoles.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_deflations.f90 \
	-I../include -o $(OBJDIR)/$@

d_swappoles12.o: u_parameters.o u_activeparts.o d_ctransformations.o \
	               d_memorymgmt.o z_ctransformations.o d_computepoles.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_swappoles12.f90 \
	-I../include -o $(OBJDIR)/$@

d_swappoles22.o: u_parameters.o u_activeparts.o d_ctransformations.o \
	               d_memorymgmt.o z_ctransformations.o d_computepoles.o \
								 d_swappoles12.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_swappoles22.f90 \
	-I../include -o $(OBJDIR)/$@

d_swappolesmk.o: u_parameters.o u_activeparts.o d_ctransformations.o \
	               d_memorymgmt.o z_ctransformations.o d_computepoles.o \
								 d_swappoles12.o d_swappoles22.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_swappolesmk.f90 \
	-I../include -o $(OBJDIR)/$@

d_RQZ.o: u_parameters.o u_activeparts.o d_setpoles.o d_swappoles12.o \
	       d_swappoles22.o d_computepoles.o d_deflations.o d_memorymgmt.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_RQZ.f90 \
	-I../include -o $(OBJDIR)/$@

d_aed.o: u_activeparts.o u_parameters.o d_memorymgmt.o d_swappoles12.o d_RQZ.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_aed.f90 \
	-I../include -o $(OBJDIR)/$@

d_RQZtp.o: u_parameters.o u_activeparts.o d_setpoles.o d_swappoles12.o \
	         d_swappoles22.o d_swappolesmk.o d_computepoles.o d_aed.o \
					 d_deflations.o d_RQZ.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(DDIR)/d_RQZtp.f90 \
	-I../include -o $(OBJDIR)/$@

# z (complex double implementation)
ZDIR := z

z_memorymgmt.o: u_parameters.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_memorymgmt.f90 \
	-I../include -o $(OBJDIR)/$@

z_ctransformations.o: u_parameters.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_ctransformations.f90 \
	-I../include -o $(OBJDIR)/$@

z_setpoles.o: u_parameters.o z_memorymgmt.o z_ctransformations.o z_swappoles.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_setpoles.f90 \
	-I../include -o $(OBJDIR)/$@

z_swappoles.o: u_parameters.o u_activeparts.o z_ctransformations.o \
	             z_memorymgmt.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_swappoles.f90 \
	-I../include -o $(OBJDIR)/$@

z_computepoles.o: u_parameters.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_computepoles.f90 \
	-I../include -o $(OBJDIR)/$@

z_deflations.o: u_activeparts.o u_parameters.o z_ctransformations.o \
	              z_swappoles.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_deflations.f90 \
	-I../include -o $(OBJDIR)/$@

z_RQZ.o: u_parameters.o u_activeparts.o z_setpoles.o z_swappoles.o \
	       z_computepoles.o z_deflations.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_RQZ.f90 \
	-I../include -o $(OBJDIR)/$@

z_aed.o: u_activeparts.o u_parameters.o z_memorymgmt.o z_swappoles.o z_RQZ.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_aed.f90 \
	-I../include -o $(OBJDIR)/$@

z_RQZtp.o: u_parameters.o u_activeparts.o z_setpoles.o z_swappoles.o \
	         z_computepoles.o z_aed.o z_deflations.o z_RQZ.o
	$(FC) $(FFLAGS) -J$(INCDIR) -cpp -fPIC -c $(ZDIR)/z_RQZtp.f90 \
	-I../include -o $(OBJDIR)/$@

clean:
	@rm -f $(OBJDIR)/*.o
