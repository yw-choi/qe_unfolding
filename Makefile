# Makefile for postprocessing (PP)

include make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)$(TOPDIR)/PW/src \
          $(MOD_FLAG)$(TOPDIR)/dft-d3/

OBJS = unfold.o

PWOBJS = $(TOPDIR)/PW/src/libpw.a
QEOBJS = $(TOPDIR)/Modules/libqemod.a $(TOPDIR)/KS_Solvers/libks_solvers.a \
         $(TOPDIR)/FFTXlib/libqefft.a $(TOPDIR)/LAXlib/libqela.a $(TOPDIR)/UtilXlib/libutil.a \
         $(TOPDIR)/dft-d3/libdftd3qe.a $(TOPDIR)/PP/src/libpp.a

MODULES = $(PWOBJS) $(QEOBJS)

all : unfold.x

unfold.x : $(OBJS) $(MODULES) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ \
		$(OBJS) $(MODULES) $(LIBOBJS)  $(QELIBS)

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L

