############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          SLmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   March 20, 1995
#
#  Modified:
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################
#
ERNST		=	/people/ehaunsch
#
SHELL         = /bin/sh
#
#  The complete path to the top level of ScaLAPACK directory, usually
#  $(HOME)/SCALAPACK
#
home          = $(ERNST)/SCALAPACK
#
#  The platform identifier to suffix to the end of library names
#
PLAT          = IRIX64
#
LAPACKdir     = $(ERNST)/LAPACK
#
#  BLACS setup.  All version need the debug level (0 or 1),
#  and the directory where the BLACS libraries are
#
BLACSDBGLVL   = 0
BLACSdir      = $(ERNST)/BLACS/LIB
#
#  MPI setup; uncomment and tailor to your system if using MPIBLACS
#  Will need to comment out the default native BLACS setup below below
#
USEMPI        = -DUsingMpiBlacs
SMPLIB        = /src/icl/MPI/mpich1.1/lib/IRIX64/ch_p4/libmpi.a
SMPLIB        = -lmpi
BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLAT)-$(BLACSDBGLVL).a
BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a
BLACSFINIT    = $(home)/libmpiblacs64.a
BLACSCINIT    = $(home)/libmpiblacs64.a
BLACSLIB      = $(home)/libmpiblacs64.a
TESTINGdir    = $(home)/TESTING

#
#  PVMBLACS setup, comment out next 6 lines if using MPI
#
#USEMPI        =
#SMPLIB        = $(PVM_ROOT)/lib/$(PLAT)/libpvm3.a
#BLACSFINIT    =
#BLACSCINIT    =
#BLACSLIB      = $(BLACSdir)/blacs_PVM-$(PLAT)-$(BLACSDBGLVL).a
#TESTINGdir    = $(HOME)/pvm3/bin/$(PLAT)

CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

BLASLIB       = -lblas
#
#  The directories to find the various pieces of ScaLapack
#
PBLASdir      = $(home)/PBLAS
SRCdir        = $(home)/SRC
TESTdir       = $(home)/TESTING
PBLASTSTdir   = $(TESTINGdir)
TOOLSdir      = $(home)/TOOLS
REDISTdir     = $(home)/REDIST
REDISTTSTdir  = $(TESTINGdir)
#
#  The fortran and C compilers, loaders, and their flags
#
F77           = f77
CC            = cc
NOOPT         = -64 -mips4 -r10000
#NOOPT         = -n32 -mips4 -r10000
F77FLAGS      = -O2 $(NOOPT)
CCFLAGS       = $(F77FLAGS)
SRCFLAG       =
F77LOADER     = $(F77)
CCLOADER      = $(CC)
F77LOADFLAGS  = $(F77FLAGS)
CCLOADFLAGS   = $(CCFLAGS)
#
#  C preprocessor defs for compilation 
#  (-DNoChange, -DAdd_, -DUpCase, or -Df77IsF2C)
#
CDEFS         = -DAdd_ $(USEMPI)
#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#
ARCH          = ar
ARCHFLAGS     = cr
RANLIB        = echo
#
#  The name of the libraries to be created/linked to
#
PBLASLIB      = $(home)/pblas_$(PLAT).a
SCALAPACKLIB  = $(home)/scalapack_$(PLAT).a
TOOLSLIB      = $(home)/tools_$(PLAT).a
REDISTLIB     = $(home)/redist_$(PLAT).a
#
PBLIBS        = $(PBLASLIB) $(TOOLSLIB) $(FBLACSLIB) $(BLASLIB) $(SMPLIB)
PRLIBS        = $(REDISTLIB) $(TOOLSLIB) $(CBLACSLIB) $(SMPLIB)
RLIBS         = $(SCALAPACKLIB) $(PBLASLIB) $(REDISTLIB) $(TOOLSLIB) \
                $(FBLACSLIB) $(CBLACSLIB) $(BLASLIB) $(SMPLIB)
LIBS          = $(SCALAPACKLIB) $(PBLIBS)
