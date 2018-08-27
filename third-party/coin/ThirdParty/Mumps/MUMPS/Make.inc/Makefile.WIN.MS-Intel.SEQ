#
#  This file is part of MUMPS 4.10.0, built on Tue May 10 12:56:32 UTC 2011
#
# We are grateful to Evgenii Rudnyi for his help and suggestions
# regarding Windows installation.
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. If you would like to
# use other orderings, you need to obtain the corresponding package and modify
# the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS, PATHMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis/ParMetis and SCOTCH/PT-SCOTCH (ver 5.1 and later) orderings are now available for MUMPS.
#

#SCOTCHDIR  = ${HOME}/scotch_5.1_esmumps
#PATHSCOTCH = -LIBPATH:$(SCOTCHDIR)/lib
#ISCOTCH    = -I$(SCOTCHDIR)/include 
# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dscotch in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dptscotch in the ORDERINGSF variable below)

#LSCOTCH    = libscotch.lib libscotcherr.lib
#LSCOTCH    = libptscotch.lib libptscotcherr.lib


LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
PATHPORD = -LIBPATH:$(LPORDDIR)
LPORD    = libpord.lib

#LMETISDIR = /local/metis/
#PATHMETIS = -LIBPATH:$(PATHMETIS)
#IMETIS    = # Should be provided if you use parmetis

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dmetis in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dparmetis in the ORDERINGSF variable below)

#LMETIS    = libmetis.lib
#LMETIS    = libparmetis.lib libmetis.lib

# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
#ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis
ORDERINGSF  = -Dpord
ORDERINGSC  = $(ORDERINGSF)

PATHORDERINGS = $(PATHMETIS) $(PATHPORD) $(PATHSCOTCH)
LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)

# FOr Windows
#End orderings
########################################################################
########################################################################
################################################################################

PLAT    =
LIBEXT  = .lib
OUTC    = -Fo
OUTF    = -Fo
RM = /bin/rm -f
CC = cl
FC = ifort
FL = ifort
AR = lib -out:
#RANLIB = ranlib
RANLIB  = echo
INCSEQ = -I$(topdir)/libseq
LIBSEQ = $(topdir)/libseq/libmpiseq.lib
#LIBBLAS = -L/usr/lib/xmm/ -lf77blas -latlas
#LIBBLAS = -L/local/BLAS -lblas
LIBBLAS = mkl_intel_c.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib
#LIBOTHERS = -lpthread
LIBOTHERS = -link $(PATHORDERINGS)
#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O -MD -Dintel_ -DALLOW_NON_INIT -fpp
OPTL    = 
OPTC    = -O2 -MD
#End Optimized options
INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
