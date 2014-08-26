# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../..
-include ${AMUSE_DIR}/config.mk

export MESASDK_ROOT:=$(shell pwd)/mesasdk
export PATH:=${MESASDK_ROOT}/bin:${PATH}
export MANPATH:=${MESASDK_ROOT}/share/man:${MANPATH}
export PGPLOT_DIR:=${MESASDK_ROOT}/pgplot

MPIF90 ?= mpif90
FORTRAN ?= gfortran
PYTHON ?= python
FC      = $(FORTRAN)
MPIFC+= $(FCFLAGS)

INCLUDE_DIR=./src/star/make
MESA_DIR = ./src
MESA_LIB=$(MESA_DIR)/lib/libstar.a

include makefile_header

CODE_GENERATOR = $(AMUSE_DIR)/build.py

FS_LIBDIR = ${AMUSE_DIR}/lib/forsockets
LOCAL_LIBDIR = ./lib
LOCAL_FS_LIBDIR = ./lib/forsockets
FS_FLAGS_MESA=-I$(LOCAL_FS_LIBDIR)
FS_LIBS_MESA=-L$(LOCAL_FS_LIBDIR) -lforsockets -lforsocketsf
FS_LIB=$(LOCAL_FS_LIBDIR)/libforsocketsf.a $(LOCAL_FS_LIBDIR)/libforsockets.a

all:	mesa_worker create_zams

$(MESA_DIR)/install:
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@echo "Regardless of whether DOWNLOAD_CODES is set or not, MESA2 will not"
		@echo "be downloaded and build automatically. If you do want MESA2, "
		@echo "download mesa-r6794.zip from"
		@echo "http://sourceforge.net/projects/mesa/files/releases/"
		@echo "and unpack it to the src directory. Also download the mesasdk from"
		@echo "http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk"
		@echo "and unpack it to the mesasdk directory."
		@echo "Note: MESA requires about 10 GB of disk space."
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@make -s --no-print-directory -C . raise_error

clean:
	$(RM) -f *.pyc *.mod *.o worker_code.f90 $(MESA_LIB)
	$(RM) -f mesa_worker
	$(RM) -rf $(LOCAL_FS_LIBDIR)
	$(RM) -rf $(LOCAL_LIBDIR)
	$(RM) -f config.mk

veryclean: clean
	cd $(MESA_DIR); ./clean

distclean: clean
	$(RM) -Rf src
	$(RM) -Rf src.*
	$(RM) -Rf .pc

mesa_worker:	worker_code.f90 interface.o $(FS_LIB)
	$(MPIF90) $(FCFLAGS) $(FFLAGS) $(FCopenmp) $(FS_FLAGS_MESA) $^ -o $@ $(LOAD_MESA_STAR) $(FS_LIBS_MESA)  $(LIBS)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MESAInterface -o $@

$(MESA_LIB):  $(MESA_DIR)/install
	cp -f makefile_header src/utils/makefile_header
	-cd $(MESA_DIR); ./install PATH=$(PATH) MPIFC="$(MPIFC) $(FCFLAGS)"

interface.o: interface.f $(MESA_LIB)
	$(WORK_COMPILE) -I$(INCLUDE_DIR) -o $@ $<

$(FS_LIB): $(FS_LIBDIR)/Makefile $(FS_LIBDIR)/forsocketsf.f90 $(FS_LIBDIR)/forsockets.c $(FS_LIBDIR)/forsockets.h
	-mkdir $(LOCAL_LIBDIR)
	-mkdir $(LOCAL_FS_LIBDIR)
	cp $^ $(LOCAL_FS_LIBDIR)
	cp ${AMUSE_DIR}/config.mk ./
	make -C $(LOCAL_FS_LIBDIR)

create_zams:
	-cd $(MESA_DIR)/star/test_suite/create_zams; ./mk PATH=$(PATH) MPIFC="$(MPIFC) $(FCFLAGS)"
