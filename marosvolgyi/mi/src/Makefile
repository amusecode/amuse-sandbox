CC = g++
CCC = mpicxx
PROF = -pg -Wall -O3 -ffast-math -funroll-loops
CFLAGS = -Wall -O3 -ffast-math -funroll-loops -fomit-frame-pointer $(MUSE_INCLUDE_DIR)

#SAPPATH = /disks/botlek1/iwasawa/work/code/others/sapporo2
#INC = -I/disks/botlek1/iwasawa/work/code/others/sapporo2/libheaders/

OBJS = mpi_interface.o external_field.o merge.o

#all:Nbody_emu.out Nbody_gpu.out Nbody_gpu_pec2.out Nbody_emu_pec2.out

AR = ar ruv
RANLIB = ranlib
RM = rm

all: Nbody_emu_amuse.o Nbody_emu_pec2.out libmi.a

libmi.a: Nbody_emu_amuse.o 
	$(RM) -f $@
	$(AR) $@ Nbody_emu_amuse.o 
	$(RANLIB) $@

Nbody_emu_amuse.o:Nbody_amuse.cc  evolve_emu.o merge.o external_field.o mpi_interface.o sapporo2_dummy.o
	$(CCC) $(CFLAGS) -DEMU -c -o $@ $< $(OBJS) evolve_emu.o sapporo2_dummy.o 

Nbody_emu.out:Nbody.cc  evolve_emu.o merge.o external_field.o mpi_interface.o sapporo2_dummy.o
	$(CCC) $(CFLAGS) -DEMU -o $@ $< $(OBJS) evolve_emu.o sapporo2_dummy.o 

evolve_emu.o:evolve.cc evolve.h
	$(CCC) -DEMU -c $(CFLAGS)  -o $@ $<

Nbody_emu_pec2.out:Nbody.cc  evolve_emu_pec2.o merge.o external_field.o mpi_interface.o sapporo2_dummy.o
	$(CCC) $(CFLAGS) -DPECEC -DEMU -o $@ $< $(OBJS) evolve_emu_pec2.o sapporo2_dummy.o 

evolve_emu_pec2.o:evolve.cc evolve.h
	$(CCC) -DPECEC -DEMU -c $(CFLAGS)  -o $@ $<

#Nbody_gpu_pec2.out:Nbody.cc  evolve_gpu_pec2.o merge.o external_field.o mpi_interface.o 
#	$(CCC) $(CFLAGS) -DPECEC -DSAP -o $@ $< $(OBJS) evolve_gpu_pec2.o -L $(SAPPATH) -lsapporo -lcuda -fopenmp

#evolve_gpu_pec2.o:evolve.cc evolve.h
#	$(CCC) -DPECEC -DSAP -c $(CFLAGS)  -o $@ $<

#Nbody_gpu.out:Nbody.cc  evolve_gpu.o merge.o external_field.o mpi_interface.o 
#	$(CCC) $(CFLAGS) -DSAP -o $@ $< $(OBJS) evolve_gpu.o -L $(SAPPATH) -lsapporo -lcuda -fopenmp

#evolve_gpu.o:evolve.cc evolve.h
#	$(CCC) -DSAP -c $(CFLAGS)  -o $@ $<

merge.o:merge.cc merge.h 
	$(CCC) -c $(CFLAGS)  -o $@ $<

external_field.o:external_field.cc external_field.h
	$(CCC) -c $(CFLAGS)  -o $@ $<

mpi_interface.o:mpi_interface.cc mpi_interface.h
	$(CCC) -c $(CFLAGS)  -o $@ $<

sapporo2_dummy.o:sapporo2_dummy.cc
	$(CCC) -c $(CFLAGS)  -o $@ $<

clean:
	rm *.o *.s Nbody_emu.out Nbody_gpu.out Nbody_emu_pec2.out Nbody_gpu_pec2.out 



