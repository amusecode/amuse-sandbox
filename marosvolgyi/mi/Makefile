MPICXX   ?= mpicxx

CFLAGS   += -Wall -g 
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

SCLIBS ?= -L../../../lib/stopcond -lstopcond
SCINC ?= -I../../../lib/stopcond

CODELIB = src/libmi.a

AMUSE_DIR?=../../../

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: mi_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ mi_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py MiInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py MiInterface -o $@

mi_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) -I./src  $< $(OBJS) $(CODELIB) -o $@ -L./src -lmi $(SCLIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(SCINC) -I./src -c -o $@ $< 
