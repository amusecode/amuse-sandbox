AMUSE_DIR ?= ../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

CXXFLAGS = -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
LDFLAGS = -lm $(MUSE_LD_FLAGS)

CODEPATH = src
$(warning make does: : $(MAKE))

OBJS = src/PACO.o
CXX = g++

all: compile worker_code

compile: PACO

PACO: 
	$(MAKE) -C $(CODEPATH) PACO.o

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py PACOInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py PACOInterface -o $@

worker_code: worker_code.cc $(OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(LIBS) $(CODELIB) $(OBJS) -o $@

cleanall: clean
	$(RM) muse_worker *~

clean:
	rm -f *.so *.o *.pyc muse_worker.cc 
	$(MAKE) -C $(CODEPATH) clean	

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.c.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
