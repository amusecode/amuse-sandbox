# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../..
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: estars_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ estars_worker worker_code.cc
	make -C src clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py eStarsInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py eStarsInterface -o $@

#estars_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
#	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

src/lib/estars.jar:
	make -C src

codeinterface.jar: amuse/code/CodeInterface.java
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -d tmp amuse/code/CodeInterface.java
	@echo Building jar file
	$(JAR) -cf codeinterface.jar -C tmp .
	$(RM) -r tmp

worker.jar: amuse/code/Worker.java
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -cp "src/external/*:src/lib/*:codeinterface.jar" -d tmp amuse/code/Worker.java
	@echo Building jar file
	$(JAR) -cf worker.jar -C tmp .
	$(RM) -r tmp

estars_worker: worker_code.h codeinterface.jar src/lib/estars.jar worker_code.cc worker.jar
	cp run_worker.sh estars_worker	

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
