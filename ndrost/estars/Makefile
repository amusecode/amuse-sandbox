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

RM ?= rm
JAVAC ?= javac
JAR ?= jar

all: estars_worker

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ estars_worker worker_code.cc
	$(RM) -f *.jar *.java
	$(RM) -rf tmp
	make -C src clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py eStarsInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py eStarsInterface -o $@

CodeInterface.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=interface interface.py eStarsInterface -o $@

Worker.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=class interface.py eStarsInterface -o $@

#estars_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
#	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

src/lib/estars.jar: src
	make -C src

codeinterface.jar: CodeInterface.java
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -g -d tmp CodeInterface.java
	@echo Building jar file
	$(JAR) -cf codeinterface.jar -C tmp .
	$(RM) -rf tmp

worker.jar: Worker.java
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -g -cp "src/external/*:src/lib/*:codeinterface.jar" -d tmp Worker.java
	@echo Building jar file
	$(JAR) -cf worker.jar -C tmp .
	$(RM) -rf tmp

estars_worker: codeinterface.jar src/lib/estars.jar worker.jar
	$(CODE_GENERATOR) -x --type=java --mode=script interface.py eStarsInterface -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
