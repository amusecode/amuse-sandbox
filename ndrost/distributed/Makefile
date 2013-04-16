# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../..
-include $(AMUSE_DIR)/config.mk

CODE_GENERATOR = $(AMUSE_DIR)/build.py

RM ?= rm
JAVAC ?= javac
JAR ?= jar
MV ?= mv
CP ?= cp

all: java_check distributed_worker

java_check:
#Makefile requires these checks to be at the beginning of the line
ifeq ($(JAVA_ENABLED),no)
        $(error "Java support not configured in AMUSE")
endif
ifeq (,$(findstring $(JAVA_VERSION), 1.7))
        $(error "Java 1.7 required for the Distributed code, but AMUSE configured with Java $(JAVA_VERSION)")
endif

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ distributed_worker worker_code.cc
	$(RM) -f *.jar CodeInterface.java Worker.java
	$(RM) -rf tmp
	make -C src clean

CodeInterface.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=interface interface.py DistributedAmuseInterface -o $@

Worker.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=class interface.py DistributedAmuseInterface -o $@

#src/lib/distributed.jar: $(shell find src)
#	@echo Compiling Asterisk Lib
#	make -C src

worker.jar: Worker.java CodeInterface.java Code.java 
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -g -cp "src/lib/*:src/lib/jogl/*:src/dist/*" -d tmp *.java
	@echo Building jar file
	$(JAR) -cf worker.jar -C tmp .
	$(RM) -rf tmp

distributed_worker: worker.jar
	$(CODE_GENERATOR) -x --type=java --mode=script interface.py DistributedAmuseInterface -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
