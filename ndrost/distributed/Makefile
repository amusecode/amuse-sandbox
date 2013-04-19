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

all: java_check distributed_worker distributed-server.jar

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

CodeInterface.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=interface interface.py DistributedAmuseInterface -o $@

Worker.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=class interface.py DistributedAmuseInterface -o $@

distributed-server.jar: $(shell find lib -name "*.jar") distributed.jar worker.jar

	$(RM) -r tmp
	mkdir tmp
	unzip -q -o -d tmp "lib/*.jar"
	unzip -q -o -d tmp "lib/ipl/*.jar"
	unzip -q -o -d tmp distributed.jar
	unzip -q -o -d tmp worker.jar
	jar -cf distributed-server.jar -C tmp/ .

distributed.jar: $(shell find src -name "*.java")
	@echo Compiling Distributed Code
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -g -cp "lib/*:lib/ipl/*" -d tmp $(shell find src -name "*.java")
	@echo Building jar file
	$(JAR) -cf distributed.jar -C tmp .
	$(RM) -rf tmp

worker.jar: Worker.java CodeInterface.java Code.java distributed.jar
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -g -cp "lib/*:lib/ipl/*:distributed.jar" -d tmp Worker.java CodeInterface.java Code.java
	@echo Building jar file
	$(JAR) -cf worker.jar -C tmp .
	$(RM) -rf tmp

distributed_worker: worker.jar
	$(CODE_GENERATOR) -x --type=java --mode=script interface.py DistributedAmuseInterface -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
