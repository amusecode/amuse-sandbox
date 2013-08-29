AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

all: dist/distributed-amuse.jar

RM ?= rm
JAVAC ?= javac
JAR ?= jar

dist/distributed-amuse.jar: $(shell find src -name "*.java")
	@echo Compiling Distributed Code
	$(RM) -r tmp dist
	mkdir tmp
	mkdir dist
	$(JAVAC) -g -classpath "lib/*:lib/ipl/*:lib/jetty/*:lib/jersey/*:lib/octopus/*"  -d tmp $(shell find src -name "*.java")

	@echo Building jar file
	$(JAR) -cf dist/distributed-amuse.jar -C tmp .
	#$(JAR) -uf dist/asterisk.jar images shaders
	$(RM) -r tmp
	@echo done compiling Distributed code

clean:
	$(RM) -r tmp dist


