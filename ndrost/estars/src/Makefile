AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

all: lib/estars.jar

RM ?= rm
JAVAC ?= javac
JAR ?= jar

lib/estars.jar: src/*/*/*/*.java
	@echo Compiling Amuse/Ibis lib
	$(RM) -r lib tmp
	mkdir tmp
	mkdir lib
	$(JAVAC) -classpath "external/*:../codeinterface.jar"  -d tmp src/*/*/*/*.java
	#$(JDK)/bin/java -classpath $(CLASSPATH) ibis.compile.Ibisc tmp

	@echo Building jar file
	$(JAR) -cf lib/estars.jar -C tmp .
	#$(JDK)/bin/jar -uf lib/ibis-amuse.jar images/*
	$(RM) -r tmp
	@echo done compiling Stars lib

clean:
	$(RM) -r lib tmp


