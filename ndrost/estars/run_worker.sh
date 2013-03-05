#!/bin/sh
cd ${AMUSE_DIR}/sandbox/ndrost/estars/src
echo java -cp "log4j.properties:.:lib/*:external/*:../*:external/jogl/*" amuse.code.Worker $@
java -cp "log4j.properties:.:lib/*:external/*:../*:external/jogl/*" amuse.code.Worker $@
