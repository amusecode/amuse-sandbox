#!/bin/bash
cd ${AMUSE_DIR}/sandbox/ndrost/estars/src
source ${AMUSE_DIR}/config.mk
echo ${JAVA} -cp "log4j.properties:.:lib/*:external/*:../*:external/jogl/*" amuse.code.Worker $@
${JAVA} -cp "log4j.properties:.:lib/*:external/*:../*:external/jogl/*" amuse.code.Worker $@
