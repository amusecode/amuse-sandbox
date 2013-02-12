#!/bin/sh
cd ${AMUSE_DIR}/sandbox/ndrost/estars
java -cp "src/log4j.properties:src:src/lib/*:src/external/*:*" amuse.code.Worker $@
