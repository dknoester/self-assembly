#!/bin/bash
trap "exit 255" SIGINT SIGTERM
PROG=${HOME}/bin/self-assembly-ca-2d-fsm

function testfunction {
    ${PROG} -l $1 --analyze ca_dom_1000x
    ${PROG} -l $1 --analyze ca_movie
}

if [[ $# != 1 ]]; then
    echo "Usage: $0 <path>"
    echo "    path: location at which to search for checkpoints."
    exit -1
fi

for i in `find $1 -name "checkpoint-*.xml.gz"`; do
    testfunction $i
done
