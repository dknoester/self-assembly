#!/bin/bash
trap "exit 255" SIGINT SIGTERM

function analyze {
    # 1: executable
    # 2: checkpoint
    $1 -l $2 --analyze ca_dom_adapt
}

if [[ $# != 2 ]]; then
    echo "Usage: $0 <executable> <path>"
    echo "    executable: binary to run on checkpoints"
    echo "    path: location at which to search for checkpoints."
    exit -1
fi

for i in `find $2 -name "checkpoint-*.xml.gz" | grep -v tc0 | grep -v tb0`; do
    DIR=`dirname $i`
    CHECKPOINT=`basename $i`
    pushd ${DIR} >/dev/null
    analyze $1 ${CHECKPOINT}
    popd >/dev/null
done
