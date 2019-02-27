#!/bin/bash
if [ $# != 2 ]
then
    echo "usage: $0 <output_dir> <worker_id>"
    exit 1
fi

OUTPUT_DIR=$1
WORKER_ID=$2



udocker run --rm -v$HOME:$HOME  $HOME/simexp/psom/psom_worker.py -d $OUTPUT_DIR -w $WORKER_ID
