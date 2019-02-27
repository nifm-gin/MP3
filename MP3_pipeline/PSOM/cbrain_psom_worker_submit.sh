#!/bin/bash

if [ $# != 2 ]
then
    echo "usage: $0 <output_dir> <worker_id>"
    exit 1
fi

OUTPUT_DIR=$1
WORKER_ID=$2

FILE=`mktemp .new-task-XXXX.json`

cat << NEWTASK > ${FILE}
{
  "tool-class": "CbrainTask::PSOMWorker",
  "description": "A PSOM worker submitted by PSOM through cbrain-psom-worker-submit.",
  "parameters": {
          "output_dir": "${OUTPUT_DIR}",
          "worker_id": "${WORKER_ID}"
      }
}

NEWTASK
