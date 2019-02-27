#!/usr/bin/env bash
# Script config
# docker command
DOCKER_EXE=docker
# docker options DEFAULT: " -v $HOME:$HOME --user $UID "
DOCKER_OPTIONS=" -v $HOME:$HOME --user $UID "
# image that docker will use to run agent, image is pulled before it is run.
DOCKER_IMAGE=simexp/niak-boss:latest
# Worker init script
WORKER_INIT_SCRIPT=psom_worker.py


if [ $# != 2 ]
then
    echo "usage: $0 <psom_dir> <output_dir> <worker_id>"
    exit 1
fi

PSOM_DIR=$1s
OUTPUT_DIR=$2
WORKER_ID=$3


${DOCKER_EXE} pull ${DOCKER_IMAGE}
${DOCKER_EXE} run --rm ${DOCKER_OPTIONS} ${DOCKER_IMAGE} \
    /bin/bash -lic "${PSOM_DIR}/${WORKER_INIT_SCRIPT} -d ${OUTPUT_DIR} -w ${WORKER_ID}"