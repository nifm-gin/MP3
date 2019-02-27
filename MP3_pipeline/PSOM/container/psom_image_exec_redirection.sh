#!/bin/bash
# This script is used to execute container command (typically qsub) on the host

echo ${PSOM_FIFO}

[[ -p ${PSOM_FIFO} ]] || echo No FIFO to redirect to ;

cmd=$@

echo "$cmd" > ${PSOM_FIFO};
