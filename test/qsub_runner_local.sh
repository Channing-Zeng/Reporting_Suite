#!/bin/bash
set -x

DONE_MARKER_FILE=$1
ERROR_MARKER_FILE=$2
CMDLINE=$3

date
hostname
echo
echo "${CMDLINE}"
echo
echo
eval "${CMDLINE}"
status=$?
if [ "${status}" -ne 0 ]; then
    echo "${status}">${ERROR_MARKER_FILE}
    echo "Error: command returned code ${status}" >&2
    exit 1
else
    echo "${status}">${DONE_MARKER_FILE}
fi
echo
date

set +x