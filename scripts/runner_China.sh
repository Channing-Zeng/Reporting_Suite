#!/bin/bash
#set -x

DONE_MARKER_FILE=$1
ERROR_MARKER_FILE=$2
CMDLINE="${@:3}"

date >&2
hostname >&2
source /etc/profile.d/modules.sh >&2
module unload python >&2 2>&2
module unload gcc >&2 2>&2
module load gcc/4.7 python/2.7.5 java perl >&2 2>&2
source /ngs/RDI/SCRIPTS/load_postproc.sh
echo >&2
echo "${CMDLINE}" >&2
echo >&2
echo >&2
eval "${CMDLINE}"
status=$?
if [ "${status}" -ne 0 ]; then
    echo "ERROR: command exit with code ${status}" >&2
    echo "${status}">${ERROR_MARKER_FILE}
else
    echo "$?">${DONE_MARKER_FILE}
fi
echo >&2
date >&2

#set +x