#!/bin/bash
#set -x

DONE_MARKER_FILE=$1
ERROR_MARKER_FILE=$2
CMDLINE="${@:3}"
PATH=/home/vsaveliev/NGS_Reporting/venv_ngs_reporting/bin:/home/vsaveliev/bcbio/anaconda/bin:/home/vsaveliev/bcbio_tools/bin:/home/vsaveliev/.linuxbrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/usr/lib/jvm/java-8-oracle/bin:/usr/lib/jvm/java-8-oracle/db/bin:/usr/lib/jvm/java-8-oracle/jre/bin:/home/vsaveliev/src/reporting_suite:/home/vsaveliev/src/reporting_suite/tools:/home/vsaveliev/src/reporting_suite/tools/bed_processing:$PATH

date >&2
hostname >&2
echo >&2
echo "${CMDLINE}" >&2
echo >&2
echo >&2
eval "${CMDLINE}"
status=$?
if [ "${status}" -ne 0 ]; then
    echo "${status}">${ERROR_MARKER_FILE}
    echo "Error: command returned code ${status}" >&2
    exit 1
else
    echo "${status}">${DONE_MARKER_FILE}
fi
echo >&2
date >&2

#set +x