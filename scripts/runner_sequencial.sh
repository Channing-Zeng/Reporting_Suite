#!/bin/bash
#set -x

DONE_MARKER_FILE=$1
ERROR_MARKER_FILE=$2
CMDLINE="${@:3}"
PATH=$PATH:/Users/vlad/vagrant/NGS_Reporting/venv_ngs_reporting/bin:/Users/vlad/vagrant/miniconda2/bin:/Users/vlad/vagrant/miniconda2/bin:/Users/vlad/bin/:/Users/vlad/vagrant/reporting_suite/tools/:/Users/vlad/vagrant/reporting_suite/tools/bed_processing/:/Users/vladsaveliev/perl5/bin:/usr/local/bin/python:/usr/local/bin:/usr/local/Cellar/gcc/4.7.2/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/usr/local/mysql/bin:/usr/local/opt:/usr/local/sbin

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