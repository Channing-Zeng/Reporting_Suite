#!/bin/bash
#set -x

runner=${20}
done_marker=${21}
err_marker=${22}
params=${23}
out=${13}
err=${15}
cmdline="${runner} ${done_marker} ${err_marker} \"${params}\""
echo "${params}"
echo "${out}"
eval "${cmdline}" 2>&1 | tee ${out}
status=$?

if [ "${status}" -ne 0 ]; then
    exit 1
fi

#set +x