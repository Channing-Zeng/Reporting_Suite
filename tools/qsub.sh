#!/bin/bash
#set -x

runner=${18}
done_marker=${19}
err_marker=${20}
params=${21}
out=${11}
err=${13}
cmdline="${runner} ${done_marker} ${err_marker} \"${params}\""
echo "${params}"
echo "${out}"
eval "${cmdline}" 2>&1 | tee ${out}
status=$?

if [ "${status}" -ne 0 ]; then
    exit 1
fi

#cat ${err}

#set +x