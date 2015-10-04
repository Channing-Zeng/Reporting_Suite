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
eval "${cmdline}" | tee ${out}
#cat ${err}

#set +x