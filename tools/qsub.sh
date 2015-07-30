#!/bin/bash
set -x

runner=${18}
marker=${19}
params=${20}
out=${11}
err=${13}
cmdline="${runner} ${marker} \"${params}\""
echo "${params}"
eval "${cmdline}" &  #>${out} 2>${err}
#cat ${err}

set +x