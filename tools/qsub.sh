#!/bin/bash
set -x

runner=${20}
marker=${21}
params=${22}
out=${13}
err=${15}
cmdline="${runner} ${marker} \"${params}\""
echo "${params}"
eval "${cmdline}" &  #>${out} 2>${err}
#cat ${err}

set +x