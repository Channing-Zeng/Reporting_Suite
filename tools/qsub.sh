#!/bin/bash
set -x
runner=${18}
params=${19}
out=${11}
err=${13}
cmdline="${runner} \"${params}\" "
echo "${params}"
eval "${cmdline}"
#>${out} 2>${err}
set +x