#!/bin/bash
runner=${18}
params=${19}
out=${11}
err=${13}
echo 'Submit job'
cmdline="${runner} \"${params}\" "
echo ${cmdline}
eval "${cmdline}"
#>${out} 2>${err}