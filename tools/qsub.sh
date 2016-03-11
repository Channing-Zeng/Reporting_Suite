#!/bin/bash
#set -x

out=${13}
err=${15}
runner=${20}
done_marker=${21}
err_marker=${22}
final_cmdline=${@:23}
runner_cmdline="${runner} ${done_marker} ${err_marker} ${final_cmdline}"
echo "${final_cmdline}"
echo "${out}"
eval "${runner_cmdline}" 2>&1 | tee ${out}
status=$?

if [ "${status}" -ne 0 ]; then
    exit 1
fi

#set +x