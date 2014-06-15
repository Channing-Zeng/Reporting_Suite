#!/bin/bash
runner=${18}
params=${19}
out=${11}
err=${13}
eval "${runner} ${params}" > ${out} 2> ${err}