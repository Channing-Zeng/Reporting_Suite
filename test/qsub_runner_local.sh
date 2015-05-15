#!/bin/bash
set -x
date >&2
eval $@
date >&2
set +x