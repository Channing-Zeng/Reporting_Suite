#!/bin/bash
#set -x

date >&2
hostname >&2
echo >&2
echo "$2" >&2
echo >&2
echo >&2
bash -c "$2"
echo "$?">$1
echo >&2
date >&2

#set +x