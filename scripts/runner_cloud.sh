#!/bin/bash
#set -x

date >&2
hostname >&2
echo >&2
echo "$2" >&2
echo >&2
echo >&2
eval $2
echo >&2
date >&2
touch $1 >&2

#set +x