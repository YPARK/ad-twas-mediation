#!/bin/bash -l

command="$@"

source /broad/software/scripts/useuse > /dev/null
reuse -q R-3.3
reuse -q GCC-5.2
reuse -q .icc-2015

export MKROOT=/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/intel/mkl
export LD_LIBRARY_PATH=${LIBRARY_PATH}:${LD_LIBRARY_PATH}
export LD_PRELOAD=${MKLROOT}/lib/intel64/libmkl_core.so:${MKLROOT}/lib/intel64/libmkl_sequential.so

printf "[%s] Running ... \n$command" "$(date)" > /dev/stderr
printf "\n\n" > /dev/stderr

$command

printf "[%s] \nDone\n\n" "$(date)" > /dev/stderr
