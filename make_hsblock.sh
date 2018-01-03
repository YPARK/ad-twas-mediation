#!/bin/bash -l

source /broad/software/scripts/useuse > /dev/null
reuse -q GCC-5.2
reuse -q Python-2.7

[ -f $2/final.txt.gz ] || /home/unix/yjpark/work/common/bin/hsblock_cluster.py $1 -d 1 -o $2 -i 1000 -z 10 -r 20 -v -v -v
