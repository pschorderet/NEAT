#!/bin/bash
FIRST=`qsub -N Iterate -o /data/schorderet/projects/EXAMPLE/scripts/qsub -e /data/schorderet/projects/EXAMPLE/scripts/qsub /data/schorderet/projects/EXAMPLE/scripts/ChIPseq.sh`
