#!/bin/bash

#./etc/sysconfig/pssc

#PBS -S /bin/bash
#PBS JOB_NAME="QSH_$(whoami)"
#PBS NODE_NUM="1"
#PBS NODE_PPN="${NODE_NCPUS}"
#PBS HOURS="24"
#PBS MINUTES="00"
#PBS SECONDS="00"
#PBS WALLTIME=${HOURS}:${MINUTES}:${SECONDS}
#PBS RES_LIST="nodes=${NODE_NUM}:ppn=${NODE_PPN}"
#PBS DIR_WORK="${PBS_O_WORKDIR}"
#PBS QUEUE="high"

#PBS cd ${DIR_WORK}

