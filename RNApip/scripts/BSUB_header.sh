#!/bin/bash
#
#BSUB -a poe                  # set parallel operating environment
#BSUB -P project_code         # project code
#BSUB -J hybrid_job_name      # job name
#BSUB -W 24:00                # wall-clock time (hrs:mins)
#BSUB -n 32                   # number of tasks in job
#BSUB -R "span[ptile=4]"      # run four MPI tasks per node
#BSUB -q regular              # queue
#BSUB -e errors.%J.hybrid     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.hybrid     # output file name in which %J is replaced by the job ID
