#!/bin/bash
#SBATCH --job-name="dmapp_put"
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --time=0:30:00

export XT_SYMMETRIC_HEAP_SIZE=2000M

cd ${SLURM_SUBMIT_DIR}

#aprun -n 2 -N 2 ${JOB_DIR}/dmapp_put.x  1> ${SLURM_SUBMIT_DIR}/dmapp_put.${SLURM_JOBID}.out 2> ${SLURM_SUBMIT_DIR}/dmapp_put.${SLURM_JOBID}.err
aprun -n 2 -N 2 dmapp_put.x
