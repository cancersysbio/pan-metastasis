#!/bin/bash -l
#
# set the name of the job
#SBATCH -J met-clonality
#
# set the maximum memory usage (per slot)
#SBATCH --mem=50000

# set the maximum run time
#SBATCH -t 20:00:00

# specify an email address
#SBATCH --mail-user="email"
#SBATCH --mail-type=ALL
#SBATCH --account="account"
#SBATCH -o met-clonality.out
#SBATCH -e met-clonality.error
#
#SBATCH --array 1-100


WORKDIR='dir'

cd $WORKDIR

##monclonal seeding
./clonality_simulation.py $SLURM_ARRAY_TASK_ID 1 > metastatic_clonality_${SLURM_ARRAY_TASK_ID}.out

##polyclonal seeding
#./clonality_simulation.py $SLURM_ARRAY_TASK_ID 10 > metastatic_clonality_${SLURM_ARRAY_TASK_ID}.out
