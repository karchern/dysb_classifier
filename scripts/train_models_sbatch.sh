#!/bin/bash -e
#SBATCH -A zeller
#SBATCH --time=0-96:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
### The following number might need to be updated!
#SBATCH --array=1-4608
#SBATCH --error=/scratch/karcher/dysbiosis/err/%a.err
#SBATCH --output=/scratch/karcher/dysbiosis/out/%a.out
execLine=$(cat ${1} | head -n $SLURM_ARRAY_TASK_ID | tail -1)
echo $execLine
eval $execLine
