#!/bin/bash -e
#SBATCH -A zeller
#SBATCH --time=0-96:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
### The following number might need to be updated!
#SBATCH --array=1009,1020,1042,107,118,129,184,195,19,206,217,228,239,250,261,30,316,327,338,349,360,371,382,393,41,448,459,470,481,514,52,536,547,558,580,591,602,613,624,635,63,646,657,712,723,734,745,74,756,767,778,789,844,855,85,866,877,888,899,8,910,921,96,976,987,998
#SBATCH --error=/scratch/karcher/dysbiosis/err/%a.err
#SBATCH --output=/scratch/karcher/dysbiosis/out/%a.out

execLine=$(cat ${1} | head -n $SLURM_ARRAY_TASK_ID | tail -1)
eval $execLine
