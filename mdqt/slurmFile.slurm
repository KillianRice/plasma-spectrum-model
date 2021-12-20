#!/bin/bash
#SBATCH --job-name=sims
#SBATCH --array=0-9
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=1024M

echo "My job ran on:"
echo $SLURM_NODELIST
if [[ -d $SHARED_SCRATCH/$USER && -w $SHARED_SCRATCH/$USER ]]
then
    srun main -p 10.27.21 -a $SLURM_ARRAY_TASK_ID -s ucnp.settings
fi
