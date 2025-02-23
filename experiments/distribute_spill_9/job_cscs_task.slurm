#!/bin/bash -l
#SBATCH --job-name="distribute_spill_9"
#SBATCH --account="ich030"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zahrakh@kth.se
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=31
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --output=/scratch/snx3000/bp000380/Plasticity.data/distribute_spill_9.o
#SBATCH --error=/scratch/snx3000/bp000380/Plasticity.data/distribute_spill_9.e
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1
# set up environment
source ~/snudda_env/bin/activate
module swap PrgEnv-cray PrgEnv-gnu
module load daint-mc
rm -rf x86_64
nrnivmodl  /users/bp000380/projects/plasticityfinal/mechanisms >/dev/null

srun -n 31 python paral_task.py