#!/bin/bash
# HoreKA GPU job: one node, one MPI rank per GPU.
#
#   cd <case directory with hst.in>
#   sbatch ~/hst/jobs/horeka_gpu.slurm            (A100 node, 4 GPUs)
#   sbatch -p accelerated-h100 ~/hst/jobs/horeka_gpu.slurm
#
# For a test use  -p dev_accelerated -t 00:10:00 .  Set HST_ROOT if the
# repository is not in ~/hst.  The binary must have been built with the
# matching GPU_ARCH (cc80 A100, cc90 H100).
#
#SBATCH --job-name=hst
#SBATCH --account=hk-project-exasim
#SBATCH --partition=accelerated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=hst-%j.out

HST_ROOT=${HST_ROOT:-$HOME/hst}
source $HST_ROOT/env/horeka.sh gpu

echo "job $SLURM_JOB_ID on $(hostname), $SLURM_NTASKS ranks, deck: $(pwd)/hst.in"
nvidia-smi -L

mpirun -np $SLURM_NTASKS --bind-to none $HST_ROOT/build-gpu/hst hst.in
