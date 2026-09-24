#!/bin/bash
# HoreKA CPU job: one node, one MPI rank per core (the rank count must
# divide nx+1 and nzd; 64 fits the default deck: nx+1 = 64, nzd = 192).
#   cd <case directory with hst.in> && sbatch ~/hst/jobs/horeka_cpu.slurm
#SBATCH --job-name=hst-cpu
#SBATCH --account=hk-project-exasim
#SBATCH --partition=cpuonly
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=hst-%j.out

HST_ROOT=${HST_ROOT:-$HOME/hst}
source $HST_ROOT/env/horeka.sh cpu
echo "job $SLURM_JOB_ID on $(hostname), $SLURM_NTASKS ranks, deck: $(pwd)/hst.in"
mpirun -np $SLURM_NTASKS $HST_ROOT/build-cpu/hst hst.in
