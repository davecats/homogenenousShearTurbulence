#!/bin/bash
# HoreKA benchmark: bench_64 on 1 GPU, bench_256 on 1 and 4 GPUs, bench_512 on 1 and 4 GPUs.
#   cd ~/hst && sbatch jobs/horeka_bench_all.slurm      (HST_ROOT, RUNS overridable)
#SBATCH --job-name=hst-bench
#SBATCH --account=hk-project-exasim
#SBATCH --partition=dev_accelerated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=hst-bench-%j.out

HST_ROOT=${HST_ROOT:-$HOME/hst}
RUNS=${RUNS:-$HOME/hst-runs}
source $HST_ROOT/env/horeka.sh gpu
nvidia-smi -L | head -1
run() {  # run <ranks> <deck> <nstep>
  d=$RUNS/bench-$2-np$1; mkdir -p $d; cd $d; rm -f Dati.cart.out Runtimedata
  sed -e "s/nstep = [0-9]*/nstep = $3/" $HST_ROOT/examples/$2.in > hst.in
  echo "=== $2 on $1 GPU(s), $3 steps"
  mpirun -np $1 --bind-to none $HST_ROOT/build-gpu/hst hst.in | grep -E "step 5:|s/step; writing|rror|nxB|line solver"
  nvidia-smi --query-gpu=memory.used --format=csv,noheader | head -1
}
run 1 bench_64 50
run 1 bench_256 20
run 4 bench_256 20
run 1 bench_512 10
run 4 bench_512 10
