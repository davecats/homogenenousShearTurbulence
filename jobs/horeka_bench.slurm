#!/bin/bash
# HoreKA benchmark: examples/bench_256.in on 1 and 4 GPUs, bench_512.in on 4 GPUs.
#   cd ~/hst && sbatch jobs/horeka_bench.slurm
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
source $HST_ROOT/env/horeka.sh gpu
nvidia-smi -L | head -1
run() {  # run <ranks> <deck>
  d=$HOME/hst-runs/bench-$2-np$1; mkdir -p $d; cd $d; rm -f Dati.cart.out Runtimedata
  cp $HST_ROOT/examples/$2.in hst.in
  echo "=== $2 on $1 GPU(s)"
  mpirun -np $1 --bind-to none $HST_ROOT/build-gpu/hst hst.in | grep -E "step|end of run|rror|nxB"
}
run 1 bench_256
run 4 bench_256
run 4 bench_512
