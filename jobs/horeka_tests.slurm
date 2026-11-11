#!/bin/bash
# HoreKA GPU test job: runs the test suite on 4 GPUs, then the small deck on
# 4 GPUs, and compares its final field with a CPU run made on the login node
# (tests/decks/small.in, expected agreement ~1e-13).  NPY > 1 runs both with
# that many y slabs (NPY=2 is 2 x-z pencils x 2 slabs).
#   cd ~/hst && sbatch jobs/horeka_tests.slurm           (HST_ROOT, NPY overridable)
#SBATCH --job-name=hst-tests
#SBATCH --account=hk-project-exasim
#SBATCH --partition=dev_accelerated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --output=hst-tests-%j.out

HST_ROOT=${HST_ROOT:-$HOME/hst}
NPY=${NPY:-1}
source $HST_ROOT/env/horeka.sh gpu
cd $HST_ROOT
echo "=== test suite, build-gpu, 4 ranks, npy $NPY"
tests/run_tests.sh build-gpu 4 $NPY
echo "=== small deck, build-gpu, 4 ranks, npy $NPY"
mkdir -p $HOME/hst-runs/tests-gpu && cd $HOME/hst-runs/tests-gpu && rm -f Dati.cart.out statistics/Runtimedata
sed "s/&mesh /\&mesh npy = $NPY, /" $HST_ROOT/tests/decks/small.in > hst.in
mpirun -np 4 --bind-to none $HST_ROOT/build-gpu/hst hst.in | tail -2
if [ -f $HOME/hst-runs/tests-cpu/Dati.cart.out ]; then
  echo "=== GPU (4 ranks) vs CPU (login node) final field"
  python3 $HST_ROOT/tests/compare_fields.py Dati.cart.out $HOME/hst-runs/tests-cpu/Dati.cart.out
fi
