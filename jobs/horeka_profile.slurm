#!/bin/bash
# HoreKA profile: bench_64 and bench_256 on one A100 with the per-phase timer,
# then nsys on five steps of bench_256 (kernel summary in the job output).
#   cd ~/hst && sbatch jobs/horeka_profile.slurm
#SBATCH --job-name=hst-prof
#SBATCH --account=hk-project-exasim
#SBATCH --partition=dev_accelerated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:20:00
#SBATCH --output=hst-prof-%j.out

HST_ROOT=${HST_ROOT:-$HOME/hst}
RUNS=${RUNS:-$HOME/hst-runs}
source $HST_ROOT/env/horeka.sh gpu
NSYS=$(which nsys 2>/dev/null || ls $NVHPC_ROOT/compilers/bin/nsys $NVHPC_ROOT/profilers/*/bin/nsys 2>/dev/null | head -1)
nvidia-smi -L | head -1
run() {  # run <deck> <nstep> [nsys]
  d=$RUNS/prof-$1; mkdir -p $d; cd $d; rm -f Dati.cart.out Runtimedata prof.*
  sed -e "s/nstep = [0-9]*/nstep = $2/" -e "s/time_from_restart = .false./time_from_restart = .false., timing = .true./" $HST_ROOT/examples/$1.in > hst.in
  echo "=== $1, $2 steps, 1 GPU $3"
  if [ "$3" = nsys ]; then
    mpirun -np 1 --bind-to none $NSYS profile -o prof -f true -t cuda --sample=none --cpuctxsw=none $HST_ROOT/build-gpu/hst hst.in > run.log 2>&1
    grep -E "end of run|rror" run.log
    $NSYS stats --force-export=true -q --report cuda_gpu_kern_sum --format csv prof.nsys-rep 2>/dev/null | python3 -c "
import sys, csv
rows = list(csv.DictReader(sys.stdin))
tot = sum(float(r['Total Time (ns)']) for r in rows); n = sum(int(r['Instances']) for r in rows)
print('kernel time total %.3f s, %d launches' % (tot/1e9, n))
for r in rows[:22]:
    print('%6.1f%% %8.3f ms  n=%5s  avg %8.1f us  %s' % (100*float(r['Total Time (ns)'])/tot, float(r['Total Time (ns)'])/1e6, r['Instances'], float(r['Avg (ns)'])/1e3, r['Name'][:90]))
"
  else
    mpirun -np 1 --bind-to none $HST_ROOT/build-gpu/hst hst.in | grep -E -A9 "end of run|rror|time per phase|line solver"
  fi
}
run bench_64 50
run bench_256 20
run bench_256 5 nsys
