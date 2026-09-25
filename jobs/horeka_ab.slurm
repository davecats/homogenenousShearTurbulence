#!/bin/bash
# HoreKA A/B: the same benchmarks on two builds back to back, ~/hst (A) and
# ~/hst-exp (B, an rsync copy of the variant under test), on the same node.
#   cd ~/hst && sbatch jobs/horeka_ab.slurm             (A_ROOT, B_ROOT, RUNS overridable)
# Per build: bench_256 and bench_512 on 1 and 4 A100 with the per-phase timer
# (4 GPUs: NCCL), then the nsys kernel summary of bench_256 on 1 A100.
#SBATCH --job-name=hst-ab
#SBATCH --account=hk-project-exasim
#SBATCH --partition=dev_accelerated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=hst-ab-%j.out

A_ROOT=${A_ROOT:-$HOME/hst}
B_ROOT=${B_ROOT:-$HOME/hst-exp}
RUNS=${RUNS:-$HOME/hst-runs}
source $A_ROOT/env/horeka.sh gpu
NSYS=$(which nsys 2>/dev/null || ls $NVHPC_ROOT/compilers/bin/nsys $NVHPC_ROOT/profilers/*/bin/nsys 2>/dev/null | head -1)
nvidia-smi -L | head -1
run() {  # run <root> <ranks> <deck> <nstep> [nsys]
  d=$RUNS/ab-$3-np$2; mkdir -p $d; cd $d; rm -f Dati.cart.out Runtimedata prof.*
  sed -e "s/nstep = [0-9]*/nstep = $4/" -e "s/time_from_restart = .false./time_from_restart = .false., timing = .true./" $1/examples/$3.in > hst.in
  echo "=== $(basename $1): $3 on $2 GPU(s), $4 steps $5"
  if [ "$5" = nsys ]; then
    mpirun -np $2 --bind-to none $NSYS profile -o prof -f true -t cuda --sample=none --cpuctxsw=none $1/build-gpu/hst hst.in > run.log 2>&1
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
    mpirun -np $2 --bind-to none $1/build-gpu/hst hst.in | grep -E -A11 "s/step; writing|rror|time per phase"
  fi
}
for root in $A_ROOT $B_ROOT; do
  run $root 1 bench_256 20
  run $root 4 bench_256 20
  run $root 1 bench_512 10
  run $root 4 bench_512 10
  run $root 1 bench_256 5 nsys
done
