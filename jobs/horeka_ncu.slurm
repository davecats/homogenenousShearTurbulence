#!/bin/bash
# HoreKA Nsight Compute: hardware counters of a few launches of one kernel
# family on one A100 (DRAM bytes, achieved occupancy, registers, block size),
# the way to tell a bandwidth-bound kernel from a latency-bound one.
#   cd ~/hst && sbatch jobs/horeka_ncu.slurm
#   sbatch --export=ALL,KERNEL=regex:buildrhs,SKIP=6,COUNT=2 jobs/horeka_ncu.slurm
# KERNEL is an ncu -k filter (default: the three line-solver kernels), SKIP
# and COUNT the launches skipped and profiled among those that match.
#SBATCH --job-name=hst-ncu
#SBATCH --account=hk-project-exasim
#SBATCH --partition=dev_accelerated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --output=hst-ncu-%j.out

HST_ROOT=${HST_ROOT:-$HOME/hst}
RUNS=${RUNS:-$HOME/hst-runs}
DECK=${DECK:-bench_256}
KERNEL=${KERNEL:-regex:linsolve}
SKIP=${SKIP:-9}
COUNT=${COUNT:-6}
source $HST_ROOT/env/horeka.sh gpu
NCU=$(which ncu 2>/dev/null || ls $NVHPC_ROOT/compilers/bin/ncu $NVHPC_ROOT/profilers/*/Nsight_Compute/ncu 2>/dev/null | head -1)
nvidia-smi --query-gpu=name,memory.total,clocks.max.memory,clocks.max.sm --format=csv
d=$RUNS/ncu-$DECK; mkdir -p $d; cd $d; rm -f Dati.cart.out Runtimedata
sed -e "s/nstep = [0-9]*/nstep = 2/" $HST_ROOT/examples/$DECK.in > hst.in
echo "=== $DECK, kernels $KERNEL, launches $SKIP+1..$SKIP+$COUNT"
mpirun -np 1 --bind-to none $NCU -k "$KERNEL" --launch-skip $SKIP --launch-count $COUNT \
  --section SpeedOfLight --section Occupancy --section WarpStateStats --section LaunchStats \
  --metrics gpu__time_duration.sum,dram__bytes_read.sum,dram__bytes_write.sum,dram__throughput.avg.pct_of_peak_sustained_elapsed,lts__t_bytes.sum,l1tex__t_bytes.sum,launch__grid_size,launch__block_size,launch__registers_per_thread,launch__occupancy_limit_registers,sm__warps_active.avg.pct_of_peak_sustained_active,sm__throughput.avg.pct_of_peak_sustained_elapsed,smsp__inst_executed.sum,smsp__average_warp_latency_per_inst_issued.ratio \
  $HST_ROOT/build-gpu/hst hst.in 2>&1 | grep -v -E "^==PROF==|^ +[a-z]+ +=|^ *$|Authorization|hcoll"
