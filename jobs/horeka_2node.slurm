#!/bin/bash
# HoreKA two-node measurement: the same benchmarks on 4 A100 of one node and
# on 8 A100 of two nodes (InfiniBand between them), with the per-phase timer.
# Configurations: transport nccl / mpi, ranks bound to the NUMA domain of
# their GPU (jobs/bind_numa.sh) or not, npy y slabs (0 = one per node, 1 =
# x-z pencils only, 2 = one slab per node of two).  A run marked info prints
# NCCL's network setup (NCCL_DEBUG=INFO: adapter and GDRDMA per
# connection), and the regression decks run on 8 ranks with NPY_REG slabs
# against tests/reference.
#   cd ~/hst && sbatch --partition=accelerated jobs/horeka_2node.slurm     (HST_ROOT, RUNS, CONFIGS, NPY_REG overridable;
#                                                                            dev_accelerated has 3 nodes and rarely two free)
#   sbatch ... --export=ALL,HST_ROOT=$HOME/hst-y,NPY_REG=2,CONFIGS="4:bench_256:20:nccl:none:1 4:bench_256:20:nccl:none:2 8:bench_256:20:nccl:none:2:info" jobs/horeka_2node.slurm
#   sbatch ... --export=ALL,ROOTS="$HOME/hst-y $HOME/hst-exp",... jobs/horeka_2node.slurm     (A/B: the builds back to back)
#   sbatch --partition=accelerated --nodes=4 --time=00:15:00 --export=ALL,HST_ROOT=$HOME/hst-y,RUNS=$HOME/hst-runs/scal4,CONFIGS="16:bench_512:10:nccl:none:0 16:bench_1024:5:nccl:none:0:mem" jobs/horeka_2node.slurm
#   sbatch --partition=accelerated-h100 --export=ALL,HST_ROOT=$HOME/hst-y,BUILD=build-h100,... jobs/horeka_2node.slurm   (the H100 build, make GPU_ARCH=cc90 BUILD=build-h100)
# CONFIGS lists the runs as ranks:deck:steps:transport:bind:npy[:info|:lc<N>|:mem] (a deck
# of examples/ or of tests/decks/, e.g. small for the fixed cost of the y exchange;
# lc<N> sets line_chunk = N, the line solver's batch of x columns; mem samples
# nvidia-smi on the first node during the run and prints the peak memory per
# GPU, for the largest grid that fits); NCCL_* variables in the environment
# reach the runs (e.g. NCCL_NET_GDR_LEVEL=SYS lets NCCL use GPUDirect RDMA for
# GPUs that are not next to the adapter).  Four nodes: --nodes=4 on the sbatch
# line (npy = 0 chooses one y slab per node); two jobs of the
# same configs must not share RUNS.
#SBATCH --job-name=hst-2node
#SBATCH --account=hk-project-exasim
#SBATCH --partition=dev_accelerated
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --time=00:10:00
# The limit is the expected length plus a margin, not a round number: shorter
# jobs backfill sooner.  One root with the default CONFIGS takes 2-3.5 min, a
# two-root A/B of 7 configs 6.5 min: for that, sbatch --time=00:12:00.
#SBATCH --output=hst-2node-%j.out

HST_ROOT=${HST_ROOT:-$HOME/hst}
ROOTS=${ROOTS:-$HST_ROOT}         # the builds to run, back to back, in the same job
BUILD=${BUILD:-build-gpu}         # the build directory under each root (build-h100 on the H100 partition)
RUNS=${RUNS:-$HOME/hst-runs}
source $HST_ROOT/env/horeka.sh gpu
export NCCL_DEBUG=WARN
CONFIGS=${CONFIGS:-"4:bench_256:20:nccl:none:1 8:bench_256:20:nccl:none:1:info 8:bench_256:20:mpi:none:1 4:bench_512:10:nccl:none:1 8:bench_512:10:nccl:none:1"}
NPY_REG=${NPY_REG:-1}
env | grep '^NCCL_' | grep -v NCCL_DEBUG
scontrol show hostnames $SLURM_JOB_NODELIST | paste -sd' '
echo "=== topology of $(hostname)"; nvidia-smi topo -m | head -12; numactl -H | grep -E "available|cpus" | cut -c1-80
nvidia-smi --query-gpu=name,memory.total --format=csv | head -2
run() {  # run <ranks> <deck> <nstep> <transport> <bind: none|numa> <npy> [info|lc<N>]
  d=$RUNS/2node-$2-np$1-npy$6; mkdir -p $d; cd $d; rm -f Dati.cart.out statistics/Runtimedata
  deck=$HST_ROOT/examples/$2.in; [ -f $deck ] || deck=$HST_ROOT/tests/decks/$2.in     # a benchmark or a test deck
  lc=0; case "$7" in lc*) lc=${7#lc};; esac                                             # lc<N>: line_chunk = N
  sed -e "s/nstep = [0-9]*/nstep = $3/" -e "s/time_from_restart = .false./time_from_restart = .false., timing = .true./" \
      -e "s/&mesh /\&mesh transport = '$4', npy = $6, line_chunk = $lc, /" $deck > hst.in
  wrap=; [ "$5" = numa ] && wrap=$HST_ROOT/jobs/bind_numa.sh
  echo "=== $2 on $1 GPU(s), $3 steps, transport $4, bind $5, npy $6 $7"
  smi=; if [ "$7" = mem ]; then nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -l 1 > mem.log & smi=$!; fi
  if [ "$7" = info ]; then
    NCCL_DEBUG=INFO NCCL_DEBUG_SUBSYS=INIT,NET mpirun -np $1 --map-by ppr:4:node --bind-to none $wrap $HST_ROOT/$BUILD/hst hst.in > run.log 2>&1
    grep -E "NCCL version|NET/IB : Using|Using network|via NET|GDRDMA|NCCL INFO Channel 00/0" run.log | sed 's/.*NCCL INFO //' | sort | uniq -c | sort -rn | head -12
    grep -E -A12 "s/step; writing|rror|time per phase|ranks =|transport|bind_numa" run.log
  else
    mpirun -np $1 --map-by ppr:4:node --bind-to none $wrap $HST_ROOT/$BUILD/hst hst.in 2>&1 | grep -E -A12 "s/step; writing|rror|time per phase|ranks =|transport|bind_numa|MB per rank|out of memory|failed"
  fi
  if [ -n "$smi" ]; then kill $smi; awk -F', ' '$2 > m[$1] {m[$1] = $2} END {for (i in m) printf "   peak memory GPU %s: %d MB\n", i, m[i]}' mem.log | sort; fi
}
for HST_ROOT in $ROOTS; do
  echo "=== build $HST_ROOT"
  for c in $CONFIGS; do run ${c//:/ }; done
  echo "=== regression decks on 8 ranks, transport nccl, npy $NPY_REG"
  for deck in small small_s2x small_stokes; do
    d=$RUNS/2node-reg-$deck; mkdir -p $d; cd $d; rm -f Dati.cart.out
    sed -e "s/&mesh /\&mesh transport = 'nccl', npy = $NPY_REG, /" $HST_ROOT/tests/decks/$deck.in > hst.in
    mpirun -np 8 --map-by ppr:4:node --bind-to none $HST_ROOT/$BUILD/hst hst.in > run.log 2>&1 || echo "$deck: run failed"
    echo "--- $deck"; python3 $HST_ROOT/tests/compare_fields.py Dati.cart.out $HST_ROOT/tests/reference/$deck.fld 1e-10 | tail -1
  done
done
