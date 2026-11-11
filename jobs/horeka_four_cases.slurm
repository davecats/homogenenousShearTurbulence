#!/bin/bash
# Run four independent single-GPU simulations concurrently on one HoreKA node.
# Submit from a directory containing sim1/, ..., sim4/, each with hst.in.
#SBATCH --job-name=hst4
#SBATCH --account=hk-project-exasim
#SBATCH --partition=dev_accelerated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:10:00
#SBATCH --output=hst-%j.out

set -euo pipefail

# Slurm executes a submitted batch script from a spool copy under
# /var/spool/slurmd, so do not derive HST_ROOT from BASH_SOURCE.  Set HST_ROOT
# to the actual checkout explicitly when submitting the job.
HST_ROOT=${HST_ROOT:?Set HST_ROOT to the repository root containing env/ and build-gpu/}
HST_EXE=${HST_EXE:-$HST_ROOT/build-gpu/hst}
CASE_ROOT=${CASE_ROOT:-$PWD}

if [[ ! -r "$HST_ROOT/env/horeka.sh" ]]; then
    echo "ERROR: cannot find $HST_ROOT/env/horeka.sh" >&2
    echo "       Set HST_ROOT to the repository root containing env/ and build-gpu/." >&2
    exit 2
fi
source "$HST_ROOT/env/horeka.sh" gpu

if [[ ! -f "$HST_EXE" ]]; then
    echo "ERROR: executable does not exist: $HST_EXE" >&2
    echo "       Set HST_EXE to the GPU build executable." >&2
    exit 2
fi
if [[ ! -x "$HST_EXE" ]]; then
    echo "ERROR: executable is not executable: $HST_EXE" >&2
    echo "       On HoreKA, run: chmod u+x $HST_EXE" >&2
    ls -l "$HST_EXE" >&2 || true
    exit 2
fi

for i in 1 2 3 4; do
    case_dir="$CASE_ROOT/sim$i"
    if [[ ! -d "$case_dir" || ! -r "$case_dir/hst.in" ]]; then
        echo "ERROR: missing simulation input: $case_dir/hst.in" >&2
        exit 2
    fi
done

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "MPI launcher: $(command -v mpirun)"
echo "OPAL_PREFIX=${OPAL_PREFIX:-unset}"
echo "HST_ROOT=$HST_ROOT"
echo "HST_EXE=$HST_EXE"
nvidia-smi -L

run_sim() {
    local dir=$1
    echo "Starting $dir"
    # The job allocation contains four Slurm tasks.  This single-task srun
    # step gets one GPU and one CPU allocation; --exclusive prevents overlap.
    # hst is MPI-enabled, so launch its one rank through the matching OpenMPI
    # rather than executing the MPI binary directly under srun.  --mpi=none
    # prevents Slurm from trying to provide a different MPI wireup layer.
    srun --exclusive --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
        --gpus-per-task=1 --gpu-bind=closest --mpi=none \
        --chdir="$dir" \
        mpirun -np 1 --bind-to none "$HST_EXE" hst.in >"$dir/hst.out" 2>&1
    local rc=$?
    if (( rc == 0 )); then
        echo "$dir finished successfully"
    else
        echo "$dir FAILED with exit code $rc"
    fi
    return "$rc"
}

pids=()
for i in 1 2 3 4; do
    dir="$CASE_ROOT/sim$i"
    run_sim "$dir" &
    pids+=("$!")
done

status=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        status=1
    fi
done
exit "$status"