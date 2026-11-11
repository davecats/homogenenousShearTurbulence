# Environment for HoreKA (KIT).
# Usage:  source env/horeka.sh gpu    then   make GPU=1 GPU_ARCH=cc80 NCCL=1   (A100)
#                                                        GPU_ARCH=cc90          (H100)
#         source env/horeka.sh cpu    then   make
#
# The GPU build uses NVHPC's own HPC-X OpenMPI (CUDA-aware), not the system
# mpi/openmpi module, exactly as the channel code is built there.

module purge
case "${1:-gpu}" in
  gpu)
    module load toolkit/nvidia-hpc-sdk/25.3
    export PATH=$NVHPC_ROOT/comm_libs/mpi/bin:$PATH
    # NVHPC ships a relocated HPC-X OpenMPI.  Without its real prefix,
    # mpirun/srun can find the libraries but MPI_Init looks for the OpenMPI
    # help files under the uninstalled HPC-X build path.
    export OPAL_PREFIX=$NVHPC_ROOT/comm_libs/mpi
    export UCX_MEMTYPE_CACHE=n
    export GPU_ARCH=${GPU_ARCH:-cc80}
    ;;
  cpu)
    module load compiler/gnu/13 mpi/openmpi/5.0 numlib/fftw/3.3_serial
    export FFTW_DIR=/software/all/numlib/fftw/3.3_serial_gnu_13
    ;;
  *)
    echo "usage: source env/horeka.sh [gpu|cpu]"
    ;;
esac
