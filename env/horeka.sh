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
