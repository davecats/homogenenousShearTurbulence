# Environment for the ISTM machines (istmio2, istmcetus, istmcorax).
# Usage:  source env/istm.sh          then   make GPU=1   or   make
#
# CPU build needs nothing beyond the system mpifort and FFTW.
# GPU build uses the shared NVHPC tree; its HPC-X OpenMPI is CUDA-aware.

NVHPC_DIR=/opt/Nvidia/nvhpc/Linux_x86_64/25.9
export PATH=$NVHPC_DIR/compilers/bin:$NVHPC_DIR/comm_libs/mpi/bin:$PATH
export LD_LIBRARY_PATH=$NVHPC_DIR/compilers/lib:$NVHPC_DIR/math_libs/lib64:$NVHPC_DIR/cuda/lib64:$LD_LIBRARY_PATH

# Device architecture for `make GPU=1`
case "$(hostname)" in
  ISTM-corax) export GPU_ARCH=cc120 ;;   # RTX 5090
  *)          export GPU_ARCH=cc86  ;;   # RTX 3060 (io2), RTX A6000 (cetus)
esac

# Multi-rank GPU runs on these boxes need this (UCX bug, see channel/HPC_SESSION.md)
export UCX_MEMTYPE_CACHE=n
