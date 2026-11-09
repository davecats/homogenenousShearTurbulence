#!/usr/bin/env bash

set -euo pipefail

make clean

##############################
# COMPILE AND TEST FOR CPU
##############################
BUILD=${BUILD:-build-cpu}
RANKS=${RANKS:-4}
make BUILD="$BUILD"				    # CPU: gfortran + MPI + FFTW  -> build-cpu/hst
make test BUILD="$BUILD"
tests/run_tests.sh "$BUILD" "$RANKS"  # second argument: ranks

##############################
# COMPILE AND TEST FOR GPU
# (single-GPU nodes)
##############################
# source env/istm.sh          	        # istmio2 / istmcetus / istmcorax
# BUILD=${BUILD:-build-gpu-istmvolans}
# RANKS=${RANKS:-1}
# make GPU=1 BUILD="$BUILD" 		    # GPU: nvfortran + cuFFT -> $BUILD/hst
# make test GPU=1 BUILD="$BUILD"		# the test programs, same build directory
# tests/run_tests.sh "$BUILD" "$RANKS"  # second argument: ranks

##############################
# COMPILE AND TEST FOR GPU
# (multi-GPU nodes)
##############################
#source env/istm.sh           	        # istmio2 / istmcetus / istmcorax
#BUILD=${BUILD:-build-gpu-2}
#RANKS=${RANKS:-2}
#make GPU=1 NCCL=1 BUILD="$BUILD"        # the same with NCCL for the alltoall (multi-GPU nodes)
#make test GPU=1 NCCL=1 BUILD="$BUILD"   # the test programs, same build directory; NCCL=1 keeps the test link of hst_mpi.o (-DHAVE_NCCL) consistent
#tests/run_tests.sh "$BUILD" "$RANKS"    # second argument: ranks
