# hst

Direct numerical simulation of homogeneous shear turbulence on CPUs and
NVIDIA GPUs.  Derived from the [`channel`](https://github.com/davecats/channel)
DNS code and from the CPL HST code; the design and the work plan are in
[PLAN.md](PLAN.md).

**Status: work package 0** (skeleton: build system, input, MPI and device
setup).  The solver itself is not in yet.

## What it solves

Incompressible Navier-Stokes fluctuations about the mean flow `U = S*y` in a
box that is periodic in `x` (streamwise) and `z` (spanwise) and
shear-periodic in `y`.  Fourier in `x` and `z`, compact 4th-order finite
differences in `y`, velocity-vorticity (`v`, `eta`) formulation, RK3 /
Crank-Nicolson time stepping with the mean-shear advection integrated
analytically.

Note the axis convention: the shear direction is `y`, as in the channel
code.  The CPL code calls that direction `z`.

## Build

```bash
source env/istm.sh            # istmio2 / istmcetus / istmcorax
make                          # CPU: gfortran + MPI + FFTW  -> build-cpu/hst
make GPU=1                    # GPU: nvfortran + cuFFT       -> build-gpu/hst
```

On HoreKA:

```bash
source env/horeka.sh gpu;  make GPU=1 GPU_ARCH=cc80     # A100 (cc90 for H100)
source env/horeka.sh cpu;  make
```

`GPU_ARCH` names the device (cc86 RTX 3060 / A6000, cc120 RTX 5090, cc80
A100, cc90 H100); `BUILD=<dir>` keeps builds for different machines apart
on a shared home.  No CMake, no fypp, no Python.

## Run

From a directory containing `hst.in` (a Fortran namelist, see the commented
example in this repository):

```bash
mpirun -np 4 /path/to/build-cpu/hst          # or build-gpu/hst, one rank per GPU
sbatch jobs/horeka_gpu.slurm                 # HoreKA, 1 node, 4 GPUs
```

The number of ranks must divide both `nx+1` and `nzd` (`3*nz` rounded up to
a power of two times at most one factor three).
