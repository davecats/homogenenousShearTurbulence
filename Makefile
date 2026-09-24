# hst -- DNS of homogeneous shear turbulence
#
#   make            CPU build: gfortran (mpifort) + FFTW
#   make GPU=1      GPU build: nvfortran (NVHPC mpif90) + cuFFT, OpenMP offload
#   make clean
#
# Source `env/istm.sh` or `env/horeka.sh` first so the right compilers are on
# PATH.  GPU_ARCH selects the device: cc86 (RTX 3060, A6000), cc80 (A100),
# cc90 (H100), cc120 (RTX 5090).  Override on the command line:
#   make GPU=1 GPU_ARCH=cc90
# BUILD=<dir> keeps builds for different machines apart on a shared home:
#   make GPU=1 GPU_ARCH=cc120 BUILD=build-corax

# Sources in dependency order (each file only uses the ones above it).
SRC = src/hst_params.f90 \
      src/hst_input.f90 \
      src/hst.f90

# The MPI wrapper is chosen here unless FC is given on the command line.
# (make has a built-in default FC=f77, and compiler modules export FC=nvfortran
# or gfortran; neither is what we want, so only `make FC=...` wins.)
ifneq ($(origin FC),command line)
  FC = $(if $(GPU),mpif90,mpifort)
endif
# Same for FFLAGS and LIBS: modules export FFLAGS too.
ifneq ($(origin FFLAGS),command line)
  undefine FFLAGS
endif
ifneq ($(origin LIBS),command line)
  undefine LIBS
endif

ifeq ($(GPU),1)
  GPU_ARCH ?= cc86
  FFLAGS   ?= -cpp -O3 -mp=gpu -gpu=$(GPU_ARCH) -cuda -Minfo=mp -DHAVE_CUDA
  LIBS     ?= -cudalib=cufft
  MODFLAG   = -module $(BUILD)
  BUILD    ?= build-gpu
else
  FFLAGS   ?= -cpp -O2 -g -ffree-line-length-none -fbacktrace -DHAVE_FFTW
  LIBS     ?= -lfftw3 -lm
  MODFLAG   = -J$(BUILD)
  BUILD    ?= build-cpu
endif

OBJ = $(patsubst src/%.f90,$(BUILD)/%.o,$(SRC))
EXE = $(BUILD)/hst

all: $(EXE)
	@echo "built $(EXE)"

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS)

$(BUILD)/%.o: src/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -I$(BUILD) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

# Module dependencies (so that `make -j` stays correct).
$(BUILD)/hst_input.o: $(BUILD)/hst_params.o
$(BUILD)/hst.o:       $(BUILD)/hst_params.o $(BUILD)/hst_input.o

clean:
	rm -rf build-cpu build-gpu

.PHONY: all clean
