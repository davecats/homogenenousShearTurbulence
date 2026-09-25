# hst -- DNS of homogeneous shear turbulence
#
#   make            CPU build: gfortran (mpifort) + FFTW     -> build-cpu/hst
#   make GPU=1      GPU build: nvfortran (NVHPC mpif90) + cuFFT, OpenMP offload
#                                                            -> build-gpu/hst
#   make test       also build the test programs (tests/*.f90, on tests/test_common.f90) into the build dir
#   make clean
#
# Source `env/istm.sh` or `env/horeka.sh` first so the right compilers are on
# PATH.  GPU_ARCH selects the device: cc86 (RTX 3060, A6000), cc80 (A100),
# cc90 (H100), cc120 (RTX 5090).  Override on the command line:
#   make GPU=1 GPU_ARCH=cc90
# BUILD=<dir> keeps builds for different machines apart on a shared home:
#   make GPU=1 GPU_ARCH=cc120 BUILD=build-corax
# FFTW_DIR is where the CPU build finds include/fftw3.f03 and lib/libfftw3
# (default /usr; the env scripts set it where needed).

# Sources in dependency order (each file only uses the ones above it).
SRC = src/hst_params.f90 \
      src/hst_input.f90 \
      src/hst_mpi.f90 \
      src/hst_fft.f90 \
      src/hst_setup.f90 \
      src/hst_transforms.f90 \
      src/hst_initial.f90 \
      src/hst_derivatives.f90 \
      src/hst_io.f90 \
      src/hst_linsolve.f90 \
      src/hst_stokes.f90 \
      src/hst_equations.f90 \
      src/hst_stats.f90 \
      src/hst_pressure.f90

TESTS = tests/test_roundtrip.f90 tests/test_linsolve.f90 tests/test_kelvin.f90 tests/test_pressure.f90 tests/test_taylorgreen.f90 tests/test_forcing.f90 tests/test_conservation.f90 tests/test_stokes.f90

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
  FFTW_DIR ?= /usr
  FFLAGS   ?= -cpp -O2 -g -ffree-line-length-none -fbacktrace -DHAVE_FFTW -I$(FFTW_DIR)/include
  LIBS     ?= -L$(FFTW_DIR)/lib -lfftw3 -lm
  MODFLAG   = -J$(BUILD)
  BUILD    ?= build-cpu
endif

OBJ     = $(patsubst src/%.f90,$(BUILD)/%.o,$(SRC))
TESTOBJ = $(BUILD)/test_common.o
EXE     = $(BUILD)/hst
TESTEXE = $(patsubst tests/%.f90,$(BUILD)/%,$(TESTS))

all: $(EXE)
	@echo "built $(EXE)"

test: $(TESTEXE)
	@echo "built $(TESTEXE)"

$(EXE): $(OBJ) $(BUILD)/hst.o
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(BUILD)/hst.o $(LIBS)

$(BUILD)/%: $(OBJ) $(TESTOBJ) $(BUILD)/%.o
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(TESTOBJ) $(BUILD)/$*.o $(LIBS)

$(BUILD)/%.o: src/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -I$(BUILD) -c $< -o $@

$(BUILD)/%.o: tests/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -I$(BUILD) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

# Module dependencies (so that `make -j` stays correct).
$(BUILD)/hst_input.o:      $(BUILD)/hst_params.o
$(BUILD)/hst_mpi.o:        $(BUILD)/hst_params.o
$(BUILD)/hst_fft.o:        $(BUILD)/hst_params.o
$(BUILD)/hst_setup.o:      $(BUILD)/hst_params.o
$(BUILD)/hst_transforms.o: $(BUILD)/hst_params.o $(BUILD)/hst_mpi.o $(BUILD)/hst_fft.o
$(BUILD)/hst_initial.o:    $(BUILD)/hst_params.o
$(BUILD)/hst_io.o:         $(BUILD)/hst_params.o $(BUILD)/hst_mpi.o $(BUILD)/hst_initial.o $(BUILD)/hst_derivatives.o
$(BUILD)/hst_derivatives.o: $(BUILD)/hst_params.o
$(BUILD)/hst_linsolve.o:   $(BUILD)/hst_params.o $(BUILD)/hst_derivatives.o
$(BUILD)/hst_stokes.o:     $(BUILD)/hst_params.o
$(BUILD)/hst_equations.o:  $(BUILD)/hst_params.o $(BUILD)/hst_derivatives.o $(BUILD)/hst_linsolve.o $(BUILD)/hst_fft.o $(BUILD)/hst_transforms.o $(BUILD)/hst_stokes.o
$(BUILD)/hst_stats.o:      $(BUILD)/hst_params.o $(BUILD)/hst_linsolve.o $(BUILD)/hst_derivatives.o $(BUILD)/hst_stokes.o
$(BUILD)/hst_pressure.o:   $(BUILD)/hst_params.o $(BUILD)/hst_fft.o $(BUILD)/hst_transforms.o $(BUILD)/hst_linsolve.o $(BUILD)/hst_io.o $(BUILD)/hst_derivatives.o
$(BUILD)/hst.o:            $(OBJ)
$(TESTOBJ):                $(OBJ)
$(patsubst tests/%.f90,$(BUILD)/%.o,$(TESTS)): $(OBJ) $(TESTOBJ)

clean:
	rm -rf build-cpu build-gpu

.PHONY: all test clean
