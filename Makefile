# NuJuliet Fortran port ... build with `make` (optimized) or `make debug` (runtime checks) ...
# Toolchain: COMPILER=intel (default; ifx + MKL, `module load intel-oneapi-hpc-compilers`)
#            COMPILER=gnu   (gfortran + reference LAPACK/BLAS)
COMPILER ?= intel
BUILD    ?= build
EXE      ?= nujuliet

ifeq ($(COMPILER),intel)
ifeq ($(origin FC),default)
FC = ifx
endif
# MKL shipped with the oneAPI HPC toolkit (located relative to ifx unless MKLROOT is set)
MKLROOT ?= $(abspath $(dir $(shell which $(FC)))/../../../mkl/latest)
MKLLIB   = $(MKLROOT)/lib
# rpath to the compiler runtime (libiomp5, libimf) too, so the binary runs without the module loaded
FCLIB    = $(abspath $(dir $(shell which $(FC)))/../lib)
MKLTHR  ?= intel_thread
LAPACK  ?= -Wl,--disable-new-dtags -Wl,-rpath,$(FCLIB) -L$(MKLLIB) -Wl,-rpath,$(MKLLIB) -Wl,--no-as-needed \
           -lmkl_intel_lp64 -lmkl_$(MKLTHR) -lmkl_core -lpthread -lm -ldl
FFLAGS  ?= -O2 -march=core-avx2 -qopenmp   # portable: Intel & AMD nodes (-xHost binaries refuse to run on AMD)
DBGFLAGS = -O0 -g -qopenmp -check all -traceback -warn all -warn nounused
COMMON   = -module $(BUILD) -I $(BUILD)
else
ifeq ($(origin FC),default)
FC = gfortran
endif
LAPACK  ?= -llapack -lblas
FFLAGS  ?= -O2 -fopenmp
DBGFLAGS = -O0 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wno-unused-dummy-argument -Wno-maybe-uninitialized
COMMON   = -ffree-line-length-none -J $(BUILD) -I $(BUILD)
endif

# Sources in dependency order ... chain fragments define SRC_RPA (chain A), SRC_QP (chain B) & SRC_EMPM (chain C) ...
include src/rpa/files.mk
include src/qp/files.mk
include src/empm/files.mk
# Stub dispatchers (src/*/stub) allow building one chain without the other, e.g.
#   make BUILD=build_qp EXE=nujuliet_qp SRC_RPA=src/rpa/stub/nj_dispatch_rpa.f90

SRC_CORE = src/core/nj_kinds.f90 \
           src/core/nj_utils.f90 \
           src/core/nj_wigner.f90 \
           src/core/nj_hashmap.f90 \
           src/core/nj_linalg.f90 \
           src/core/nj_types.f90 \
           src/core/nj_params.f90 \
           src/core/nj_functions.f90 \
           src/core/nj_orbitals.f90 \
           src/core/nj_o1b.f90 \
           src/core/nj_o2b.f90 \
           src/core/nj_v3b.f90 \
           src/core/nj_v2b.f90 \
           src/core/nj_obdm.f90 \
           src/core/nj_tr1b.f90

SRC_HF = src/hf/nj_hf.f90 \
         src/hf/nj_hf_mbpt.f90 \
         src/hf/nj_hf_efa.f90 \
         src/hf/nj_hf_driver.f90

SRC_NCSM = src/ncsm/nj_ncsm.f90

SRC = $(SRC_CORE) $(SRC_HF) src/app/nj_input.f90 $(SRC_RPA) $(SRC_QP) $(SRC_EMPM) $(SRC_NCSM) src/app/nujuliet.f90

OBJ = $(patsubst src/%.f90,$(BUILD)/%.o,$(SRC))

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@ $(LAPACK)

# Every object depends on all preceding ones (module order) ... simple & robust ...
define make_rule
$(BUILD)/$(1:src/%.f90=%.o): $(1) $(2)
	@mkdir -p $$(dir $$@)
	$$(FC) $$(FFLAGS) $$(COMMON) -c $$< -o $$@
endef
$(foreach i,$(shell seq 1 $(words $(SRC))),$(eval $(call make_rule,$(word $(i),$(SRC)),$(patsubst src/%.f90,$(BUILD)/%.o,$(wordlist 1,$(shell expr $(i) - 1),$(SRC))))))

debug:
	$(MAKE) BUILD=build_debug FFLAGS="$(DBGFLAGS)" nujuliet_debug

nujuliet_debug: $(patsubst src/%.f90,build_debug/%.o,$(SRC))
	$(FC) $(DBGFLAGS) $^ -o $@ $(LAPACK)

clean:
	rm -rf build build_debug nujuliet nujuliet_debug

.PHONY: debug clean
