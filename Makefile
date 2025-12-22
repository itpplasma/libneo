BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

# NVHPC configuration - auto-detect version
NVHPC_BASE := /opt/nvidia/hpc_sdk/Linux_x86_64
NVHPC_VERSION := $(shell ls -1 $(NVHPC_BASE) 2>/dev/null | grep -E '^[0-9]+\.[0-9]+$$' | sort -V | tail -1)
NVHPC_DIR := $(NVHPC_BASE)/$(NVHPC_VERSION)
NVHPC_HPCX := $(NVHPC_DIR)/comm_libs/13.0/hpcx/latest
NVHPC_MPIFORT := $(NVHPC_HPCX)/ompi/bin/mpifort
NVHPC_NVC := $(NVHPC_DIR)/compilers/bin/nvc
NVHPC_BUILD_DIR := build_nvfortran

.PHONY: all ninja test install clean nvfortran test-nvfortran clean-nvfortran
all: ninja tools/h5merge/build/h5merge.x

tools/h5merge/build/h5merge.x:
	if [ ! -d "tools/h5merge/build" ] ; then \
    	echo "Creating 'build' directory..."; \
    	mkdir -p tools/h5merge/build; \
	fi
	cd tools/h5merge/build && cmake .. && make

$(BUILD_NINJA):
	cmake --preset default

ninja: $(BUILD_NINJA)
	cmake --build --preset default

test: ninja
	cd $(BUILD_DIR) && ctest

install: ninja
	cd $(BUILD_DIR) && ninja install

fpm:
	fpm build

clean:
	rm -rf $(BUILD_DIR)
	rm -rf tools/h5merge/build

# NVHPC/nvfortran targets
# Requires: NVHPC SDK installed at /opt/nvidia/hpc_sdk, CUDA at /opt/cuda
nvfortran:
	@if [ -z "$(NVHPC_VERSION)" ]; then \
		echo "ERROR: NVHPC SDK not found at $(NVHPC_BASE)"; \
		exit 1; \
	fi
	@echo "Using NVHPC $(NVHPC_VERSION) from $(NVHPC_DIR)"
	@echo "HPC-X MPI: $(NVHPC_MPIFORT)"
	@if [ ! -f "$(NVHPC_MPIFORT)" ]; then \
		echo "ERROR: HPC-X mpifort not found at $(NVHPC_MPIFORT)"; \
		echo "Check that comm_libs/13.0/hpcx exists in your NVHPC installation"; \
		exit 1; \
	fi
	@# Source HPC-X and run cmake
	bash -c '\
		export NVHPC_CUDA_HOME=$${NVHPC_CUDA_HOME:-/opt/cuda}; \
		source $(NVHPC_HPCX)/hpcx-mt-init.sh && hpcx_load && \
		cmake -S . -B $(NVHPC_BUILD_DIR) -G Ninja \
			-DCMAKE_Fortran_COMPILER=$(NVHPC_MPIFORT) \
			-DCMAKE_C_COMPILER=$(NVHPC_NVC) \
			-DCMAKE_BUILD_TYPE=Release \
			-DFORCE_FETCH_DEPS=ON && \
		cmake --build $(NVHPC_BUILD_DIR) -j'

test-nvfortran: nvfortran
	cd $(NVHPC_BUILD_DIR) && ctest

clean-nvfortran:
	rm -rf $(NVHPC_BUILD_DIR)
