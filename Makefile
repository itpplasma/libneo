BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all deps ninja test install clean
all: ninja tools/h5merge/build/h5merge.x

deps:
	apt-get update
	apt-get install -y --no-install-recommends \
		gfortran \
		cmake \
		ninja-build \
		pkg-config \
		openmpi-bin \
		libopenmpi-dev \
		libfftw3-dev \
		libhdf5-dev \
		libnetcdff-dev \
		libnetcdf-dev \
		libgsl-dev \
		libopenblas-dev \
		liblapack-dev
	rm -rf /var/lib/apt/lists/*
	python3 -m pip install --break-system-packages -r requirements.txt h5py

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
