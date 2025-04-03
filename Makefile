BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all ninja test install clean
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
	fpm build --flag "-I$(HDF5_INCLUDE) -I$(NETCDF_FORTRAN_INCLUDE) -I$(OPENBLAS_INCLUDE) -L$(OPENBLAS_LIB)"

clean:
	rm -rf $(BUILD_DIR)
	rm -rf tools/h5merge/build
