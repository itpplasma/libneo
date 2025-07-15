BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all ninja test install clean coverage coverage-build coverage-report
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

# Coverage targets
coverage: coverage-build coverage-report

coverage-build:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake .. -DENABLE_COVERAGE=ON
	cd $(BUILD_DIR) && ninja
	cd $(BUILD_DIR) && ctest --output-on-failure

coverage-report:
	gcovr --root . --exclude 'build/*' --exclude 'doc/*' --exclude 'example/*' \
	      --exclude 'test/*' --exclude 'tools/*' --exclude 'src/contrib/*' \
	      --exclude 'matlab/*' --exclude 'python/*' --exclude 'extra/*' \
	      --html --html-details -o coverage.html --print-summary
