BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all ninja install clean
all: ninja

$(BUILD_NINJA):
	cmake --preset default

ninja: $(BUILD_NINJA)
	cmake --build --preset default

install: ninja
	cd $(BUILD_DIR) && ninja install

clean:
	rm -rf $(BUILD_DIR)
