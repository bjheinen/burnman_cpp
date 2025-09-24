# Copyright (c) 2025 Benedict Heinen
#
# This file is part of burnman_cpp and is licensed under the
# GNU General Public License v3.0 or later. See the LICENSE file
# or <https://www.gnu.org/licenses/> for details.
#
# burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
#
# ------------------- END OF LICENSE SECTION -----------------
# Target Executable (to make)
TARGET_EXEC := burnman
TEST_EXEC := run_tests

# Compiler
CXX := g++

# Directories
CWD := $(shell pwd)
SRC_DIR := $(CWD)/src
INCLUDE_DIR := $(CWD)/include
BUILD_DIR := $(CWD)/build
BIN_DIR := $(CWD)/bin
TEST_DIR := $(CWD)/tests
INCLUDE_DIR_TEST := $(TEST_DIR)/include

# For GSL, etc. (Change if non-standard)
# Only defined here so they can be changed
# both should be on the system path anyway
EXTRA_INCLUDE := /usr/local/include /usr/include/eigen3
EXTRA_LIB := /usr/local/lib
LDFLAGS_COMMON := -lgsl -lgslcblas -lm
LDFLAGS_TEST := -lCatch2Main -lCatch2

# Expand to paths
INCLUDE_FLAGS_COMMON := -I$(INCLUDE_DIR) $(addprefix -I, $(EXTRA_INCLUDE))
INCLUDE_FLAGS_TEST := -I$(INCLUDE_DIR_TEST)
LIB_PATHS := $(addprefix -L, $(EXTRA_LIB))

# Compiler flags
CXXFLAGS := -Wall -Wextra -pedantic -Wshadow -Wconversion -std=c++17
CPPFLAGS := -MMD -MP $(INCLUDE_FLAGS_COMMON)

# Build type specific flags
DEBUG_FLAGS = -g -Og
RELEASE_FLAGS = -O3 -DNDEBUG
TEST_FLAGS = -g -O1 -fno-omit-frame-pointer -fsanitize=address,undefined

# Select Mode (default: release)
BUILD_MODE ?= release
ifeq ($(BUILD_MODE), debug)
  CXXFLAGS += $(DEBUG_FLAGS)
  BIN_SUFFIX = _debug
else ifeq ($(BUILD_MODE), test)
  CXXFLAGS += $(TEST_FLAGS)
  LDFLAGS_COMMON += -fsanitize=address,undefined
  BIN_SUFFIX = _test
else
  CXXFLAGS += $(RELEASE_FLAGS)
  BIN_SUFFIX =
endif
# Optional profiling (default: false)
PROFILE ?= 0
ifeq ($(PROFILE), 1)
  CXXFLAGS += -pg
  LDFLAGS_COMMON += -pg
	BIN_SUFFIX := $(BIN_SUFFIX)_profile
endif

# Store base build dir for clean
BASE_BUILD_DIR := $(BUILD_DIR)
BUILD_DIR := $(BUILD_DIR)/$(BUILD_MODE)
# Append -profile for profiling builds
ifeq ($(PROFILE), 1)
	BUILD_DIR := $(BUILD_DIR)-profile
endif

# Find source files
SOURCES := $(shell find $(SRC_DIR) -name '*.cpp')
TEST_SOURCES := $(shell find $(TEST_DIR) -name '*.cpp')

# Generate object fnames
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))
TEST_OBJECTS := $(patsubst $(TEST_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(TEST_SOURCES))

# Generate dependencies
DEPENDS := $(OBJECTS:.o=.d) $(TEST_OBJECTS:.o=.d)

# Output target
TARGET := $(BIN_DIR)/$(TARGET_EXEC)$(BIN_SUFFIX)
TEST_TARGET := $(BIN_DIR)/$(TEST_EXEC)$(BIN_SUFFIX)

# Default rule
all: $(TARGET)

# Link and build
$(TARGET) : $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(LIB_PATHS) $^ -o $@ $(LDFLAGS_COMMON)

# Link and build test executable
$(TEST_TARGET): $(OBJECTS) $(TEST_OBJECTS) | $(BIN_DIR)
	$(CXX) $(LIB_PATHS) $^ -o $@ $(LDFLAGS_COMMON) $(LDFLAGS_TEST)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Compile tests
$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp | $(BUILD_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(INCLUDE_FLAGS_TEST) $(CXXFLAGS) -c $< -o $@

# Create build dirs
$(BUILD_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

-include $(DEPENDS)

# Run tests
test: $(TEST_TARGET)
	@echo "Build complete. Run tests with: $(TEST_TARGET)"

# Clean up only a specific build-mode/profile dir/exec
clean-mode:
	rm -rf $(BUILD_DIR) $(TARGET) $(TEST_TARGET)

# Clean up everything
clean:
	rm -rf $(BASE_BUILD_DIR) $(BIN_DIR)

.PHONY: all clean clean-mode test

#NOTE: careful to test/checkout eigen opt flags, e.g. -DEIGEN_USE_THREADS -DEIGEN_DONT_ALIGN
